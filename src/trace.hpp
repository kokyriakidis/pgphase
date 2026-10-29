// Stage trace (diagnostic, off by default): PGPHASE_TRACE=DIR writes DIR/trace.tsv, one line per decision that
// drops, merges, re-keys or re-genotypes a candidate, a site or an output record, with a reason code, so a variant can
// be followed through the pipeline (pgphase-dev scripts/dev/harness). Columns: contig, chunk (RegionChunk::chunk_id,
// -1 outside a chunk), pass (the union EM's solve, 1 or 2; 0 elsewhere), stage, event, beg, end (the reference span,
// 1-based inclusive, as the record would be written), key (TYPE:REF_LEN:ALT, or the VCF REF>ALT where the caller has
// it), reason, detail (space-separated key=value tokens and bare flags; docs/TRACE.md).

#pragma once

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <mutex>
#include <string>
#include <string_view>
#include <utility>

#include <htslib/hts.h>

#include "phasing_types.hpp"

namespace pgphase_trace {

inline FILE* sink() {
    static FILE* const file = [] {
        const char* dir = std::getenv("PGPHASE_TRACE");
        if (dir == nullptr || *dir == '\0') return static_cast<FILE*>(nullptr);
        FILE* fh = std::fopen((std::string(dir) + "/trace.tsv").c_str(), "w");
        if (fh != nullptr) std::fputs("contig\tchunk\tpass\tstage\tevent\tbeg\tend\tkey\treason\tdetail\n", fh);
        return fh;
    }();
    return file;
}

// Where this thread is working, written with every line: the contig, the chunk and the union EM pass.
struct Context {
    std::string contig = "-";
    int chunk = -1;
    int pass = 0;
};

inline Context& context() {
    thread_local Context c;
    return c;
}

// Sets the context for a scope (a chunk's worker, a batch's contig, one EM solve) and restores it after.
class Scope {
public:
    Scope(const std::string& contig, int chunk) : saved_(context()) {
        context().contig = contig.empty() ? "-" : contig;
        context().chunk = chunk;
    }
    explicit Scope(int chunk) : saved_(context()) { context().chunk = chunk; }
    static Scope pass(int pass) {
        Scope s(context().chunk);
        context().pass = pass;
        return s;
    }
    Scope(const Scope&) = delete;
    Scope& operator=(const Scope&) = delete;
    Scope(Scope&& o) noexcept : saved_(o.saved_), active_(o.active_) { o.active_ = false; }
    ~Scope() {
        if (active_) context() = saved_;
    }

private:
    Context saved_;
    bool active_ = true;
};

// Whether tracing is on: callers building a detail string test this first.
inline bool on() { return sink() != nullptr; }

// DIR/reads.tsv: one row per read as the phased BAM writes it -- its haplotype and phase set, the chunks that saw
// it, and per chunk its label there and how it came about (PhasingChunk::ReadTrace).
inline FILE* reads_sink() {
    static FILE* const file = [] {
        if (!on()) return static_cast<FILE*>(nullptr);
        FILE* fh = std::fopen((std::string(std::getenv("PGPHASE_TRACE")) + "/reads.tsv").c_str(), "w");
        if (fh != nullptr) std::fputs("contig\tbeg\tend\tread\thap\tps\tcopies\tchunks\n", fh);
        return fh;
    }();
    return file;
}

// A read's row: beg and end are the span of the sites it was called at (1-based, inclusive).
inline void emit_read(hts_pos_t beg, hts_pos_t end, std::string_view name, int hap, hts_pos_t ps, int copies,
                      std::string_view chunks) {
    FILE* fh = reads_sink();
    if (fh == nullptr) return;
    static std::mutex lock;
    const std::lock_guard<std::mutex> guard(lock);
    const std::string& contig = context().contig;
    std::fprintf(fh, "%s\t%lld\t%lld\t%.*s\t%d\t%lld\t%d\t%.*s\n", contig.c_str(), static_cast<long long>(beg),
                 static_cast<long long>(end), static_cast<int>(name.size()), name.data(), hap,
                 static_cast<long long>(ps), copies, static_cast<int>(chunks.size()), chunks.data());
}

// A chunk's per-read record, grown on demand; nullptr when tracing is off.
inline pgphase_collect::PhasingChunk::ReadTrace* read_trace(pgphase_collect::PhasingChunk& chunk, size_t ri) {
    if (!on()) return nullptr;
    if (chunk.read_trace.size() <= ri) chunk.read_trace.resize(std::max(chunk.reads.size(), ri + 1));
    return &chunk.read_trace[ri];
}

inline std::string describe_read(const pgphase_collect::PhasingChunk::ReadTrace& t) {
    std::string out = t.label.empty() ? "-" : t.label;
    if (!t.fill.empty()) out += "; fill " + t.fill;
    if (t.alignment_only) out += "; alignment_only filled=" + std::to_string(t.filled);
    if (t.realigned != 0 || t.graded != 0) out += "; realigned=" + std::to_string(t.realigned) +
                                                  " graded=" + std::to_string(t.graded);
    if (t.snp_replaced != 0) out += "; snp_replaced=" + std::to_string(t.snp_replaced);
    if (t.aln_hap >= 0) out += "; aln hap=" + std::to_string(t.aln_hap) + " ps=" + std::to_string(t.aln_ps);
    return out;
}

inline void emit(std::string_view stage, std::string_view event, hts_pos_t beg, hts_pos_t end, std::string_view key,
                 std::string_view reason, std::string_view detail = {}) {
    FILE* fh = sink();
    if (fh == nullptr) return;
    static std::mutex lock;
    const Context& where = context();
    const std::lock_guard<std::mutex> guard(lock);
    std::fprintf(fh, "%s\t%d\t%d\t", where.contig.c_str(), where.chunk, where.pass);
    std::fprintf(fh, "%.*s\t%.*s\t%lld\t%lld\t%.*s\t%.*s\t%.*s\n", static_cast<int>(stage.size()), stage.data(),
                 static_cast<int>(event.size()), event.data(), static_cast<long long>(beg),
                 static_cast<long long>(end), static_cast<int>(key.size()), key.data(),
                 static_cast<int>(reason.size()), reason.data(), static_cast<int>(detail.size()), detail.data());
}

// TYPE:REF_LEN:ALT of a key (ALT cut at 60 letters; a catalog row's ALT is its site id).
inline std::string key_text(const pgphase_collect::VariantKey& key) {
    const char type = key.type == pgphase_collect::VariantType::Snp ? 'S'
                      : key.type == pgphase_collect::VariantType::Insertion ? 'I' : 'D';
    std::string alt = key.alt.size() > 60 ? key.alt.substr(0, 60) + "..." : key.alt;
    return std::string(1, type) + ":" + std::to_string(key.ref_len) + ":" + (alt.empty() ? "-" : alt);
}

// The reference span a candidate's record covers: from its VCF POS (vcf_record_pos) to its last REF base.
inline std::pair<hts_pos_t, hts_pos_t> span(const pgphase_collect::CandidateVariant& c) {
    const hts_pos_t beg = pgphase_collect::vcf_record_pos(c);
    const hts_pos_t last = c.key.pos + std::max(1, c.key.ref_len) - 1;
    return {beg, std::max(beg, last)};
}

// One event about a candidate.
inline void candidate(std::string_view stage, std::string_view event, const pgphase_collect::CandidateVariant& c,
                      std::string_view reason, std::string_view detail = {}) {
    if (!on()) return;
    const auto [beg, end] = span(c);
    emit(stage, event, beg, end, key_text(c.key), reason, detail);
}

// "h1/h2 ps=N cat=C" of a candidate, for details.
// gt= (the alleles haplotypes 1 and 2 carry, -1 none), ps= and cat= of a candidate; prefix names another one's
// (to_, kept_).
inline std::string genotype_text(const pgphase_collect::CandidateVariant& c, const std::string& prefix = "") {
    return prefix + "gt=" + std::to_string(c.hap_to_cons_alle[1]) + "/" + std::to_string(c.hap_to_cons_alle[2]) +
           " " + prefix + "ps=" + std::to_string(c.phase_set) + " " + prefix + "cat=" +
           std::to_string(static_cast<int>(c.counts.category));
}

}  // namespace pgphase_trace
