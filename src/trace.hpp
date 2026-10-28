// Stage trace (diagnostic, off by default): PGPHASE_TRACE=DIR writes DIR/trace.tsv, one line per decision that
// drops, merges, re-keys or re-genotypes a candidate, a site or an output record, with a reason code, so a variant can
// be followed through the pipeline (pgphase-dev scripts/dev/harness). Columns: stage, event, beg, end (the reference
// span, 1-based inclusive, as the record would be written), key (TYPE:REF_LEN:ALT, or the VCF REF>ALT where the
// caller has it), reason, detail.

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
        if (fh != nullptr) std::fputs("stage\tevent\tbeg\tend\tkey\treason\tdetail\n", fh);
        return fh;
    }();
    return file;
}

// Whether tracing is on: callers building a detail string test this first.
inline bool on() { return sink() != nullptr; }

inline void emit(std::string_view stage, std::string_view event, hts_pos_t beg, hts_pos_t end, std::string_view key,
                 std::string_view reason, std::string_view detail = {}) {
    FILE* fh = sink();
    if (fh == nullptr) return;
    static std::mutex lock;
    const std::lock_guard<std::mutex> guard(lock);
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
inline std::string genotype_text(const pgphase_collect::CandidateVariant& c) {
    return std::to_string(c.hap_to_cons_alle[1]) + "/" + std::to_string(c.hap_to_cons_alle[2]) +
           " ps=" + std::to_string(c.phase_set) + " cat=" + std::to_string(static_cast<int>(c.counts.category));
}

}  // namespace pgphase_trace
