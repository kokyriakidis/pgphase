// Union gap phasing, stage 1: the alignment solve's heterozygous sites enter
// the graph chunk (injection, hidden complement sites).

#include "union_internal.hpp"

#include "noise_filter.hpp"
#include "trace.hpp"

#include <algorithm>
#include <set>

namespace pgphase_collect {

namespace {

using CandKey = CandidateIdentityKey;

CandKey cand_key_of(const CandidateVariant& cand) {
    return CandKey{cand.key.sort_pos(), static_cast<int>(cand.key.type), cand.key.ref_len,
                   cand.key.alt};
}

}  // namespace

// An MSA-verified site is admitted only if each of its two alleles is seen by
// at least this many reads mapped at kDefaultMinMapq or better (a clean het is
// exempt; see the selection loop).
constexpr int kUnionMinConfidentAlleleReads = 2;
// A site read by more than this multiple of its chunk's median depth is a
// collapsed duplication: its two "alleles" separate paralogous copies, not
// haplotypes. Measured on chr12 at 10.31 Mb: ~150 such sites at twice the local
// depth agreed with each other perfectly, separated the parents at 0.52-0.66,
// and joined the flanking blocks in the wrong orientation.
// 1.5 rather than 1.6: a paralog run at 1.45-1.52 times the depth (chr18
// 73.19 Mb, 70/35 alleles, reads split by parent at 0.53) joined two blocks
// the wrong way (546 reads) once the graph's catalog calls were checked
// against the alignment's (confirm_catalog_snp_calls). 1.4 cost real sites
// (chr12: 125 correct reads to unphased).
constexpr double kUnionMaxDepthRatio = 1.5;
// The longest voting row an alignment call can be another description of:
// the alignment's own noisy-region MSA covers up to 2 kb (max_noisy_reg_len).
// A call inside a longer row -- a snarl spanning a large SV, up to 1 Mb on
// chr20 -- is a variant the graph does not describe at that scale, i.e.
// sample-specific signal, and is injected as a site of its own.
constexpr hts_pos_t kUnionMaxComplementRowSpan = 2000;

/// VCF-form metadata for an alignment candidate the catalog never held, in the
/// convention the graph writer emits and vcf_to_variant_key reads back. An
/// empty REF makes the writer skip the record while index-parallel arrays stay
/// aligned.
static std::pair<GraphSiteMeta, std::vector<int>> synthesize_alignment_site_meta(
        const CandidateVariant& cand, const PhasingChunk& source, const char* contig_name) {
    const auto ref_at = [&source](hts_pos_t pos, int len) -> std::string {
        if (len <= 0 || source.ref_seq.empty()) return std::string();
        const hts_pos_t off = pos - source.ref_beg;
        if (off < 0 || static_cast<size_t>(off) + static_cast<size_t>(len) > source.ref_seq.size())
            return std::string();
        return source.ref_seq.substr(static_cast<size_t>(off), static_cast<size_t>(len));
    };
    const VariantKey& key = cand.key;
    const auto anchored = [&](const std::string& allele) -> std::string {
        if (key.type == VariantType::Snp) return allele;
        if (key.type == VariantType::Insertion && key.ref_len != 0) return allele;
        return ref_at(key.pos - 1, 1) + allele;
    };
    hts_pos_t vcf_pos = key.pos;
    std::string vcf_ref;
    if (key.type == VariantType::Snp) {
        vcf_ref = ref_at(key.pos, std::max(1, key.ref_len));
    } else if (key.type == VariantType::Insertion && key.ref_len == 0) {
        vcf_pos = key.pos - 1;
        vcf_ref = ref_at(vcf_pos, 1);
    } else if (key.type == VariantType::Insertion) {
        vcf_ref = ref_at(key.pos, key.ref_len);
    } else {
        vcf_pos = key.pos - 1;
        vcf_ref = ref_at(vcf_pos, 1) + ref_at(key.pos, key.ref_len);
    }
    const std::string vcf_alt = anchored(key.alt);
    bool usable = !vcf_ref.empty() && !vcf_alt.empty();
    if (usable) {
        const VariantKey round_trip = vcf_to_variant_key(key.tid, vcf_pos, vcf_ref, vcf_alt);
        usable = round_trip.type == key.type && round_trip.pos == key.pos &&
                 round_trip.ref_len == key.ref_len && round_trip.alt == key.alt;
    }
    // An MSA indel whose anchor letter is not the reference's on its haplotype (a SNP right before an insertion or a
    // deletion, aligned as the mismatch then the gap): the consensus's letter there leads the ALT, one complex record,
    // as collect-bam-variation writes it (alt_ref_base). With the reference's letter the SNP was lost: chr12:6,396,535
    // C>CG beside A>G written as A>AG, another haplotype. HG002 HiFi chr12/18/20: true het regions 29,979 -> 30,200,
    // wrong 681 -> 487; reads unchanged.
    std::string written_alt = vcf_alt;
    if (usable && key.type != VariantType::Snp && cand.msa_insertion_alts.empty() &&
        cand.alt_ref_base < 4 && !written_alt.empty()) {
        const char consensus = "ACGT"[cand.alt_ref_base];
        if (std::toupper(static_cast<unsigned char>(written_alt[0])) != consensus) written_alt[0] = consensus;
    }
    std::vector<std::string> alts;
    if (cand.msa_insertion_alts.size() >= 2) {
        for (size_t ai = 0; ai < cand.msa_insertion_alts.size(); ++ai) {
            std::string a = anchored(cand.msa_insertion_alts[ai]);
            if (a.empty()) { alts.clear(); break; }
            // Each allele's own consensus letter leads it, as for a single ALT above.
            const uint8_t anchor = msa_allele_anchor(cand, ai);
            if (key.type != VariantType::Snp && !(key.type == VariantType::Insertion && key.ref_len != 0) && anchor < 4)
                a[0] = "ACGT"[anchor];
            alts.push_back(a);
        }
    }
    if (alts.empty() && usable) alts.push_back(written_alt);
    GraphSiteMeta meta;
    meta.chrom = contig_name != nullptr ? contig_name : std::string();
    if (!usable || alts.empty()) pgphase_trace::candidate("inject", "drop", cand, "unusable_vcf_meta");
    if (usable && !alts.empty()) {
        meta.pos = vcf_pos;
        meta.ref = vcf_ref;
        meta.alts = alts;
    }
    std::vector<int> orig;
    for (int i = 0; i <= static_cast<int>(alts.size()); ++i) orig.push_back(i);
    if (orig.size() < 2) orig = std::vector<int>{0, 1};
    return {std::move(meta), std::move(orig)};
}

// Left-aligned identity of a simple insertion or deletion, so placements of
// one allele at different repeat copies compare equal. Other rows keep their key.
static CandKey left_aligned_indel_key(const CandidateVariant& c, const PhasingChunk& source) {
    VariantKey k = c.key;
    const auto base = [&source](hts_pos_t pos) -> char {
        const hts_pos_t off = pos - source.ref_beg;
        if (off < 0 || static_cast<size_t>(off) >= source.ref_seq.size()) return 'N';
        return static_cast<char>(std::toupper(static_cast<unsigned char>(source.ref_seq[static_cast<size_t>(off)])));
    };
    if (!c.msa_insertion_alts.empty()) return cand_key_of(c);
    if (k.type == VariantType::Insertion && k.ref_len == 0 && !k.alt.empty()) {
        // The inserted bases sit after pos - 1; shift while the base before
        // the insertion point equals the last inserted base.
        while (true) {
            const char prev = base(k.pos - 1);
            if (prev == 'N' || std::toupper(static_cast<unsigned char>(k.alt.back())) != prev) break;
            k.alt = std::string(1, prev) + k.alt.substr(0, k.alt.size() - 1);
            --k.pos;
        }
    } else if (k.type == VariantType::Deletion && k.ref_len > 0 && k.alt.empty()) {
        while (true) {
            const char prev = base(k.pos - 1);
            if (prev == 'N' || prev != base(k.pos + k.ref_len - 1)) break;
            --k.pos;
        }
    }
    return CandKey{k.sort_pos(), static_cast<int>(k.type), k.ref_len, k.alt};
}

namespace {

// Rows that already vote describe this locus. An alignment call overlapping
// one is another description of it: it joins as a hidden complement site
// (no VCF record) carrying the alignment's calls, because the graph often
// leaves reads at such a locus without any allele. Only rows of local span
// (kUnionMaxComplementRowSpan) qualify: a 373 kb snarl at chr18:20,440,896
// swallowed every alignment het for 35 kb (dozens of truth hets LongPhase
// phases), and a paralog cluster then joined across the gap the wrong way.
class VotingRows {
public:
    explicit VotingRows(const PhasingChunk& graph) {
        for (size_t ci = 0; ci < graph.candidates.size(); ++ci) {
            const CandidateVariant& row = graph.candidates[ci];
            if ((row.lcd_var_i_to_cate & kCandGermlineVarCate) == 0) continue;
            if (std::max(1, row.key.ref_len) > kUnionMaxComplementRowSpan) continue;
            spans_.emplace_back(row.key.pos - 1, row.key.pos + std::max(1, row.key.ref_len), ci);
        }
        std::sort(spans_.begin(), spans_.end());
        for (const auto& span : spans_) widest_ = std::max(widest_, std::get<1>(span) - std::get<0>(span));
    }

    // Voting rows overlapping [beg, end).
    std::vector<size_t> over(hts_pos_t beg, hts_pos_t end) const {
        std::vector<size_t> rows;
        auto it = std::lower_bound(spans_.begin(), spans_.end(),
                                   std::make_tuple(beg - widest_, hts_pos_t{0}, size_t{0}));
        for (; it != spans_.end() && std::get<0>(*it) < end; ++it)
            if (std::get<1>(*it) > beg) rows.push_back(std::get<2>(*it));
        return rows;
    }

private:
    std::vector<std::tuple<hts_pos_t, hts_pos_t, size_t>> spans_;
    hts_pos_t widest_ = 0;
};

// Each alignment candidate's read support: reads per allele, of them reads
// mapped at kDefaultMinMapq or better, depth and the chunk's median depth.
struct AlignmentSupport {
    std::vector<std::map<int, int>> observed;
    std::vector<std::map<int, int>> confident;
    // Diagnostic counts: reads the alignment solve labelled, and of those, the
    // ones carrying the allele the solve assigned to their haplotype.
    std::vector<std::array<int, 2>> labelled;
    std::vector<int> depth;
    double median_depth = 0.0;

    explicit AlignmentSupport(const PhasingChunk& bam)
        : observed(bam.candidates.size()), confident(bam.candidates.size()),
          labelled(bam.candidates.size(), {0, 0}), depth(bam.candidates.size(), 0) {
        for (const ReadVariantProfile& prof : bam.read_var_profile) {
            if (prof.read_id < 0 || prof.start_var_idx < 0) continue;
            const ReadRecord& read = bam.reads[static_cast<size_t>(prof.read_id)];
            if (read.is_skipped) continue;
            const int hap = static_cast<size_t>(prof.read_id) < bam.haps.size()
                                ? bam.haps[static_cast<size_t>(prof.read_id)] : 0;
            for (size_t k = 0; k < prof.alleles.size(); ++k) {
                const int a = prof.alleles[k];
                if (a < 0) continue;
                const size_t ci = static_cast<size_t>(prof.start_var_idx) + k;
                ++observed[ci][a];
                if (read.mapq >= kDefaultMinMapq) ++confident[ci][a];
                if (hap == 1 || hap == 2) {
                    ++labelled[ci][0];
                    const auto& cons = bam.candidates[ci].hap_to_cons_alle;
                    if (cons.size() > 2 && cons[static_cast<size_t>(hap)] == a) ++labelled[ci][1];
                }
            }
        }
        std::vector<int> depths;
        for (size_t bi = 0; bi < bam.candidates.size(); ++bi) {
            for (const auto& [a, n] : observed[bi]) depth[bi] += n;
            // Calls kept from noisy regions too long for MSA sit largely in
            // repeats of unusual depth (a VNTR at depth 13 on chr18); counted
            // here they pulled the median down and gated out ordinary sites
            // (chr12:133,056,148 at depth 81).
            if (depth[bi] > 0 && !bam.candidates[bi].unresolved_noisy) depths.push_back(depth[bi]);
        }
        if (!depths.empty()) {
            std::nth_element(depths.begin(), depths.begin() + depths.size() / 2, depths.end());
            median_depth = depths[depths.size() / 2];
        }
    }

    // Confident reads on the two most observed alleles, in that order.
    std::array<int, 2> confident_on_leading(size_t bi) const {
        std::vector<std::pair<int, int>> top;
        for (const auto& [a, n] : observed[bi]) top.emplace_back(n, a);
        std::sort(top.rbegin(), top.rend());
        std::array<int, 2> conf{0, 0};
        for (size_t i = 0; i < 2 && i < top.size(); ++i) {
            const auto it = confident[bi].find(top[i].second);
            if (it != confident[bi].end()) conf[i] = it->second;
        }
        return conf;
    }

    // Reads admitted at low MAPQ can represent a divergent haplotype but
    // cannot create a site on their own: each of the site's two alleles (a 1/2
    // site carries no REF reads) needs confidently mapped support.
    bool confidently_supported(size_t bi) const {
        if (observed[bi].size() < 2) return false;
        const std::array<int, 2> conf = confident_on_leading(bi);
        return conf[0] >= kUnionMinConfidentAlleleReads && conf[1] >= kUnionMinConfidentAlleleReads;
    }
};

// Which alignment candidates enter the graph chunk, and as what.
struct Selection {
    std::vector<size_t> selected;
    std::set<size_t> island_only;    // selected only as island sites
    std::set<size_t> low_mapq_only;  // clean hets admitted by the exemption
    std::unordered_map<size_t, std::vector<size_t>> complement_rows;  // bam candidate -> graph rows it overlaps
    std::set<size_t> replaced_rows;  // graph rows an MSA call of an allele they do not spell replaces in the VCF
    std::set<size_t> unwritten_rows;  // graph rows whose allele the consensus sequences of their MSA region lack
    std::set<size_t> hom_rows;        // graph rows both consensus sequences of their MSA region carry: written 1/1
    std::set<size_t> hom_records;     // MSA homozygous calls in a het MSA region: written 1/1, no vote
    std::unordered_map<size_t, const char*> row_reason;  // why a graph row is unwritten or 1/1 (stage trace)
    std::vector<const char*> decision;  // per candidate (diagnostic only)
};

// Whether an alignment candidate lies strictly inside a region the MSA resolved into two different consensus
// sequences (MsaRegion::het).
static bool in_het_msa_region(const PhasingChunk& bam, const CandidateVariant& c) {
    const hts_pos_t lo = c.key.pos - 1;
    const hts_pos_t hi = c.key.pos + std::max(1, c.key.ref_len);
    for (const PhasingChunk::MsaRegion& region : bam.msa_regions)
        if (region.het() && region.beg < lo && hi < region.end) return true;
    return false;
}

// Clean heterozygotes and MSA-verified noisy heterozygotes, keeping the
// alignment's own category bits. And the MSA's homozygous calls inside a region whose two consensus sequences differ:
// the haplotypes there are spelled only with them. Union phasing otherwise writes no injected homozygous call (an
// alignment-discovered subset of the genome's), and at chr12:17,633,561 the region's consensus sequences (+ATCT, and
// a T run one letter apart) came out as the run lengths alone, both haplotypes missing the insertion.
// The rows the MSA regions already genotyped (sel.unwritten_rows, sel.hom_rows) come in set.
void select_sites(const PhasingChunk& graph, const PhasingChunk& bam, const AlignmentSupport& support,
                  const VotingRows& voting, const std::vector<GraphSiteMeta>& graph_meta, bool record_decisions,
                  Selection& sel) {
    if (record_decisions) sel.decision.assign(bam.candidates.size(), "");
    const auto decide = [&](size_t bi, const char* why) {
        if (record_decisions) sel.decision[bi] = why;
        if (!pgphase_trace::on()) return;
        std::string observed;
        for (const auto& [a, n] : support.observed[bi]) {
            const auto it = support.confident[bi].find(a);
            observed += (observed.empty() ? " reads=" : ",") + std::to_string(a) + ":" + std::to_string(n) + "/" +
                        std::to_string(it == support.confident[bi].end() ? 0 : it->second);
        }
        pgphase_trace::candidate("inject", "decision", bam.candidates[bi], why,
                                 "cat=" + std::to_string(static_cast<int>(bam.candidates[bi].counts.category)) +
                                 (bam.candidates[bi].msa_verified ? " msa" : "") + observed);
    };
    for (size_t bi = 0; bi < bam.candidates.size(); ++bi) {
        const CandidateVariant& c = bam.candidates[bi];
        const VariantCategory cat = c.counts.category;
        const bool clean = cat == VariantCategory::CleanHetSnp || cat == VariantCategory::CleanHetIndel;
        const bool msa = cat == VariantCategory::NoisyCandHet && c.msa_verified;
        if (cat == VariantCategory::NoisyCandHom && c.msa_verified && in_het_msa_region(bam, c) &&
            c.key.pos >= graph.ref_beg && c.key.pos <= graph.ref_end) {
            decide(bi, "hom");
            sel.hom_records.insert(bi);
            sel.selected.push_back(bi);
            continue;
        }
        if (!clean && !msa) { decide(bi, "category"); continue; }
        if ((c.lcd_var_i_to_cate & kCandGermlineVarCate) == 0) { decide(bi, "not_germline"); continue; }
        if (c.key.pos < graph.ref_beg || c.key.pos > graph.ref_end) { decide(bi, "outside_chunk"); continue; }
        // An MSA-verified het these two rules drop enters as an island-only
        // site (island_only): its reads may phase locally, but it never votes
        // or joins. Admitted as an ordinary site instead, the same calls cost
        // -1,258 correct / +1,269 discordant reads (chr12, depth gate: a
        // paralog split joined a block the wrong way) or doubled chr20's arm
        // switches (confident-read rule).
        // A clean het is exempt from the confident-read rule: the depth gate
        // alone screens it. Where a region has a second copy elsewhere in the
        // reference every read loses its MAPQ, yet the reads are this copy's
        // (depth at or below the median, not a collapse): at chr18:21.09-21.10
        // Mb six truth hets read only at MAPQ < 30 were dropped and 25 reads
        // LongPhase phases stayed unphased. Such a site enters flagged
        // low_mapq_only: phased, but it joins no blocks and labels reads only
        // as an unverified call does (union EM). Measured on chr12/18/20 with
        // the echo sites, unverified fallback and per-pair stitching: +155
        // correct / -88 discordant / -67 unphased, +1,569 truth hets phased;
        // arm switches 32 -> 37 (sparse pericentromeric, satellite and VNTR
        // blocks: chr12 34.65-37.19 Mb, chr20 32.76 and 65.51 Mb).
        const bool confident = support.confidently_supported(bi);
        const char* dropped = !clean && !confident ? "few_confident_reads"
                            : (support.median_depth > 0.0 &&
                               support.depth[bi] > kUnionMaxDepthRatio * support.median_depth) ? "depth_gate"
                            : nullptr;
        std::vector<size_t> rows = voting.over(c.key.pos - 1, c.key.pos + std::max(1, c.key.ref_len));
        if (dropped != nullptr) {
            if (msa && rows.empty()) {
                decide(bi, dropped[0] == 'f' ? "island:few_confident_reads" : "island:depth_gate");
                sel.island_only.insert(bi);
                sel.selected.push_back(bi);
            } else {
                decide(bi, dropped);
            }
            continue;
        }
        if (msa && !rows.empty() && graph_meta.size() == graph.candidates.size()) {
            // An MSA het over catalog rows is written itself, the rows (still voting) not: the MSA's haplotypes are
            // the region's reads' own, the rows' alleles what the graph's reads were forced onto (chr20:509,755, a
            // 4 bp deletion the catalog describes only as a SNP beside it) or genotyped from them (chr12:6,797,743, the MSA's het T>A, 33 against 31 reads, under a catalog row
            // the graph's reads called 14 against 50 and wrote 1/1). Hiding a call whose allele the rows spelled, as a
            // complement site beside them, left these written wrong: HG002 HiFi chr12/18/20 true het regions
            // 29,870 -> 29,979, wrong 795 -> 681, missed 539 -> 531; arm reads +25 discordant.
            decide(bi, "replaces");
            // Every row: keeping one an MSA call reproduces (by VCF allele) put a co-located copy beside the call and
            // left 162 truth hets on chr12/18/20 unphased.
            sel.replaced_rows.insert(rows.begin(), rows.end());
            if (!confident) sel.low_mapq_only.insert(bi);
            sel.selected.push_back(bi);
            continue;
        }
        decide(bi, rows.empty() ? "injected" : "complement");
        if (!rows.empty()) sel.complement_rows.emplace(bi, std::move(rows));
        if (!confident) sel.low_mapq_only.insert(bi);
        sel.selected.push_back(bi);
    }
    // A complement whose rows another MSA call replaced, or that write no record (unwritten_rows): no row writes it
    // any longer, so it is written itself, its rows replaced as a replacing call's are (an unwritten one keeps its
    // REF, and its reads' calls). Repeated: each promotion can replace rows another complement leans on
    // (chr12:322,418-419, C>G and T>G under one catalog row; the MSA's 322,418 replaced the row and its 322,419,
    // hidden behind that row, went unwritten). HG002 HiFi chr12/18/20: true het regions 29,702 -> 29,870, wrong
    // 1,014 -> 795; reads unchanged. An unwritten row: chr20:5,460,922, the MSA's 1/2 deletion of TT and T under a
    // catalog SNP neither consensus carries.
    for (bool changed = true; changed;) {
        changed = false;
        for (auto it = sel.complement_rows.begin(); it != sel.complement_rows.end();) {
            bool orphaned = false;
            for (const size_t r : it->second)
                orphaned = orphaned || sel.replaced_rows.count(r) != 0 || sel.unwritten_rows.count(r) != 0;
            if (!orphaned) {
                ++it;
                continue;
            }
            decide(it->first, "promoted");
            for (const size_t r : it->second)
                if (sel.unwritten_rows.count(r) == 0 && sel.replaced_rows.insert(r).second)
                    changed = true;
            it = sel.complement_rows.erase(it);
        }
    }
}

// Diagnostic, off by default: PGPHASE_DIAG=DIR writes each alignment
// candidate's injection decision (inj.TID.BEG.tsv), each read's calls in the
// alignment solve (alnreads.TID.BEG.tsv) and the graph chunk's rows before the
// merge (graph.TID.BEG.tsv).
void write_injection_diagnostics(const char* diag_dir, const GraphChunkBuildResult& graph_chunk,
                                 const PhasingChunk& bam, const AlignmentSupport& support, const Selection& sel) {
    const PhasingChunk& graph = graph_chunk.chunk;
    const std::string tag = std::to_string(graph.region.tid) + "." + std::to_string(graph.region.beg) + ".tsv";
    const auto type_of = [](const CandidateVariant& c) {
        return c.key.type == VariantType::Snp ? 'S' : c.key.type == VariantType::Insertion ? 'I' : 'D';
    };
    if (FILE* fh = std::fopen((std::string(diag_dir) + "/inj." + tag).c_str(), "w")) {
        // VCF-convention position, key position, type, REF length, ALT,
        // category code, MSA-verified, unresolved_noisy, depth, chunk median
        // depth, confident reads on the two leading alleles, decision, reads
        // the alignment solve labelled, and of those the ones carrying their
        // haplotype's consensus allele.
        for (size_t bi = 0; bi < bam.candidates.size(); ++bi) {
            const CandidateVariant& c = bam.candidates[bi];
            const std::array<int, 2> conf = support.confident_on_leading(bi);
            std::fprintf(fh, "%lld\t%lld\t%c\t%d\t%s\t%d\t%d\t%d\t%d\t%.0f\t%d,%d\t%s\t%d\t%d\n",
                         static_cast<long long>(c.key.sort_pos()), static_cast<long long>(c.key.pos), type_of(c),
                         c.key.ref_len, c.key.alt.substr(0, 40).c_str(), static_cast<int>(c.counts.category),
                         c.msa_verified ? 1 : 0, c.unresolved_noisy ? 1 : 0, support.depth[bi], support.median_depth,
                         conf[0], conf[1], sel.decision[bi], support.labelled[bi][0], support.labelled[bi][1]);
        }
        std::fclose(fh);
    }
    if (FILE* fh = std::fopen((std::string(diag_dir) + "/alnreads." + tag).c_str(), "w")) {
        // Name, MAPQ, skipped flag, the solve's haplotype and phase set, then
        // VCF-convention position:type:allele for every het, hom or
        // repeat-indel candidate the read has an allele at.
        constexpr uint32_t kListed = kCandCleanHetSnp | kCandCleanHetIndel | kCandCleanHom |
                                     kCandNoisyCandHet | kCandNoisyCandHom | kCandNonAnchorHet;
        for (const ReadVariantProfile& prof : bam.read_var_profile) {
            if (prof.read_id < 0) continue;
            const size_t ri = static_cast<size_t>(prof.read_id);
            const ReadRecord& read = bam.reads[ri];
            std::fprintf(fh, "%s\t%d\t%d\t%d\t%lld", read.qname.c_str(), read.mapq, read.is_skipped ? 1 : 0,
                         ri < bam.haps.size() ? bam.haps[ri] : 0,
                         static_cast<long long>(ri < bam.phase_sets.size() ? bam.phase_sets[ri] : 0));
            for (size_t k = 0; prof.start_var_idx >= 0 && k < prof.alleles.size(); ++k) {
                if (prof.alleles[k] < 0) continue;
                const CandidateVariant& c = bam.candidates[static_cast<size_t>(prof.start_var_idx) + k];
                const bool repeat = c.counts.category == VariantCategory::RepeatHetIndel;
                if ((c.lcd_var_i_to_cate & kListed) == 0 && !repeat) continue;
                std::fprintf(fh, "\t%lld:%c:%d", static_cast<long long>(c.key.sort_pos()), type_of(c), prof.alleles[k]);
            }
            std::fputc('\n', fh);
        }
        std::fclose(fh);
    }
    if (FILE* fh = std::fopen((std::string(diag_dir) + "/graph." + tag).c_str(), "w")) {
        // VCF-convention position, key position, type, REF length, ALT,
        // category code, voting (germline), catalog id, allele coverages.
        for (size_t ci = 0; ci < graph.candidates.size(); ++ci) {
            const CandidateVariant& c = graph.candidates[ci];
            std::string covs;
            for (const int n : c.counts.alle_covs) covs += (covs.empty() ? "" : ",") + std::to_string(n);
            const std::string id = ci < graph_chunk.site_ids.size() ? graph_chunk.site_ids[ci] : std::string();
            std::fprintf(fh, "%lld\t%lld\t%c\t%d\t%s\t%d\t%d\t%s\t%s\n",
                         static_cast<long long>(c.key.sort_pos()), static_cast<long long>(c.key.pos), type_of(c),
                         c.key.ref_len, c.key.alt.substr(0, 40).c_str(), static_cast<int>(c.counts.category),
                         (c.lcd_var_i_to_cate & kCandGermlineVarCate) != 0 ? 1 : 0,
                         id.empty() ? "." : id.c_str(), covs.empty() ? "." : covs.c_str());
        }
        std::fclose(fh);
    }
}

// Catalog rows inside a region the MSA resolved where neither of its consensus sequences differs from the reference:
// the graph's reads put there an allele the region's own reads, aligned as the MSA aligns them, do not carry
// (chr12:5,090,686, a SNP of 37 against 26 graph reads, 7 bp from the MSA's insertion, no truth variant;
// chr12:7,225,347, a 6-letter catalog allele ending a base before the MSA's SNP at 7,225,350). Written no record
// (GraphSiteMeta::no_record), still voting and serving the reads' calls. HG002 HiFi chr12/18/20, plegma: true het
// regions 29,446 -> 29,702, wrong 1,265 -> 1,014, false 205 -> 184, missed 533 -> 539; reads unchanged. Clearing
// the rows' REF instead, as replaced rows are, cost 36 discordant arm reads (their read calls went too).
void suppress_rows_msa_reads_as_reference(const PhasingChunk& graph, const PhasingChunk& bam, Selection& sel) {
    if (bam.msa_regions.empty()) return;
    for (size_t ci = 0; ci < graph.candidates.size(); ++ci) {
        const CandidateVariant& c = graph.candidates[ci];
        const hts_pos_t lo = c.key.pos - 1;
        const hts_pos_t hi = c.key.pos + std::max(1, c.key.ref_len);
        for (const PhasingChunk::MsaRegion& region : bam.msa_regions) {
            if (lo <= region.beg || hi >= region.end) continue;
            bool differs = false;
            for (const Interval& d : region.differs) differs = differs || (d.beg <= hi && lo <= d.end);
            if (!differs && sel.unwritten_rows.insert(ci).second) sel.row_reason[ci] = "msa_reads_as_reference";
            break;
        }
    }
}

// Catalog substitution rows (a SNP, or an MNP the writer splits into its SNPs) against what the MSA's consensus
// sequences carry at each base the row changes strictly inside a resolved region (MsaRegion::spelled). A change no
// consensus carries: the row is no haplotype of the region's reads, written no record and still voting, as
// suppress_rows_msa_reads_as_reference does for rows inside a region read as reference -- here also rows that run
// past the region and rows where a consensus differs otherwise (chr12:8,726,713-753, a 41 bp catalog walk with 4 of
// 16 graph reads written as 27 het SNPs, its first 11 bases the region's, reference in both consensus sequences).
// Every change carried by both consensus sequences of a region whose two differ (MsaRegion::het): written 1/1, still
// voting. The graph's reads, forced onto the catalog's walks beside an indel, call such a site het (chr12:15,541,161
// A>T, 41 against 32 graph reads, 79 against 8 alignments, the truth 1|1; chr18:8,187,001 A>G, 23 against 54 and 35
// against 0). Only in a het region: the MSA's other homozygous calls are written there too (select_sites), and a SNP
// alone can be half of an indel's alignment (chr20:40,719,136 C>T, the truth's deletion of that C, aligned as the
// mismatch and a T deleted further along the run).
void genotype_substitution_rows(const GraphChunkBuildResult& graph_chunk, const PhasingChunk& bam, Selection& sel) {
    const PhasingChunk& graph = graph_chunk.chunk;
    if (bam.msa_regions.empty() || graph_chunk.site_meta.size() != graph.candidates.size()) return;
    const auto region_at = [&bam](hts_pos_t p) -> const PhasingChunk::MsaRegion* {
        for (const PhasingChunk::MsaRegion& region : bam.msa_regions)
            if (region.beg < p && p < region.end && !region.spelled.empty()) return &region;
        return nullptr;
    };
    const auto upper = [](char c) { return static_cast<char>(std::toupper(static_cast<unsigned char>(c))); };
    for (size_t ci = 0; ci < graph.candidates.size(); ++ci) {
        const GraphSiteMeta& meta = graph_chunk.site_meta[ci];
        if (graph.candidates[ci].counts.alle_covs.size() != 2 || meta.ref.empty() || sel.unwritten_rows.count(ci) != 0)
            continue;
        const std::vector<int>* orig =
            ci < graph_chunk.site_allele_orig_idx.size() ? &graph_chunk.site_allele_orig_idx[ci] : nullptr;
        const int alt_i = (orig != nullptr && orig->size() > 1 ? (*orig)[1] : 1) - 1;
        if (alt_i < 0 || alt_i >= static_cast<int>(meta.alts.size())) continue;
        hts_pos_t pos = meta.pos;
        std::string ref = meta.ref, alt = meta.alts[static_cast<size_t>(alt_i)];
        trim_to_minimal_vcf(pos, ref, alt);
        if (ref.size() != alt.size() || ref == alt) continue;
        int changed = 0, by_none = 0, by_both = 0;
        for (size_t i = 0; i < ref.size(); ++i) {
            if (upper(ref[i]) == upper(alt[i])) continue;
            const hts_pos_t p = pos + static_cast<hts_pos_t>(i);
            const PhasingChunk::MsaRegion* region = region_at(p);
            if (region == nullptr) continue;
            ++changed;
            const int n = region->carriers_of_snp(p, upper(alt[i]));
            if (n == 0) ++by_none;
            else if (n == 2 && region->het()) ++by_both;
        }
        if (changed == 0) continue;
        if (by_none > 0) {
            if (sel.unwritten_rows.insert(ci).second) sel.row_reason[ci] = "msa_lacks_change";
        } else if (by_both == changed && sel.hom_rows.insert(ci).second) {
            sel.row_reason[ci] = "msa_homozygous";
        }
    }
}

// Catalog indel rows against the MSA's consensus sequences, as genotype_substitution_rows checks substitutions: a row
// whose haplotype, over every place in the repeat tract its indel can sit, neither consensus spells is no haplotype of
// the region's reads; written no record, still voting. The SNP check alone left such rows written beside the MSA's
// own calls (chr20:14,371,298 GA>G, phased, in a region whose two consensus sequences carry no deletion there).
void unwrite_indel_rows_msa_lacks(const GraphChunkBuildResult& graph_chunk, const PhasingChunk& bam, Selection& sel) {
    const PhasingChunk& graph = graph_chunk.chunk;
    if (bam.msa_regions.empty() || graph_chunk.site_meta.size() != graph.candidates.size()) return;
    const auto upper = [](char c) { return static_cast<char>(std::toupper(static_cast<unsigned char>(c))); };
    const auto base = [&](hts_pos_t pos) -> char {
        const hts_pos_t off = pos - bam.ref_beg;
        if (off < 0 || static_cast<size_t>(off) >= bam.ref_seq.size()) return 'N';
        return upper(bam.ref_seq[static_cast<size_t>(off)]);
    };
    for (size_t ci = 0; ci < graph.candidates.size(); ++ci) {
        const GraphSiteMeta& meta = graph_chunk.site_meta[ci];
        if (graph.candidates[ci].counts.alle_covs.size() != 2 || meta.ref.empty() || sel.unwritten_rows.count(ci) != 0)
            continue;
        const std::vector<int>* orig =
            ci < graph_chunk.site_allele_orig_idx.size() ? &graph_chunk.site_allele_orig_idx[ci] : nullptr;
        const int alt_i = (orig != nullptr && orig->size() > 1 ? (*orig)[1] : 1) - 1;
        if (alt_i < 0 || alt_i >= static_cast<int>(meta.alts.size())) continue;
        hts_pos_t pos = meta.pos;
        std::string ref = meta.ref, alt = meta.alts[static_cast<size_t>(alt_i)];
        for (char& ch : ref) ch = upper(ch);
        for (char& ch : alt) ch = upper(ch);
        trim_to_minimal_vcf(pos, ref, alt);
        // A simple indel after its padding base.
        if (ref.size() == alt.size() || std::min(ref.size(), alt.size()) != 1 || ref[0] != alt[0]) continue;
        const bool deletion = ref.size() > alt.size();
        std::string unit = deletion ? ref.substr(1) : alt.substr(1);
        const hts_pos_t len = deletion ? static_cast<hts_pos_t>(unit.size()) : 0;
        // Its leftmost and rightmost placements (q: the first deleted base, or the base the insertion precedes).
        hts_pos_t lq = pos + 1, rq = pos + 1;
        std::string left_unit = unit;
        while (base(lq - 1) != 'N' && base(lq - 1) == left_unit.back()) {
            left_unit = std::string(1, left_unit.back()) + left_unit.substr(0, left_unit.size() - 1);
            --lq;
        }
        for (std::string u = unit; base(rq + len) != 'N' && base(rq + (deletion ? len : 0)) == u.front();) {
            u = u.substr(1) + u.front();
            ++rq;
        }
        // The window: the padding base of the leftmost placement to the last base the rightmost one touches.
        const hts_pos_t w0 = lq - 1, w1 = deletion ? rq + len - 1 : rq - 1;
        const PhasingChunk::MsaRegion* region = nullptr;
        for (const PhasingChunk::MsaRegion& r : bam.msa_regions)
            if (r.beg < w0 && w1 < r.end && !r.spelled.empty()) { region = &r; break; }
        if (region == nullptr) continue;
        std::string row(1, base(w0));
        if (!deletion) row += left_unit;
        for (hts_pos_t x = deletion ? lq + len : lq; x <= w1; ++x) row.push_back(base(x));
        bool carried = false;
        for (const std::vector<std::string>& cons : region->spelled) {
            std::string seq;
            for (hts_pos_t x = w0; x <= w1; ++x) seq += cons[static_cast<size_t>(x - region->beg)];
            for (char& ch : seq) ch = upper(ch);
            carried = carried || seq == row;
        }
        if (!carried && sel.unwritten_rows.insert(ci).second) sel.row_reason[ci] = "msa_lacks_indel";
    }
}

// One site per left-aligned allele: other placements become aliases whose
// read calls fold into the canonical row. Regular sites come first: an
// island-only call never absorbs one, and one at an allele a regular site
// already carries is dropped. Returns alias -> canonical candidate.
std::unordered_map<size_t, size_t> fold_aliases(Selection& sel, const PhasingChunk& bam) {
    std::unordered_map<size_t, size_t> alias_of;
    std::map<CandKey, size_t> canonical;
    std::vector<size_t> kept;
    std::stable_partition(sel.selected.begin(), sel.selected.end(),
                          [&](size_t bi) { return sel.island_only.count(bi) == 0; });
    for (const size_t bi : sel.selected) {
        const auto [it, inserted] = canonical.emplace(left_aligned_indel_key(bam.candidates[bi], bam), bi);
        if (inserted) {
            kept.push_back(bi);
        } else if (sel.island_only.count(bi) == 0) {
            alias_of.emplace(bi, it->second);
            pgphase_trace::candidate("inject", "merge", bam.candidates[bi], "alias_folded",
                                     "into " + pgphase_trace::key_text(bam.candidates[it->second].key));
        } else {
            pgphase_trace::candidate("inject", "drop", bam.candidates[bi], "island_duplicate_dropped",
                                     "of " + pgphase_trace::key_text(bam.candidates[it->second].key));
        }
    }
    sel.selected = std::move(kept);
    return alias_of;
}

}  // namespace

size_t inject_alignment_sites(GraphChunkBuildResult& graph_chunk, const PhasingChunk& bam,
                              const char* contig_name) {
    PhasingChunk& graph = graph_chunk.chunk;
    const VotingRows voting(graph);
    const AlignmentSupport support(bam);
    static const char* const diag_dir = std::getenv("PGPHASE_DIAG");
    Selection sel;
    suppress_rows_msa_reads_as_reference(graph, bam, sel);
    genotype_substitution_rows(graph_chunk, bam, sel);
    unwrite_indel_rows_msa_lacks(graph_chunk, bam, sel);
    select_sites(graph, bam, support, voting, graph_chunk.site_meta, diag_dir != nullptr, sel);
    for (const size_t ci : sel.unwritten_rows)
        if (ci < graph_chunk.site_meta.size()) graph_chunk.site_meta[ci].no_record = true;
    for (const size_t ci : sel.hom_rows)
        if (ci < graph_chunk.site_meta.size()) graph_chunk.site_meta[ci].msa_hom = true;
    if (pgphase_trace::on() && graph_chunk.site_meta.size() == graph.candidates.size()) {
        const auto row = [&](size_t ci, std::string_view event, std::string_view reason) {
            const GraphSiteMeta& m = graph_chunk.site_meta[ci];
            std::string key = m.ref + ">";
            for (size_t i = 0; i < m.alts.size(); ++i) key += (i ? "," : "") + m.alts[i].substr(0, 60);
            pgphase_trace::emit("inject", event, m.pos, m.pos + std::max<hts_pos_t>(1, m.ref.size()) - 1, key, reason,
                                ci < graph_chunk.site_ids.size() ? graph_chunk.site_ids[ci] : std::string());
        };
        for (const auto& [ci, reason] : sel.row_reason) row(ci, sel.hom_rows.count(ci) ? "regenotype" : "drop", reason);
        for (const size_t ci : sel.replaced_rows) row(ci, "drop", "row_replaced_by_msa_call");
    }
    if (diag_dir != nullptr) write_injection_diagnostics(diag_dir, graph_chunk, bam, support, sel);
    if (sel.selected.empty()) return 0;
    const std::unordered_map<size_t, size_t> alias_of = fold_aliases(sel, bam);

    // Merge in key order; every index-parallel array is rebuilt with the table.
    struct Slot { CandKey key; long graph_index; long bam_index; };
    std::vector<Slot> slots;
    slots.reserve(graph.candidates.size() + sel.selected.size());
    for (size_t ci = 0; ci < graph.candidates.size(); ++ci)
        slots.push_back(Slot{cand_key_of(graph.candidates[ci]), static_cast<long>(ci), -1});
    for (const size_t bi : sel.selected) {
        CandidateVariant relabeled = bam.candidates[bi];
        relabeled.key.tid = graph.region.tid;
        slots.push_back(Slot{cand_key_of(relabeled), -1, static_cast<long>(bi)});
    }
    std::stable_sort(slots.begin(), slots.end(),
                     [](const Slot& a, const Slot& b) { return a.key < b.key; });
    const bool had_meta = graph_chunk.site_meta.size() == graph.candidates.size();
    CandidateTable merged;
    std::vector<std::string> merged_ids;
    std::vector<GraphSiteMeta> merged_meta;
    std::vector<std::vector<int>> merged_orig;
    std::vector<long> graph_to_new(graph.candidates.size(), -1);
    std::unordered_map<size_t, size_t> bam_to_new;
    for (const Slot& slot : slots) {
        const size_t ni = merged.size();
        if (slot.graph_index >= 0) {
            const size_t ci = static_cast<size_t>(slot.graph_index);
            graph_to_new[ci] = static_cast<long>(ni);
            merged.push_back(graph.candidates[ci]);
            merged_ids.push_back(ci < graph_chunk.site_ids.size() ? graph_chunk.site_ids[ci] : std::string());
            merged_meta.push_back(had_meta ? graph_chunk.site_meta[ci] : GraphSiteMeta{});
            if (sel.replaced_rows.count(ci) != 0) {
                merged_meta.back().ref.clear();  // empty REF: the writer emits no record
                merged_meta.back().alts.clear();
            }
            merged_orig.push_back(ci < graph_chunk.site_allele_orig_idx.size()
                                      ? graph_chunk.site_allele_orig_idx[ci] : std::vector<int>{});
            continue;
        }
        const size_t bi = static_cast<size_t>(slot.bam_index);
        bam_to_new.emplace(bi, ni);
        CandidateVariant cand = bam.candidates[bi];
        cand.key.tid = graph.region.tid;
        cand.bam_injected = true;
        cand.alignment_verified = true;
        // Phase comes from the union's own solve, not the alignment's gauge.
        cand.phase_set = kUnsetCandidatePhaseSet;
        cand.hap_to_cons_alle[1] = -1;
        cand.hap_to_cons_alle[2] = -1;
        cand.hap_alt = 0;
        cand.hap_ref = 0;
        if (cand.counts.alle_covs.size() < 2)
            cand.counts.alle_covs = {cand.counts.ref_cov, cand.counts.alt_cov};
        auto [meta, orig] = synthesize_alignment_site_meta(bam.candidates[bi], bam, contig_name);
        cand.low_mapq_only = sel.low_mapq_only.count(bi) != 0;
        if (sel.island_only.count(bi) != 0) {
            cand.island_only = true;
            cand.lcd_var_i_to_cate &= ~kCandGermlineVarCate;
        }
        if (sel.hom_records.count(bi) != 0) {
            cand.lcd_var_i_to_cate &= ~kCandGermlineVarCate;
            cand.hap_to_cons_alle[1] = 1;
            cand.hap_to_cons_alle[2] = 1;
            meta.msa_hom = true;
        }
        if (sel.complement_rows.count(bi) != 0 || sel.island_only.count(bi) != 0) {
            meta = GraphSiteMeta{};  // empty REF: the writer emits no record
            meta.chrom = contig_name != nullptr ? contig_name : std::string();
        }
        merged.push_back(std::move(cand));
        merged_ids.emplace_back();
        merged_meta.push_back(std::move(meta));
        merged_orig.push_back(std::move(orig));
    }

    // Alignment observations at the injected sites, keyed by read name. Reads
    // the GAF never had are added: they are what carries a catalog gap.
    std::unordered_map<std::string_view, size_t> graph_read_by_qname;
    graph_read_by_qname.reserve(graph.reads.size());
    for (size_t ri = 0; ri < graph.reads.size(); ++ri)
        graph_read_by_qname.try_emplace(graph.reads[ri].qname, ri);
    std::vector<std::vector<std::tuple<size_t, int, int>>> added_obs(graph.reads.size());
    const size_t old_profiles = graph.reads.size();
    std::unordered_map<size_t, const std::vector<size_t>*> complement_of_new;  // merged index -> graph rows
    for (const auto& [bi, rows] : sel.complement_rows) {
        const auto it = bam_to_new.find(bi);
        if (it != bam_to_new.end()) complement_of_new.emplace(it->second, &rows);
    }
    for (const ReadVariantProfile& profile : bam.read_var_profile) {
        if (profile.read_id < 0 || profile.start_var_idx < 0) continue;
        const ReadRecord& source_read = bam.reads[static_cast<size_t>(profile.read_id)];
        if (source_read.is_skipped) continue;
        std::vector<std::tuple<size_t, int, int>> calls;
        for (size_t k = 0; k < profile.alleles.size(); ++k) {
            if (profile.alleles[k] < 0) continue;
            size_t source_i = static_cast<size_t>(profile.start_var_idx) + k;
            const auto alias = alias_of.find(source_i);
            if (alias != alias_of.end()) source_i = alias->second;
            const auto found = bam_to_new.find(source_i);
            if (found == bam_to_new.end()) continue;
            if (sel.hom_records.count(source_i) != 0) continue;  // written 1/1, no read calls
            const int qi = k < profile.alt_qi.size() ? profile.alt_qi[k] : 0;
            // An alias and its canonical row describe one allele: ALT wins.
            const auto prior = std::find_if(calls.begin(), calls.end(),
                [&](const auto& call) { return std::get<0>(call) == found->second; });
            if (prior == calls.end())
                calls.emplace_back(found->second, profile.alleles[k], qi);
            else if (profile.alleles[k] > std::get<1>(*prior))
                *prior = std::make_tuple(found->second, profile.alleles[k], qi);
        }
        if (calls.empty()) continue;
        auto read_i = graph_read_by_qname.find(source_read.qname);
        size_t ri;
        if (read_i != graph_read_by_qname.end()) {
            ri = read_i->second;
        } else {
            ReadRecord read;
            read.tid = graph.region.tid;
            read.qname = source_read.qname;
            read.mapq = source_read.mapq;
            read.beg = source_read.beg;
            read.end = source_read.end;
            read.reverse = source_read.reverse;
            ri = graph.reads.size();
            graph.reads.push_back(std::move(read));
            graph.read_var_profile.push_back(ReadVariantProfile{});
            added_obs.emplace_back();
            graph_read_by_qname.try_emplace(graph.reads.back().qname, ri);
        }
        for (auto& call : calls) {
            const auto comp = complement_of_new.find(std::get<0>(call));
            // A hidden complement site takes calls only from reads the graph
            // left silent: the read speaks through the graph row already.
            if (comp != complement_of_new.end() && ri < old_profiles) {
                const ReadVariantProfile& old = graph.read_var_profile[ri];
                bool graph_call = false;
                for (const size_t row : *comp->second) {
                    if (old.start_var_idx < 0 || static_cast<int>(row) < old.start_var_idx ||
                        static_cast<int>(row) > old.end_var_idx) continue;
                    if (old.alleles[row - static_cast<size_t>(old.start_var_idx)] >= 0) graph_call = true;
                }
                if (graph_call) continue;
            }
            added_obs[ri].push_back(call);
        }
    }

    // Rebuild each profile over the merged index space, every channel included.
    for (size_t ri = 0; ri < graph.reads.size(); ++ri) {
        const ReadVariantProfile& old = graph.read_var_profile[ri];
        std::map<size_t, std::array<int, 5>> by_index;  // allele, qi, graph, bam, bam_qi
        const auto slot = [&by_index](size_t i) -> std::array<int, 5>& {
            return by_index.try_emplace(i, std::array<int, 5>{-1, 0, -1, -1, 0}).first->second;
        };
        std::map<size_t, uint8_t> qualities;
        if (old.start_var_idx >= 0) {
            for (size_t k = 0; k < old.alleles.size(); ++k) {
                const size_t ci = static_cast<size_t>(old.start_var_idx) + k;
                if (ci >= graph_to_new.size() || graph_to_new[ci] < 0) continue;
                const size_t ni = static_cast<size_t>(graph_to_new[ci]);
                const bool any = old.alleles[k] >= 0 ||
                    (k < old.graph_alleles.size() && old.graph_alleles[k] >= 0) ||
                    (k < old.bam_alleles.size() && old.bam_alleles[k] != -1);
                if (!any) continue;
                auto& s = slot(ni);
                s[0] = old.alleles[k];
                s[1] = k < old.alt_qi.size() ? old.alt_qi[k] : 0;
                s[2] = k < old.graph_alleles.size() ? old.graph_alleles[k] : -1;
                s[3] = k < old.bam_alleles.size() ? old.bam_alleles[k] : -1;
                s[4] = k < old.bam_qi.size() ? old.bam_qi[k] : 0;
                if (k < old.bam_base_qualities.size() && old.bam_base_qualities[k] > 0)
                    qualities.emplace(ni, old.bam_base_qualities[k]);
            }
        }
        for (const auto& [ni, allele, qi] : added_obs[ri]) {
            auto& s = slot(ni);
            s[0] = allele;
            s[1] = qi;
            s[3] = allele;
            s[4] = qi;
        }
        ReadVariantProfile prof;
        prof.read_id = static_cast<int>(ri);
        prof.bam_mapq = old.bam_mapq;
        if (!by_index.empty()) {
            prof.start_var_idx = static_cast<int>(by_index.begin()->first);
            prof.end_var_idx = static_cast<int>(by_index.rbegin()->first);
            const size_t span = static_cast<size_t>(prof.end_var_idx - prof.start_var_idx + 1);
            prof.alleles.assign(span, -1);
            prof.alt_qi.assign(span, 0);
            const bool graph_channel = !old.graph_alleles.empty();
            const bool bam_channel = !old.bam_alleles.empty() || !added_obs[ri].empty();
            if (graph_channel) prof.graph_alleles.assign(span, -1);
            if (bam_channel) { prof.bam_alleles.assign(span, -1); prof.bam_qi.assign(span, 0); }
            if (!qualities.empty()) prof.bam_base_qualities.assign(span, 0);
            for (const auto& [ni, s] : by_index) {
                const size_t off = ni - static_cast<size_t>(prof.start_var_idx);
                prof.alleles[off] = s[0];
                prof.alt_qi[off] = s[1];
                if (graph_channel) prof.graph_alleles[off] = s[2];
                if (bam_channel) { prof.bam_alleles[off] = s[3]; prof.bam_qi[off] = s[4]; }
            }
            for (const auto& [ni, q] : qualities)
                prof.bam_base_qualities[ni - static_cast<size_t>(prof.start_var_idx)] = q;
        }
        graph.read_var_profile[ri] = std::move(prof);
    }
    graph.candidates = std::move(merged);
    graph_chunk.site_ids = std::move(merged_ids);
    graph_chunk.site_meta = std::move(merged_meta);
    graph_chunk.site_allele_orig_idx = std::move(merged_orig);
    graph_chunk.candidate_allele_contrasts.assign(graph.candidates.size(), std::nullopt);
    graph.haps.resize(graph.reads.size(), 0);
    graph.phase_sets.resize(graph.reads.size(), kUnphasedReadPhaseSet);
    rebuild_read_var_cr(graph);
    return sel.selected.size();
}

std::vector<char> collapsed_duplication_sites(const PhasingChunk& bam) {
    const AlignmentSupport support(bam);
    std::vector<char> out(bam.candidates.size(), 0);
    for (size_t bi = 0; bi < bam.candidates.size(); ++bi)
        out[bi] = support.median_depth > 0.0 && support.depth[bi] > kUnionMaxDepthRatio * support.median_depth;
    return out;
}

} // namespace pgphase_collect
