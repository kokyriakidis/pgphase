// Union gap phasing, stage 1a: the injection plan. Which of the alignment solve's candidates enter the graph chunk,
// as what, and how the catalog's rows are written beside them (plan_injection); union_inject.cpp applies it.

#include "union_select.hpp"

#include "noise_filter.hpp"
#include "trace.hpp"

#include <algorithm>
#include <cstdlib>
#include <optional>
#include <map>
#include <set>

namespace pgphase_collect {

// An MSA-verified site is admitted only if each of its two alleles is seen by
// at least this many reads mapped at kConfidentReadMinMapq or better (a clean het is
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
// A depth-gated MSA het is written (unphased) only below this depth over the median: above it the reads are mostly
// collapsed paralogs. In the arms of chr12/17/18/19/20/21, islands at 2x or more were 309 false het loci against 32
// truth hets (271 of the false in chr12:10.30-10.40 Mb, the NK receptor cluster); between 1.5x and 2x, 103 against 247.
constexpr double kIslandWriteMaxDepthRatio = 2.0;
// (A collapsed copy's island the union EM phases badly is not written either: union_em.cpp.)
//
// Two islands the gates drop are written all the same (unphased) where the call's two alleles lie one on each
// haplotype the alignment solve labelled its reads with (alleles_split_by_haplotype: kIslandSplitMinReads reads each
// side, at most kIslandSplitMaxStray of a haplotype's reads on the other allele) and its ALT carries
// kIslandMinAltShare-kIslandMaxAltShare of its reads:
//   - a clean het the depth gate drops, under kIslandWriteMaxDepthRatio (taking the MSA's island path; 40aa9c5 wrote
//     every one, and the paralog runs of chr16/chr18 came out as false hets);
//   - an MSA het the few-confident-reads gate drops at full depth (kUnionHomFullDepthRatio to kUnionMaxDepthRatio of
//     the median): a region giraffe maps at low MAPQ though the reads are this copy's (chr6:26.5 Mb, chr16, chr17).
// Inside the benchmark of chr12/16/17/18/19/20/21, as an offline ledger counted them: clean depth-gated hets 61 right
// against 2 false loci, the low-MAPQ MSA hets 280 right against 10, no chromosome losing.
constexpr int kIslandSplitMinReads = 5;
constexpr double kIslandSplitMaxStray = 0.05;
constexpr double kIslandMinAltShare = 0.3;
constexpr double kIslandMaxAltShare = 0.7;

// The share of its MSA region's reads an MSA hom's clusters must take for it to be written as a fallback hom.
constexpr double kUnionFallbackHomMinClusteredShare = 0.8;
// Confident reads (MAPQ >= kConfidentReadMinMapq) an alignment hom call needs to be written 1/1.
constexpr int kUnionHomMinConfidentReads = 5;
// ... and the share of all its reads that are confident.
constexpr double kUnionHomMinConfidentShare = 0.9;
// ... or this share when no read carries REF.
constexpr double kUnionHomMinConfidentShareNoRef = 0.5;
// None of that at full depth (kUnionHomFullDepthRatio of the median or more) with at most one read on REF: both
// haplotypes' reads are there, so the call is no het with one haplotype's reads mapped elsewhere (what the rule
// guards against), however low their MAPQ. Giraffe maps reads at MAPQ 0-29 across whole regions minimap2 maps at
// 60 (chr17:0.93 Mb, 20.4 Mb, 44.1 Mb: KANSL1), and there the rule held back every hom -- inside the benchmark of
// chr17 and chr19, 570 truth homs against 1 het and 4 absent, most of the SNPs DeepVariant called and pgphase did not.
constexpr double kUnionHomFullDepthRatio = 0.75;
// None of it either, at up to kUnionMaxDepthRatio of the median, where the alignment solve's two haplotypes each carry
// the ALT on kUnionHomBothHaplotypesMinReads reads and at most kUnionHomBothHaplotypesMaxRef labelled reads carry REF
// (alt_on_both_haplotypes): both haplotypes are seen to carry it, whatever the depth or the stray REF reads. Inside the
// benchmark of chr12/16/17/18/19/20/21, as an offline ledger counted them, 165 truth homs gained against 9 false loci
// and 14 hets written 1/1; outside it 417 against 88 and 171, no chromosome losing.
constexpr int kUnionHomBothHaplotypesMinReads = 3;
constexpr int kUnionHomBothHaplotypesMaxRef = 1;
// The longest voting row an alignment call can be another description of. A
// call inside a longer row -- a snarl spanning a large SV, up to 1 Mb on chr20
// -- is a variant the graph does not describe at that scale, i.e.
// sample-specific signal, and is injected as a site of its own (VotingRows).
constexpr hts_pos_t kUnionMaxComplementRowSpan = 2000;

// An MSA het the alignment solve left unphased -- its reads too noisy to phase it, a long homopolymer read one
// short as often as not (chr12:953,298, C>CA in A24, the reads' calls 30 REF to 12 ALT where the clusters are 23 and
// 18) -- still lies on that solve's haplotypes: its two consensus sequences are two clusters of reads, and those
// reads' labels in the solve say which haplotype of which phase set each is. Its phase set and haplotype-1 allele
// from them (aln_phase_set, aln_hap1_allele), for phase_by_alignment_solve; nothing when the call is not one allele
// per consensus, or the clusters do not lie on the two haplotypes of one phase set, kMsaClusterMinReads reads each,
// at most kMsaClusterMaxStray of them on the other.
constexpr int kMsaClusterMinReads = 3;
constexpr double kMsaClusterMaxStray = 0.2;
void aln_phase_from_msa_clusters(const PhasingChunk& bam, CandidateVariant& cand) {
    const std::array<int8_t, 2>& allele = cand.msa_cons_allele;
    if (allele[0] < 0 || allele[1] < 0 || allele[0] + allele[1] != 1) return;
    for (const PhasingChunk::MsaRegion& region : bam.msa_regions) {
        if (!region.het() || region.beg >= cand.key.pos || cand.key.pos >= region.end) continue;
        // Per phase set, each cluster's reads on haplotype 1 and 2.
        std::map<hts_pos_t, std::array<std::array<int, 3>, 2>> on;
        for (size_t c = 0; c < 2; ++c)
            for (const int rid : region.cluster_reads[c]) {
                const size_t r = static_cast<size_t>(rid);
                if (rid < 0 || r >= bam.haps.size() || r >= bam.phase_sets.size()) continue;
                const int hap = bam.haps[r];
                if (bam.phase_sets[r] > 0 && (hap == 1 || hap == 2)) ++on[bam.phase_sets[r]][c][static_cast<size_t>(hap)];
            }
        hts_pos_t ps = 0;
        int most = 0;
        for (const auto& [set, n] : on) {
            const int total = n[0][1] + n[0][2] + n[1][1] + n[1][2];
            if (total > most) most = total, ps = set;
        }
        if (ps == 0) return;
        const auto& n = on[ps];
        const auto lies_on = [&n](size_t c, size_t hap) {
            return n[c][hap] >= kMsaClusterMinReads && n[c][3 - hap] <= kMsaClusterMaxStray * (n[c][1] + n[c][2]);
        };
        const bool first_on_hap1 = lies_on(0, 1) && lies_on(1, 2);
        if (!first_on_hap1 && !(lies_on(0, 2) && lies_on(1, 1))) return;
        cand.aln_phase_set = ps;
        cand.aln_hap1_allele = allele[first_on_hap1 ? 0 : 1];
        return;
    }
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
    // A row's span is the reference its catalog allele replaces, as its record is written (trimmed to minimal VCF
    // form), where that is longer than its key's: a walk the graph keys as a one-base insertion can replace hundreds of
    // reference bases (chr20:9,231,470, I:1 over 921 bp), and keyed so it was neither replaced by the MSA calls inside
    // it nor hidden behind them -- written with them, the locus could not be spelled. The allele's own span, not its
    // site's REF: a 13-allele site at chr20:23,847,472 has a 753 bp REF, and the row's allele replaces 702 bp of it.
    explicit VotingRows(const GraphChunkBuildResult& graph_chunk) {
        const PhasingChunk& graph = graph_chunk.chunk;
        for (size_t ci = 0; ci < graph.candidates.size(); ++ci) {
            const CandidateVariant& row = graph.candidates[ci];
            if ((row.lcd_var_i_to_cate & kCandGermlineVarCate) == 0) continue;
            hts_pos_t beg = row.key.pos - 1, end = row.key.pos + std::max(1, row.key.ref_len);
            // Its interior only (a base off each end): a call beside the allele, an insertion at its last base, is not
            // another description of it (chr20:23,848,189 C>CTGTT after that 702 bp walk, hidden as its complement
            // and then written by nothing when the walk lost a conflict).
            if (const std::optional<VcfAllele> a = row_vcf_allele(graph_chunk, ci); a && a->ref.size() > 2) {
                beg = std::min(beg, a->pos);
                end = std::max(end, a->pos - 2 + static_cast<hts_pos_t>(a->ref.size()));
            }
            if (end - beg > kUnionMaxComplementRowSpan) continue;
            spans_.emplace_back(beg, end, ci);
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

// Reads below this MAPQ are counted apart for the trace: a satellite's reads map at 0 or so, a region giraffe maps
// ambiguously (a copy elsewhere in the graph) anywhere below kConfidentReadMinMapq.
constexpr int kTraceNearZeroMapq = 10;

// Each alignment candidate's read support: reads per allele, of them reads
// mapped at kConfidentReadMinMapq or better, depth and the chunk's median depth.
struct AlignmentSupport {
    std::vector<std::map<int, int>> observed;
    std::vector<std::map<int, int>> confident;
    // Reads mapped below kTraceNearZeroMapq (trace only).
    std::vector<int> near_zero;
    // Per phase set, the reads the alignment solve labelled haplotype 1 or 2 carrying REF or the first ALT:
    // {h1 REF, h1 ALT, h2 REF, h2 ALT}.
    std::vector<std::map<hts_pos_t, std::array<int, 4>>> by_hap;
    // Diagnostic counts: reads the alignment solve labelled, and of those, the
    // ones carrying the allele the solve assigned to their haplotype.
    std::vector<std::array<int, 2>> labelled;
    std::vector<int> depth;
    double median_depth = 0.0;

    explicit AlignmentSupport(const PhasingChunk& bam)
        : observed(bam.candidates.size()), confident(bam.candidates.size()), near_zero(bam.candidates.size(), 0),
          by_hap(bam.candidates.size()),
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
                if (read.mapq >= kConfidentReadMinMapq) ++confident[ci][a];
                if (read.mapq < kTraceNearZeroMapq) ++near_zero[ci];
                if ((hap == 1 || hap == 2) && a <= 1 && static_cast<size_t>(prof.read_id) < bam.phase_sets.size() &&
                    bam.phase_sets[static_cast<size_t>(prof.read_id)] > 0)
                    ++by_hap[ci][bam.phase_sets[static_cast<size_t>(prof.read_id)]][static_cast<size_t>(2 * (hap - 1) + a)];
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

    // How far the call's alleles follow the haplotypes the alignment solve labelled its reads with, in the phase set
    // labelling most of them: {reads agreeing with the better of the two allele-to-haplotype pairings, reads labelled}.
    // A het's reads agree but for errors; a collapsed paralog's carry both copies' letters on each haplotype.
    std::array<int, 2> hap_concordance(size_t bi) const {
        std::array<int, 2> best{0, 0};
        for (const auto& [ps, n] : by_hap[bi]) {
            const int labelled_here = n[0] + n[1] + n[2] + n[3];
            if (labelled_here > best[1]) best = {std::max(n[0] + n[3], n[1] + n[2]), labelled_here};
        }
        return best;
    }

    // The call's reads per haplotype in the phase set labelling most of them: {h1 REF, h1 ALT, h2 REF, h2 ALT}.
    std::array<int, 4> on_haplotypes(size_t bi) const {
        std::array<int, 4> best{0, 0, 0, 0};
        int most = 0;
        for (const auto& [ps, n] : by_hap[bi])
            if (n[0] + n[1] + n[2] + n[3] > most) most = n[0] + n[1] + n[2] + n[3], best = n;
        return best;
    }

    // Whether the call's REF and ALT lie one on each haplotype (on_haplotypes): kIslandSplitMinReads reads each side,
    // at most kIslandSplitMaxStray of either haplotype's reads on the other allele.
    bool alleles_split_by_haplotype(size_t bi) const {
        const std::array<int, 4> n = on_haplotypes(bi);
        const auto split = [](int own1, int stray1, int own2, int stray2) {
            return own1 >= kIslandSplitMinReads && own2 >= kIslandSplitMinReads &&
                   stray1 <= kIslandSplitMaxStray * (own1 + stray1) && stray2 <= kIslandSplitMaxStray * (own2 + stray2);
        };
        return split(n[0], n[1], n[3], n[2]) || split(n[1], n[0], n[2], n[3]);
    }

    // Whether both haplotypes carry the ALT (kUnionHomBothHaplotypesMinReads reads each) and at most
    // kUnionHomBothHaplotypesMaxRef labelled reads REF.
    bool alt_on_both_haplotypes(size_t bi) const {
        const std::array<int, 4> n = on_haplotypes(bi);
        return n[1] >= kUnionHomBothHaplotypesMinReads && n[3] >= kUnionHomBothHaplotypesMinReads &&
               n[0] + n[2] <= kUnionHomBothHaplotypesMaxRef;
    }

    // Share of the call's reads on its first ALT.
    double alt_share(size_t bi) const {
        const auto it = observed[bi].find(1);
        return depth[bi] > 0 && it != observed[bi].end() ? static_cast<double>(it->second) / depth[bi] : 0.0;
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

// Whether an MSA call lies inside an MSA deletion of the same region carried by its other consensus only.
static bool under_other_consensus_deletion(const PhasingChunk& bam, const CandidateVariant& c) {
    const std::array<int8_t, 2>& own = c.msa_cons_allele;
    if (own[0] < 0 || own[1] < 0 || own[0] + own[1] != 1) return false;
    const hts_pos_t c0 = c.key.pos, c1 = c.key.pos + std::max(1, c.key.ref_len);
    for (const CandidateVariant& d : bam.candidates) {
        if (d.key.type != VariantType::Deletion || !d.msa_verified || &d == &c) continue;
        const std::array<int8_t, 2>& other = d.msa_cons_allele;
        if (other[0] < 0 || other[1] < 0 || other[0] == own[0] || other[1] == own[1]) continue;
        // The deletion removes bases pos .. pos + ref_len - 1 (VCF: after the padding base at pos - 1).
        if (d.key.pos <= c0 && c1 <= d.key.pos + d.key.ref_len) return true;
    }
    return false;
}

// Whether an alignment candidate lies strictly inside a region the MSA resolved into two different consensus
// sequences (MsaRegion::het).
static bool in_het_msa_region(const PhasingChunk& bam, const CandidateVariant& c) {
    const hts_pos_t lo = c.key.pos - 1;
    const hts_pos_t hi = c.key.pos + std::max(1, c.key.ref_len);
    for (const PhasingChunk::MsaRegion& region : bam.msa_regions)
        if (region.het() && region.beg < lo && hi < region.end) return true;
    return false;
}

// The pattern an alignment candidate falls in, by what the solve sees of it: its category, its depth over the chunk's
// median, its REF and ALT reads and their MAPQ, and how its reads lie on the alignment solve's haplotypes. Each pattern
// has one handling in select_sites; a call none of the measured patterns admits is held back (the default). The trace
// records each call's pattern (inject/decision pattern=), and pgphase-dev docs/patterns.md lists them with their
// evidence and the errors left in each.
enum class CallPattern : uint8_t {
    // Homozygous calls.
    HomInHetMsaRegion,  // an MSA hom inside a region whose consensus sequences differ: written (they spell it)
    HomConfident,       // kUnionHomMinConfidentReads confident ALT reads, most reads confident: written
    HomFullDepth,       // full depth, at most one read on REF, any MAPQ: written
    HomBothHaplotypes,  // ALT on both of the solve's haplotypes, up to kUnionMaxDepthRatio: written
    HomLowMapq,         // none of those (a het whose other haplotype's reads went elsewhere looks so): held back
    // Calls that are no het the union admits.
    NotHet,        // neither a clean het nor an MSA-verified one
    NotGermline,   // the alignment's categories exclude it
    OutsideChunk,  // beyond the graph chunk's reference span
    // Heterozygous calls.
    Het,                      // a clean or confidently read MSA het at ordinary depth: admitted
    CleanHetLowMapq,          // a clean het read at low MAPQ: admitted, flagged low_mapq_only
    GatedOverRows,            // a gated het a catalog row describes: held back (the row stands)
    MsaHetLowMapq,            // an MSA het without confident reads on both alleles: island, not written
    MsaHetLowMapqSplit,       // ... at full depth with its alleles one per haplotype: island, written unphased
    MsaHetDepthGated,         // an MSA het over kUnionMaxDepthRatio: island, written unphased under 2x
    MsaHetCollapsed,          // ... at kIslandWriteMaxDepthRatio or deeper (collapsed copies): island, not written
    CleanHetDepthGatedSplit,  // a clean het over kUnionMaxDepthRatio, under 2x, alleles one per haplotype: written
    CleanHetDepthGated,       // a clean het over kUnionMaxDepthRatio otherwise: held back
};

const char* call_pattern_name(CallPattern p) {
    switch (p) {
    case CallPattern::HomInHetMsaRegion: return "hom_in_het_msa_region";
    case CallPattern::HomConfident: return "hom_confident";
    case CallPattern::HomFullDepth: return "hom_full_depth";
    case CallPattern::HomBothHaplotypes: return "hom_both_haplotypes";
    case CallPattern::HomLowMapq: return "hom_low_mapq";
    case CallPattern::NotHet: return "not_het";
    case CallPattern::NotGermline: return "not_germline";
    case CallPattern::OutsideChunk: return "outside_chunk";
    case CallPattern::Het: return "het";
    case CallPattern::CleanHetLowMapq: return "clean_het_low_mapq";
    case CallPattern::GatedOverRows: return "gated_over_rows";
    case CallPattern::MsaHetLowMapq: return "msa_het_low_mapq";
    case CallPattern::MsaHetLowMapqSplit: return "msa_het_low_mapq_split";
    case CallPattern::MsaHetDepthGated: return "msa_het_depth_gated";
    case CallPattern::MsaHetCollapsed: return "msa_het_collapsed";
    case CallPattern::CleanHetDepthGatedSplit: return "clean_het_depth_gated_split";
    case CallPattern::CleanHetDepthGated: return "clean_het_depth_gated";
    }
    return "?";
}

// What classify_call reads of a candidate.
struct CallEvidence {
    bool hom_in_het_region = false;  // an MSA hom inside a het MSA region, in the chunk
    bool hom_call = false;           // a hom call (clean, MSA or spelled) in the chunk, outside het MSA regions
    bool hom_confident = false;
    bool hom_full_depth = false;
    bool hom_both_haplotypes = false;
    bool clean = false;
    bool msa = false;
    bool germline = false;
    bool in_chunk = false;
    bool confident = false;      // a het's two alleles on confident reads (or the under-deletion exception)
    bool few_confident = false;  // an MSA het without them
    bool depth_gated = false;    // over kUnionMaxDepthRatio of the median
    double depth_ratio = 0.0;
    bool split = false;  // ALT share within kIslandMinAltShare-kIslandMaxAltShare, alleles one per haplotype
    std::vector<size_t> rows;  // voting rows describing the locus (hets only)
};

CallEvidence call_evidence(const PhasingChunk& graph, const PhasingChunk& bam, const AlignmentSupport& support,
                           const VotingRows& voting, const Selection& sel, size_t bi) {
    CallEvidence e;
    const CandidateVariant& c = bam.candidates[bi];
    const VariantCategory cat = c.counts.category;
    e.clean = cat == VariantCategory::CleanHetSnp || cat == VariantCategory::CleanHetIndel;
    e.msa = cat == VariantCategory::NoisyCandHet && c.msa_verified;
    e.germline = (c.lcd_var_i_to_cate & kCandGermlineVarCate) != 0;
    e.in_chunk = c.key.pos >= graph.ref_beg && c.key.pos <= graph.ref_end;
    e.depth_ratio = support.median_depth > 0.0 ? support.depth[bi] / support.median_depth : 0.0;
    const bool het_region = in_het_msa_region(bam, c);
    e.hom_in_het_region = cat == VariantCategory::NoisyCandHom && c.msa_verified && het_region && e.in_chunk;
    const bool hom_call = (cat == VariantCategory::CleanHom && !c.msa_verified) ||
                          (cat == VariantCategory::NoisyCandHom && c.msa_verified) || sel.spelled_hom.count(bi) != 0;
    e.hom_call = !e.hom_in_het_region && hom_call && e.in_chunk && !het_region;
    if (e.hom_call) {
        // Mostly confident reads too: where a share of them map at low MAPQ (satellites), a het with one
        // haplotype's reads elsewhere reads as a hom -- chr20, hom calls under 90% confident reads: the truth a het
        // at 231, a hom at 266 (225 and 153 in satellites); at 90% or more, a het at 20, a hom at 24,224.
        // With no read on REF, half: in segmental duplications reads map at low MAPQ around a true hom.
        // chr17/19/20 (with the depth gate dropped below): SNP F1 0.9899 -> 0.9921 (hap.py; 1,047 fewer missed,
        // 2 more false), INDEL 0.9815 -> 0.9826; outside the benchmark, satellites mostly, 570 truth hets more
        // written 1/1. No depth gate: it keeps paralog runs out of the solve, and a hom call takes no part in it.
        // Held back by it, the truth's hom at 659 of 664 such calls on chr12/17/19/20 (no het).
        const auto conf = support.confident[bi].find(1);
        int n_reads = 0, n_confident = 0;
        for (const auto& [a, n] : support.observed[bi]) n_reads += n;
        for (const auto& [a, n] : support.confident[bi]) n_confident += n;
        const auto ref = support.observed[bi].find(0);
        const bool no_ref = ref == support.observed[bi].end() || ref->second == 0;
        const double share = no_ref ? kUnionHomMinConfidentShareNoRef : kUnionHomMinConfidentShare;
        e.hom_confident = conf != support.confident[bi].end() && conf->second >= kUnionHomMinConfidentReads &&
                          n_confident >= share * n_reads;
        e.hom_full_depth = (ref == support.observed[bi].end() || ref->second <= 1) && support.median_depth > 0.0 &&
                           n_reads >= kUnionHomFullDepthRatio * support.median_depth;
        e.hom_both_haplotypes = support.median_depth > 0.0 && n_reads <= kUnionMaxDepthRatio * support.median_depth &&
                                support.alt_on_both_haplotypes(bi);
        return e;
    }
    if (!(e.clean || e.msa) || !e.germline || !e.in_chunk) return e;
    // An MSA call inside an MSA deletion its other consensus carries: that haplotype has no bases there, so REF has
    // no reads to carry it (a few error reads at most) and the call's own allele alone needs confident support.
    // chr20:49,575,972 and 49,576,250, deletions of 179 and 196 bp beside the other haplotype's 2,443 bp deletion
    // (the truth's 2|1 with *): REF 3 reads, 1 of them confident, and both went unwritten.
    e.confident = support.confidently_supported(bi) ||
                  (e.msa && under_other_consensus_deletion(bam, c) &&
                   support.confident_on_leading(bi)[0] >= kUnionMinConfidentAlleleReads);
    e.few_confident = !e.clean && !e.confident;
    e.depth_gated = support.median_depth > 0.0 && support.depth[bi] > kUnionMaxDepthRatio * support.median_depth;
    const double share = support.alt_share(bi);
    e.split = share >= kIslandMinAltShare && share <= kIslandMaxAltShare && support.alleles_split_by_haplotype(bi);
    e.rows = voting.over(c.key.pos - 1, c.key.pos + std::max(1, c.key.ref_len));
    // Two SNPs overlap only at a shared base: a call beside a row SNP is not another description of it, and hidden
    // as its complement it was never written (chr20:7,706,349 C>T, 39 against 29 reads, beside the catalog SNP at
    // 7,706,350; both truth hets). Only for a call both of whose alleles confident reads carry: in satellites the
    // calls beside catalog SNPs are mostly low-MAPQ paralog reads, and written they unsettled the EM there (chr20
    // centromere -35 truth hets right, +161 genotypes wrong). An indel keeps the one-base reach: an MSA deletion
    // whose anchor letter changes (written AA>G, chr20:51,699,469) spells the catalog SNP at its anchor.
    if (c.key.type == VariantType::Snp && support.confidently_supported(bi))
        e.rows.erase(std::remove_if(e.rows.begin(), e.rows.end(), [&](size_t r) {
                         const VariantKey& k = graph.candidates[r].key;
                         return k.type == VariantType::Snp &&
                                (k.pos + k.ref_len <= c.key.pos || c.key.pos + c.key.ref_len <= k.pos);
                     }), e.rows.end());
    return e;
}

// A call's pattern from its evidence alone.
// An MSA-verified het the confident-read rule or the depth gate drops enters as an island-only site (island_only): its
// reads may phase locally, but it never votes or joins. Admitted as an ordinary site instead, the same calls cost
// -1,258 correct / +1,269 discordant reads (chr12, depth gate: a paralog split joined a block the wrong way) or doubled
// chr20's arm switches (confident-read rule). The depth gate keeps a call out of the solve, not out of the VCF: written
// unphased, as the MSA genotyped it, the MSA's pair was the truth's in 125 of 151 such chr20 regions (chr20:3,414,746,
// CT>C in T16, clusters 47 / 42). Not the few-confident-reads gate's: in satellites its MSA pairs are mostly another
// copy's (147 of 205 wrong) -- but for those at full depth with their alleles one per haplotype. Nor one at
// kIslandWriteMaxDepthRatio or deeper. A clean het is exempt from the confident-read rule: the depth gate alone screens
// it. Where a region has a second copy elsewhere in the reference every read loses its MAPQ, yet the reads are this
// copy's (depth at or below the median, not a collapse): at chr18:21.09-21.10 Mb six truth hets read only at MAPQ < 30
// were dropped and 25 reads LongPhase phases stayed unphased. Such a site enters flagged low_mapq_only: phased, but it
// joins no blocks and labels reads only as an unverified call does (union EM). Measured on chr12/18/20 with the echo
// sites, unverified fallback and per-pair stitching: +155 correct / -88 discordant / -67 unphased, +1,569 truth hets
// phased; arm switches 32 -> 37 (sparse pericentromeric, satellite and VNTR blocks: chr12 34.65-37.19 Mb, chr20 32.76
// and 65.51 Mb).
CallPattern classify_call(const CallEvidence& e) {
    if (e.hom_in_het_region) return CallPattern::HomInHetMsaRegion;
    if (e.hom_call) {
        if (e.hom_confident) return CallPattern::HomConfident;
        if (e.hom_full_depth) return CallPattern::HomFullDepth;
        if (e.hom_both_haplotypes) return CallPattern::HomBothHaplotypes;
        return CallPattern::HomLowMapq;
    }
    if (!e.clean && !e.msa) return CallPattern::NotHet;
    if (!e.germline) return CallPattern::NotGermline;
    if (!e.in_chunk) return CallPattern::OutsideChunk;
    if (!e.few_confident && !e.depth_gated) return e.confident ? CallPattern::Het : CallPattern::CleanHetLowMapq;
    if (!e.rows.empty()) return CallPattern::GatedOverRows;
    if (e.few_confident)
        return e.depth_ratio >= kUnionHomFullDepthRatio && e.depth_ratio <= kUnionMaxDepthRatio && e.split
                   ? CallPattern::MsaHetLowMapqSplit
                   : CallPattern::MsaHetLowMapq;
    if (e.msa)
        return e.depth_ratio < kIslandWriteMaxDepthRatio ? CallPattern::MsaHetDepthGated : CallPattern::MsaHetCollapsed;
    return e.depth_ratio < kIslandWriteMaxDepthRatio && e.split ? CallPattern::CleanHetDepthGatedSplit
                                                                : CallPattern::CleanHetDepthGated;
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
    std::vector<CallPattern> pattern(bam.candidates.size(), CallPattern::NotHet);
    const auto decide = [&](size_t bi, const char* why) {
        if (record_decisions) sel.decision[bi] = why;
        if (!pgphase_trace::on()) return;
        std::string observed;
        for (const auto& [a, n] : support.observed[bi]) {
            const auto it = support.confident[bi].find(a);
            observed += (observed.empty() ? " reads=" : ",") + std::to_string(a) + ":" + std::to_string(n) + "/" +
                        std::to_string(it == support.confident[bi].end() ? 0 : it->second);
        }
        // The call's depth against the chunk's median, its reads mapped near MAPQ 0, and an MSA call's region and
        // cluster reads: what tells a collapsed copy, a satellite and a region giraffe maps ambiguously apart.
        const CandidateVariant& c = bam.candidates[bi];
        char features[192];
        const std::array<int, 2> conc = support.hap_concordance(bi);
        const std::array<int, 4> on = support.on_haplotypes(bi);
        std::snprintf(features, sizeof features,
                      " depth=%d median=%.0f mq_lt%d=%d hapconc=%d/%d haps=%d,%d/%d,%d", support.depth[bi],
                      support.median_depth, kTraceNearZeroMapq, support.near_zero[bi], conc[0], conc[1], on[0], on[1],
                      on[2], on[3]);
        std::string msa_reads;
        if (c.msa_verified)
            msa_reads = " msa_region=" + std::to_string(c.msa_region_reads) + " msa_clustered=" +
                        std::to_string(c.msa_cluster_reads);
        pgphase_trace::candidate("inject", "decision", c, why,
                                 "cat=" + std::to_string(static_cast<int>(c.counts.category)) +
                                 (c.msa_verified ? " msa" : "") + observed + features + msa_reads +
                                 " pattern=" + call_pattern_name(pattern[bi]));
    };
    for (size_t bi = 0; bi < bam.candidates.size(); ++bi) {
        const CandidateVariant& c = bam.candidates[bi];
        CallEvidence e = call_evidence(graph, bam, support, voting, sel, bi);
        const CallPattern p = pattern[bi] = classify_call(e);
        switch (p) {
        case CallPattern::HomInHetMsaRegion:
            decide(bi, "hom");
            sel.hom_records.insert(bi);
            sel.selected.push_back(bi);
            continue;
        case CallPattern::HomLowMapq:
            decide(bi, "hom:few_confident_reads");
            continue;
        case CallPattern::HomConfident:
        case CallPattern::HomFullDepth:
        case CallPattern::HomBothHaplotypes: {
            // The alignment's own homozygous calls, written 1/1 as the MSA's are (no vote, no read calls), and the
            // MSA's in a region whose two consensus sequences agree (the het regions' are above; the alignment's call
            // there is pruned for the region's noise and re-called by the MSA: chr20:148,906 G>A, 62 of 64 reads). The
            // graph drops its hom sites (pair_high_af), so the VCF had no hom calls but those of het MSA regions: on
            // chr20's arms 685 of the truth's 25,031 hom SNPs in the benchmark, SNP F1 0.80 (hap.py) against
            // DeepVariant's 0.99 at more truth het SNPs called; written, 24,850 and SNP F1 0.999, INDEL 0.83 -> 0.99.
            // Where catalog rows describe the locus, their call stands if they write one: the hom is written only
            // where no record with an ALT genotype overlaps it, which the EM decides. Handed to the rows outright,
            // a hom their EM leaves out (one allele: minor_allele_below_fraction) or the MSA drops went unwritten --
            // chr20 arms, 12 truth homs on each BAM (319,724 T>C, 34 of 35 reads).
            const bool fallback = !voting.over(c.key.pos - 1, c.key.pos + std::max(1, c.key.ref_len)).empty();
            // Not an MSA's hom from a region whose clusters took under kUnionFallbackHomMinClusteredShare of its
            // reads: the other haplotype's may not span it (a longer repeat), and its consensus is one haplotype's.
            // chr12:130,683,984 T>A, 24 of 33 reads clustered, all ALT, the truth a het in a VNTR. Over chr12/16/17/
            // 19/20, such fallback homs: none the truth's hom; at 90% or more, 160 right, 10 wrong.
            if (fallback && c.msa_verified && c.msa_region_reads > 0 &&
                c.msa_cluster_reads < kUnionFallbackHomMinClusteredShare * c.msa_region_reads) {
                decide(bi, "hom:graph_row");
                continue;
            }
            const bool spelled_hom = sel.spelled_hom.count(bi) != 0;
            decide(bi, fallback ? "aln_hom:fallback" : spelled_hom ? "aln_hom:spelled" : "aln_hom");
            if (fallback) sel.hom_fallback.insert(bi);
            sel.hom_records.insert(bi);
            sel.selected.push_back(bi);
            continue;
        }
        case CallPattern::NotHet:
            decide(bi, "category");
            continue;
        case CallPattern::NotGermline:
            decide(bi, "not_germline");
            continue;
        case CallPattern::OutsideChunk:
            decide(bi, "outside_chunk");
            continue;
        case CallPattern::GatedOverRows:
            decide(bi, e.few_confident ? "few_confident_reads" : "depth_gate");
            continue;
        case CallPattern::CleanHetDepthGated:
            decide(bi, "depth_gate");
            continue;
        case CallPattern::MsaHetLowMapq:
        case CallPattern::MsaHetLowMapqSplit:
        case CallPattern::MsaHetDepthGated:
        case CallPattern::MsaHetCollapsed:
        case CallPattern::CleanHetDepthGatedSplit:
            decide(bi, e.few_confident ? "island:few_confident_reads" : "island:depth_gate");
            sel.island_only.insert(bi);
            if (p == CallPattern::MsaHetLowMapqSplit || p == CallPattern::MsaHetDepthGated ||
                p == CallPattern::CleanHetDepthGatedSplit)
                sel.island_written.insert(bi);
            sel.selected.push_back(bi);
            continue;
        case CallPattern::Het:
        case CallPattern::CleanHetLowMapq:
            break;
        }
        if (e.msa && !e.rows.empty() && graph_meta.size() == graph.candidates.size()) {
            // An MSA het over catalog rows is written itself, the rows (still voting) not: the MSA's haplotypes are
            // the region's reads' own, the rows' alleles what the graph's reads were forced onto (chr20:509,755, a
            // 4 bp deletion the catalog describes only as a SNP beside it) or genotyped from them (chr12:6,797,743, the MSA's het T>A, 33 against 31 reads, under a catalog row
            // the graph's reads called 14 against 50 and wrote 1/1). Hiding a call whose allele the rows spelled, as a
            // complement site beside them, left these written wrong: HG002 HiFi chr12/18/20 true het regions
            // 29,870 -> 29,979, wrong 795 -> 681, missed 539 -> 531; arm reads +25 discordant.
            decide(bi, "replaces");
            // Every row: keeping one an MSA call reproduces (by VCF allele) put a co-located copy beside the call and
            // left 162 truth hets on chr12/18/20 unphased.
            sel.replaced_rows.insert(e.rows.begin(), e.rows.end());
            if (!e.confident) sel.low_mapq_only.insert(bi);
            sel.selected.push_back(bi);
            continue;
        }
        if (!e.rows.empty() && sel.additive.count(bi) != 0) {
            decide(bi, "injected:adds_to_row");
            if (!e.confident) sel.low_mapq_only.insert(bi);
            sel.selected.push_back(bi);
            continue;
        }
        decide(bi, e.rows.empty() ? "injected" : "complement");
        if (!e.rows.empty()) sel.complement_rows.emplace(bi, std::move(e.rows));
        if (!e.confident) sel.low_mapq_only.insert(bi);
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
void write_diagnostics(const char* diag_dir, const GraphChunkBuildResult& graph_chunk,
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
        const std::optional<VcfAllele> allele = row_vcf_allele(graph_chunk, ci);
        if (!allele) continue;
        const hts_pos_t pos = allele->pos;
        const std::string& ref = allele->ref;
        const std::string& alt = allele->alt;
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
        const std::optional<VcfAllele> allele = row_vcf_allele(graph_chunk, ci);
        if (!allele) continue;
        const hts_pos_t pos = allele->pos;
        const std::string& ref = allele->ref;
        const std::string& alt = allele->alt;
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

// Clean het alignment calls no read spells the reference at, outside MSA regions: over the call's repeat tract and
// kSpellFlank bases each side, fewer than kLocusMinSide reads spell the reference exactly and no sequence other than
// the call's own reaches kLocusMinMinorFraction of the reads. The call is homozygous; the reads counted as REF carry
// another spelling of the same allele, a run one letter short or the insertion placed past a mismatch, which the
// pileup counts as REF (chr20:39,274,817 G>GC, 13 of 75 reads with the G run one short; chr20:2,930,988 C>CCT, 13 of
// 63 with +T and T>C). Called het at AF 0.8, then dropped by the EM (minor allele below fraction) or written 0/1. An
// MSA region, where one forms, genotypes them; where it forms depends on where the aligner put the indel. A true 1/2
// keeps its second sequence and stays het.
constexpr hts_pos_t kSpellFlank = 5;
void spell_homozygous_calls(const PhasingChunk& graph, const PhasingChunk& bam, Selection& sel) {
    if (bam.ref_seq.empty()) return;
    const ReferenceView ref{bam};
    std::optional<SpanningReadIndex> index;
    const auto upper = [](std::string x) {
        for (char& ch : x) ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
        return x;
    };
    for (size_t bi = 0; bi < bam.candidates.size(); ++bi) {
        const CandidateVariant& c = bam.candidates[bi];
        const VariantCategory cat = c.counts.category;
        if (cat != VariantCategory::CleanHetSnp && cat != VariantCategory::CleanHetIndel) continue;
        if ((c.lcd_var_i_to_cate & kCandGermlineVarCate) == 0) continue;
        if (c.key.pos < graph.ref_beg || c.key.pos > graph.ref_end) continue;
        if (c.key.alt.find_first_not_of("ACGTacgt") != std::string::npos) continue;
        const hts_pos_t c0 = c.key.pos - 1, c1 = c0 + c.key.ref_len;
        bool in_msa = false;
        for (const PhasingChunk::MsaRegion& r : bam.msa_regions) in_msa = in_msa || (r.beg <= c1 + 1 && c0 <= r.end);
        if (in_msa) continue;
        const hts_pos_t w0 = std::max<hts_pos_t>(0, c0 - ref.repeat_left(c0) - kSpellFlank);
        const hts_pos_t w1 = c1 + ref.repeat_right(c1) + kSpellFlank;
        if (w1 - w0 > kLocusMaxLength) continue;
        const std::string hap_ref = ref.slice(w0, w1);
        const std::string hap_call = ref.slice(w0, c0) + upper(c.key.alt) + ref.slice(c1, w1);
        if (hap_ref.find('N') != std::string::npos || hap_call == hap_ref) continue;
        if (!index) index.emplace(bam);
        std::map<std::string, int> spelled;
        int n = 0;
        for (auto& [ri, seq] : index->segments(w0, w1, kAlleleWindowMaxReads)) {
            ++spelled[std::move(seq)];
            ++n;
        }
        const auto count = [&spelled](const std::string& x) {
            const auto it = spelled.find(x);
            return it == spelled.end() ? 0 : it->second;
        };
        int n_other = 0;
        for (const auto& [seq, k] : spelled)
            if (seq != hap_ref && seq != hap_call) n_other = std::max(n_other, k);
        if (count(hap_ref) >= kLocusMinSide || count(hap_call) < kLocusMinSide ||
            n_other >= kLocusMinMinorFraction * n) continue;
        sel.spelled_hom.insert(bi);
        if (pgphase_trace::on())
            pgphase_trace::candidate("inject", "spell", c, "homozygous",
                                     "reads=" + std::to_string(n) + " ref=" + std::to_string(count(hap_ref)) +
                                         " call=" + std::to_string(count(hap_call)) + " other=" + std::to_string(n_other));
    }
}

// Clean alignment hets that add a change to a catalog substitution row they overlap rather than describe it again: the
// reads spell the row's and the call's changes together (kLocusMinSide reads, more than spell the row's alone) over a
// window through the repeat tracts at its ends. Hidden as the row's complement, such a call was written by nothing:
// chr20:62,977,373 C>CAT inside the catalog's CAC>CAT, whose SNP the row writes (the truth's C>CAT and C>T at 375);
// chr20:58,239,873 G>GAT. Written as an injected call, the row still written beside it.
void find_additive_calls(const GraphChunkBuildResult& graph_chunk, const PhasingChunk& bam, const VotingRows& voting,
                         Selection& sel) {
    const PhasingChunk& graph = graph_chunk.chunk;
    if (bam.ref_seq.empty() || graph_chunk.site_meta.size() != graph.candidates.size() ||
        graph_chunk.site_allele_orig_idx.size() != graph.candidates.size()) return;
    const ReferenceView ref{bam};
    std::optional<SpanningReadIndex> index;
    const auto upper = [](std::string x) {
        for (char& ch : x) ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
        return x;
    };
    for (size_t bi = 0; bi < bam.candidates.size(); ++bi) {
        const CandidateVariant& c = bam.candidates[bi];
        const VariantCategory cat = c.counts.category;
        if (cat != VariantCategory::CleanHetSnp && cat != VariantCategory::CleanHetIndel) continue;
        if ((c.lcd_var_i_to_cate & kCandGermlineVarCate) == 0) continue;
        if (c.key.pos < graph.ref_beg || c.key.pos > graph.ref_end) continue;
        if (c.key.alt.find_first_not_of("ACGTacgt") != std::string::npos) continue;
        const hts_pos_t c0 = c.key.pos - 1, c1 = c0 + c.key.ref_len;
        for (const size_t r : voting.over(c.key.pos - 1, c.key.pos + std::max(1, c.key.ref_len))) {
            if (graph.candidates[r].counts.alle_covs.size() != 2) continue;
            const std::optional<VcfAllele> allele = row_vcf_allele(graph_chunk, r);
            if (!allele) continue;
            const std::string& row_ref = allele->ref;
            const std::string& row_alt = allele->alt;
            // A substitution row: its bases one for one.
            if (row_ref.size() != row_alt.size() || row_alt.find_first_not_of("ACGT") != std::string::npos) continue;
            const hts_pos_t r0 = allele->pos - 1, r1 = r0 + static_cast<hts_pos_t>(row_ref.size());
            if (ref.slice(r0, r1) != row_ref) continue;
            const hts_pos_t lo = std::min(r0, c0), hi = std::max(r1, c1);
            const hts_pos_t w0 = std::max<hts_pos_t>(0, lo - ref.repeat_left(lo) - kAlleleWindowFlank);
            const hts_pos_t w1 = hi + ref.repeat_right(hi) + kAlleleWindowFlank;
            if (w1 - w0 > kLocusMaxLength) continue;
            // Slots: each reference base of the window, the row's bases in place, then the call applied to bases the
            // row leaves alone (an insertion before its base).
            std::vector<std::string> slot;
            for (hts_pos_t x = w0; x < w1; ++x) slot.emplace_back(1, ref.base(x));
            std::vector<char> changed(slot.size(), 0);
            for (hts_pos_t x = r0; x < r1; ++x) {
                const size_t k = static_cast<size_t>(x - w0);
                if (slot[k][0] != row_alt[static_cast<size_t>(x - r0)]) changed[k] = 1;
                slot[k] = std::string(1, row_alt[static_cast<size_t>(x - r0)]);
            }
            const auto join = [](const std::vector<std::string>& v) {
                std::string out;
                for (const std::string& x : v) out += x;
                return out;
            };
            const std::string hap_row = join(slot);
            bool clash = false;
            for (hts_pos_t x = c0; x < c1; ++x) clash = clash || changed[static_cast<size_t>(x - w0)];
            if (clash) continue;
            const std::string alt = upper(c.key.alt);
            if (c.key.ref_len == 0) {
                slot[static_cast<size_t>(c0 - w0)] = alt + slot[static_cast<size_t>(c0 - w0)];
            } else {
                for (hts_pos_t x = c0; x < c1; ++x) slot[static_cast<size_t>(x - w0)].clear();
                slot[static_cast<size_t>(c0 - w0)] = alt;
            }
            const std::string hap_both = join(slot);
            if (hap_both == hap_row) continue;
            if (!index) index.emplace(bam);
            int n_both = 0, n_row = 0;
            for (const auto& [ri, seq] : index->segments(w0, w1, kAlleleWindowMaxReads)) {
                n_both += seq == hap_both;
                n_row += seq == hap_row;
            }
            if (n_both >= kLocusMinSide && n_both > n_row) {
                sel.additive.insert(bi);
                break;
            }
        }
    }
}

// Catalog rows a clean alignment het overlaps, outside MSA regions, that no read's aligned sequence spells: as
// unwrite_indel_rows_msa_lacks checks rows against the MSA's consensus sequences, here against the reads. Over a window
// around the row and the call, through the repeat tracts at its ends, the spanning reads' sequences are counted: a row
// whose haplotype fewer than kLocusMinSide reads spell exactly, where another non-reference sequence is spelled by
// kLocusMinSide reads and kLocusMinMinorFraction of them, is no haplotype of the reads (the graph put them on it, at
// the nearest allele it has): written no record, still voting, and the calls hidden behind it as its complement are
// promoted (select_sites). Where the alignment spelled a locus as the truth does and differently from the catalog, the
// catalog's description was written and the reads' own hidden: chr20:55,942,094 CAAAT>C, one A short of the reads'
// haplotype, written 1/0 over the truth's A>G, G>T and AATT>A the alignment called; chr20:1,573,393 A>G over the truth's
// T>TG at 391. Not per read, nearest haplotype: the call is one of several changes on the reads' haplotype there, and
// the row can be nearer the reads than it is. An MSA region resolves its rows itself.
void unwrite_rows_reads_lack(const GraphChunkBuildResult& graph_chunk, const PhasingChunk& bam,
                             const AlignmentSupport& support, const VotingRows& voting, Selection& sel) {
    const PhasingChunk& graph = graph_chunk.chunk;
    if (bam.ref_seq.empty() || graph_chunk.site_meta.size() != graph.candidates.size() ||
        graph_chunk.site_allele_orig_idx.size() != graph.candidates.size()) return;
    const ReferenceView ref{bam};
    std::optional<SpanningReadIndex> index;  // built at the first pair tested
    const auto upper = [](std::string x) {
        for (char& ch : x) ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
        return x;
    };
    const auto is_acgt = [](const std::string& x) {
        return x.find_first_not_of("ACGTacgt") == std::string::npos;
    };
    const auto in_msa_region = [&bam](hts_pos_t a, hts_pos_t b) {
        for (const PhasingChunk::MsaRegion& r : bam.msa_regions)
            if (r.beg <= b && a <= r.end) return true;
        return false;
    };
    for (size_t bi = 0; bi < bam.candidates.size(); ++bi) {
        const CandidateVariant& c = bam.candidates[bi];
        const VariantCategory cat = c.counts.category;
        if (cat != VariantCategory::CleanHetSnp && cat != VariantCategory::CleanHetIndel) continue;
        if ((c.lcd_var_i_to_cate & kCandGermlineVarCate) == 0 || !support.confidently_supported(bi)) continue;
        if (c.key.pos < graph.ref_beg || c.key.pos > graph.ref_end || !is_acgt(c.key.alt)) continue;
        const hts_pos_t c0 = c.key.pos - 1, c1 = c0 + c.key.ref_len;
        if (in_msa_region(c0, c1 + 1)) continue;
        for (const size_t r : voting.over(c.key.pos - 1, c.key.pos + std::max(1, c.key.ref_len))) {
            if (sel.unwritten_rows.count(r) != 0 || graph.candidates[r].counts.alle_covs.size() != 2) continue;
            const VariantKey& rk = graph.candidates[r].key;
            // As select_sites: a SNP beside a row SNP is not another description of it.
            if (c.key.type == VariantType::Snp && rk.type == VariantType::Snp &&
                (rk.pos + rk.ref_len <= c.key.pos || c.key.pos + c.key.ref_len <= rk.pos)) continue;
            const std::optional<VcfAllele> allele = row_vcf_allele(graph_chunk, r);
            if (!allele) continue;
            const std::string& row_ref = allele->ref;
            const std::string& row_alt = allele->alt;
            const hts_pos_t r0 = allele->pos - 1, r1 = r0 + static_cast<hts_pos_t>(row_ref.size());
            if (!is_acgt(row_alt) || ref.slice(r0, r1) != row_ref || in_msa_region(r0, r1)) continue;
            const hts_pos_t lo = std::min(r0, c0), hi = std::max(r1, c1);
            const hts_pos_t w0 = std::max<hts_pos_t>(0, lo - ref.repeat_left(lo) - kAlleleWindowFlank);
            const hts_pos_t w1 = hi + ref.repeat_right(hi) + kAlleleWindowFlank;
            if (w1 - w0 > kLocusMaxLength) continue;
            const std::string hap_ref = ref.slice(w0, w1);
            const std::string hap_row = ref.slice(w0, r0) + row_alt + ref.slice(r1, w1);
            if (hap_ref.find('N') != std::string::npos || hap_row == hap_ref) continue;
            // An alignment candidate (kLocusMinSide reads) describing the row, alone or with the others inside its span,
            // is the alignment's own spelling of it. Checked first: over a window with other variants every read
            // spells some of them too, and none the row's change alone (chr20:47.29 Mb, catalog SNPs a few bases apart
            // the alignment calls too, all unwritten without it).
            bool described = false;
            std::string combined = ref.slice(w0, w1);
            hts_pos_t shift = 0, last_end = -1;
            for (const CandidateVariant& o : bam.candidates) {
                if (described) break;
                const hts_pos_t o0 = o.key.pos - 1, o1 = o0 + o.key.ref_len;
                if (o.counts.alt_cov < kLocusMinSide || o0 < w0 || o1 > w1) continue;
                if (o.key.alt.find_first_not_of("ACGTacgt") != std::string::npos) continue;
                const std::string alt = upper(o.key.alt);
                described = ref.slice(w0, o0) + alt + ref.slice(o1, w1) == hap_row;
                if (o0 >= r0 && o1 <= r1 && o0 >= last_end) {
                    combined.replace(static_cast<size_t>(o0 - w0 + shift), static_cast<size_t>(o.key.ref_len), alt);
                    shift += static_cast<hts_pos_t>(alt.size()) - o.key.ref_len;
                    last_end = o1;
                }
            }
            if (described || combined == hap_row) continue;
            if (!index) index.emplace(bam);
            std::map<std::string, int> spelled;
            int n = 0;
            for (auto& [ri, seq] : index->segments(w0, w1, kAlleleWindowMaxReads)) {
                ++spelled[std::move(seq)];
                ++n;
            }
            const auto row_it = spelled.find(hap_row);
            if (row_it != spelled.end() && row_it->second >= kLocusMinSide) continue;
            int n_other = 0;
            for (const auto& [seq, k] : spelled)
                if (seq != hap_ref) n_other = std::max(n_other, k);
            if (n_other < kLocusMinSide || n_other < kLocusMinMinorFraction * n) continue;
            if (sel.unwritten_rows.insert(r).second) sel.row_reason[r] = "aln_reads_lack_row";
        }
    }
}

// Catalog rows across the edge of a region the MSA resolved, that no read spells: the MSA's row checks
// (suppress_rows_msa_reads_as_reference, genotype_substitution_rows, unwrite_indel_rows_msa_lacks) take rows inside a
// region only, and where the region begins depends on the alignment. Written beside the region's own calls, such a
// row could not be spelled with them: chr20:9,231,470, a 648 bp catalog walk one base off the reads' haplotype (a T at
// 9,231,535 that 2 of 53 reads carry), began 22 bp before the region and was written with its 18 SNPs. The test is
// unwrite_rows_reads_lack's: no alignment candidate describes the row, fewer than kLocusMinSide reads spell it, and
// another non-reference sequence kLocusMinSide reads and kLocusMinMinorFraction of them do. Written no record, voting.
void unwrite_edge_rows_reads_lack(const GraphChunkBuildResult& graph_chunk, const PhasingChunk& bam, Selection& sel) {
    const PhasingChunk& graph = graph_chunk.chunk;
    if (bam.ref_seq.empty() || bam.msa_regions.empty() || graph_chunk.site_meta.size() != graph.candidates.size() ||
        graph_chunk.site_allele_orig_idx.size() != graph.candidates.size()) return;
    const ReferenceView ref{bam};
    std::optional<SpanningReadIndex> index;
    const auto upper = [](std::string x) {
        for (char& ch : x) ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
        return x;
    };
    const auto acgt = [](const std::string& x) { return x.find_first_not_of("ACGTacgt") == std::string::npos; };
    for (size_t ci = 0; ci < graph.candidates.size(); ++ci) {
        const CandidateVariant& row = graph.candidates[ci];
        if ((row.lcd_var_i_to_cate & kCandGermlineVarCate) == 0 || row.counts.alle_covs.size() != 2 ||
            sel.unwritten_rows.count(ci) != 0) continue;
        const std::optional<VcfAllele> allele = row_vcf_allele(graph_chunk, ci);
        if (!allele) continue;
        const std::string& row_ref = allele->ref;
        const std::string& row_alt = allele->alt;
        const hts_pos_t r0 = allele->pos - 1, r1 = r0 + static_cast<hts_pos_t>(row_ref.size());
        bool edge = false;
        for (const PhasingChunk::MsaRegion& region : bam.msa_regions)
            if (region.het() && region.beg <= r1 && r0 <= region.end && (r0 <= region.beg || r1 >= region.end)) edge = true;
        if (!edge || !acgt(row_alt) || ref.slice(r0, r1) != row_ref) continue;
        const hts_pos_t w0 = std::max<hts_pos_t>(0, r0 - ref.repeat_left(r0) - kAlleleWindowFlank);
        const hts_pos_t w1 = r1 + ref.repeat_right(r1) + kAlleleWindowFlank;
        if (w1 - w0 > kLocusMaxLength) continue;
        const std::string hap_ref = ref.slice(w0, w1);
        const std::string hap_row = ref.slice(w0, r0) + row_alt + ref.slice(r1, w1);
        if (hap_ref.find('N') != std::string::npos || hap_row == hap_ref) continue;
        bool described = false;
        std::string combined = hap_ref;
        hts_pos_t shift = 0, last_end = -1;
        for (const CandidateVariant& o : bam.candidates) {
            const hts_pos_t o0 = o.key.pos - 1, o1 = o0 + o.key.ref_len;
            if (o.counts.alt_cov < kLocusMinSide || o0 < w0 || o1 > w1 || !acgt(o.key.alt)) continue;
            const std::string alt = upper(o.key.alt);
            described = described || ref.slice(w0, o0) + alt + ref.slice(o1, w1) == hap_row;
            if (o0 >= r0 && o1 <= r1 && o0 >= last_end) {
                combined.replace(static_cast<size_t>(o0 - w0 + shift), static_cast<size_t>(o.key.ref_len), alt);
                shift += static_cast<hts_pos_t>(alt.size()) - o.key.ref_len;
                last_end = o1;
            }
        }
        if (described || combined == hap_row) continue;
        if (!index) index.emplace(bam);
        std::map<std::string, int> spelled;
        int n = 0;
        for (auto& [ri, seq] : index->segments(w0, w1, kAlleleWindowMaxReads)) {
            ++spelled[std::move(seq)];
            ++n;
        }
        const auto row_it = spelled.find(hap_row);
        if (row_it != spelled.end() && row_it->second >= kLocusMinSide) continue;
        int n_other = 0;
        for (const auto& [seq, k] : spelled)
            if (seq != hap_ref) n_other = std::max(n_other, k);
        if (n_other < kLocusMinSide || n_other < kLocusMinMinorFraction * n) continue;
        if (sel.unwritten_rows.insert(ci).second) sel.row_reason[ci] = "msa_edge_reads_lack_row";
    }
}

}  // namespace

Selection plan_injection(const GraphChunkBuildResult& graph_chunk, const PhasingChunk& bam, bool diagnostics) {
    const PhasingChunk& graph = graph_chunk.chunk;
    const VotingRows voting(graph_chunk);
    const AlignmentSupport support(bam);
    Selection sel;
    suppress_rows_msa_reads_as_reference(graph, bam, sel);
    genotype_substitution_rows(graph_chunk, bam, sel);
    unwrite_indel_rows_msa_lacks(graph_chunk, bam, sel);
    unwrite_rows_reads_lack(graph_chunk, bam, support, voting, sel);
    unwrite_edge_rows_reads_lack(graph_chunk, bam, sel);
    spell_homozygous_calls(graph, bam, sel);
    find_additive_calls(graph_chunk, bam, voting, sel);
    select_sites(graph, bam, support, voting, graph_chunk.site_meta, diagnostics, sel);
    return sel;
}

void write_injection_diagnostics(const char* diag_dir, const GraphChunkBuildResult& graph_chunk,
                                 const PhasingChunk& bam, const Selection& sel) {
    write_diagnostics(diag_dir, graph_chunk, bam, AlignmentSupport(bam), sel);
}

std::vector<char> collapsed_duplication_sites(const PhasingChunk& bam) {
    const AlignmentSupport support(bam);
    std::vector<char> out(bam.candidates.size(), 0);
    for (size_t bi = 0; bi < bam.candidates.size(); ++bi)
        out[bi] = support.median_depth > 0.0 && support.depth[bi] > kUnionMaxDepthRatio * support.median_depth;
    return out;
}

} // namespace pgphase_collect
