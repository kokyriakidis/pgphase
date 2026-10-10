// Shared by the stages of union gap phasing (union_*.cpp): their constants
// and the helpers more than one stage uses. Internal; the API is union_phase.hpp.

#ifndef PGPHASE_UNION_INTERNAL_HPP
#define PGPHASE_UNION_INTERNAL_HPP

#include "union_phase.hpp"

#include "collect_phase.hpp"
#include "collect_var.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <map>
#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <htslib/sam.h>

namespace pgphase_collect {

// Noisy calls closer than this merge into one locus; padding; the longest
// window; the MAPQ a read needs to take part; reads per window; seed sample
// per side; minimum reads per side and minor fraction; the edit-distance margin
// needed to assign a read; repeat scan bounds; the largest per-read error.
constexpr hts_pos_t kLocusMergeGap = 100;
constexpr hts_pos_t kLocusPad = 40;
constexpr hts_pos_t kLocusMaxLength = 2000;
constexpr int kLocusMinMapq = 5;
constexpr size_t kLocusMaxReads = 60;
constexpr size_t kLocusSeedSample = 10;
constexpr int kLocusMinSide = 3;
constexpr double kLocusMinMinorFraction = 0.2;
constexpr int kLocusMinMargin = 2;
constexpr int kLocusSplitRounds = 3;  // refinements of an unseeded two-way split
// Allele windows: flank around the repeat, and reads per window. A read is
// assigned at any edit margin: both sequences are exact alleles, so a one-edit
// difference is the allele itself, not a read's sequencing error.
constexpr hts_pos_t kAlleleWindowFlank = 20;
constexpr size_t kAlleleWindowMaxReads = 100;
// An indel is a candidate for a two-allele pair when at most this fraction of
// the reads at it carry the reference.
constexpr double kPairMaxRefFraction = 0.2;
// Alignment reads missing from the graph chunk are added from this MAPQ:
// confidently placed, so not a paralog the graph put elsewhere.
constexpr int kAddedReadMinMapq = 20;
// Last-resort labels: an MSA-verified indel orients itself from at least this
// many labelled reads of one phase set, at this concordance.
constexpr int kLastResortMinLabelled = 4;
constexpr double kLastResortMinAgreement = 0.9;
// Haplotype-consensus labels: an event marks a haplotype when each side has at
// least this many reads covering it (with this flank) and it is on at least
// this fraction of one side's reads and at most the complement of the other's.
// A read is labelled when this share of its votes agrees.
constexpr int kConsensusMinReads = 3;
constexpr hts_pos_t kConsensusFlank = 10;
constexpr double kConsensusCarrierFraction = 0.8;
constexpr double kConsensusMinAgreement = 0.9;
constexpr int kEventMinBaseQuality = 10;
constexpr int kLocusMaxRepeatScan = 500;
constexpr int kLocusMaxPeriod = 50;
constexpr double kLocusMaxObsError = 0.45;
// A read's label for seeding: decisive agreement with one haplotype's
// consensus at a phase set's phased sites.
constexpr int kSeedLabelMinMargin = 2;
constexpr double kSeedLabelMinAgreement = 0.9;

// Query index aligned at (or, inside a deletion, after) a 0-based reference position.
hts_pos_t locus_query_index(const bam1_t* b, hts_pos_t ref_pos);

int locus_edit_distance(const std::string& a, const std::string& b);

// The alignment chunk's reference in 0-based coordinates (ref_seq is 1-based
// from ref_beg), with tandem-repeat extents for window bounds.
struct ReferenceView {
    const PhasingChunk& chunk;
    char base(hts_pos_t p0) const {
        const hts_pos_t off = p0 + 1 - chunk.ref_beg;
        if (off < 0 || static_cast<size_t>(off) >= chunk.ref_seq.size()) return 'N';
        return static_cast<char>(std::toupper(static_cast<unsigned char>(chunk.ref_seq[static_cast<size_t>(off)])));
    }
    std::string slice(hts_pos_t p0, hts_pos_t p1) const {
        std::string out;
        for (hts_pos_t p = p0; p < p1; ++p) out += base(p);
        return out;
    }
    // Length of the longest tandem repeat starting at `start` (rightward) or
    // ending just before it (leftward).
    hts_pos_t repeat_right(hts_pos_t start) const {
        hts_pos_t best = 0;
        for (int period = 1; period <= kLocusMaxPeriod; ++period) {
            hts_pos_t k = start + period;
            while (k - start < kLocusMaxRepeatScan && base(k) != 'N' && base(k) == base(k - period)) ++k;
            if (k - start >= 2 * period) best = std::max(best, k - start);
        }
        return best;
    }
    hts_pos_t repeat_left(hts_pos_t end) const {
        hts_pos_t best = 0;
        for (int period = 1; period <= kLocusMaxPeriod; ++period) {
            hts_pos_t k = end - 1 - period;
            while (end - 1 - k < kLocusMaxRepeatScan && base(k) != 'N' && base(k) == base(k + period)) --k;
            if (end - 1 - k >= 2 * period) best = std::max(best, end - 1 - k);
        }
        return best;
    }
};

// Primary alignment-chunk reads the graph chunk also holds, by start.
struct SpanningReadIndex {
    std::vector<const ReadRecord*> reads;
    std::vector<hts_pos_t> starts;
    std::vector<size_t> graph_index;

    SpanningReadIndex(const PhasingChunk& bam, const PhasingChunk& graph) {
        std::unordered_map<std::string_view, size_t> graph_read;
        graph_read.reserve(graph.reads.size());
        for (size_t ri = 0; ri < graph.reads.size(); ++ri) graph_read.try_emplace(graph.reads[ri].qname, ri);
        std::vector<std::pair<const ReadRecord*, size_t>> kept;
        for (const ReadRecord& r : bam.reads) {
            if (!r.alignment) continue;
            if ((r.alignment->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FDUP | BAM_FUNMAP)) != 0)
                continue;
            if (r.alignment->core.qual < kLocusMinMapq) continue;
            const auto g = graph_read.find(r.qname);
            if (g == graph_read.end()) continue;
            kept.emplace_back(&r, g->second);
        }
        std::sort(kept.begin(), kept.end(), [](const auto& a, const auto& b) {
            return a.first->alignment->core.pos < b.first->alignment->core.pos;
        });
        for (const auto& [r, gi] : kept) {
            reads.push_back(r);
            starts.push_back(r->alignment->core.pos);
            graph_index.push_back(gi);
        }
    }

    // Reads covering the whole of [l0, l1), cut to it: (graph read, sequence).
    std::vector<std::pair<size_t, std::string>> segments(hts_pos_t l0, hts_pos_t l1, size_t max_reads) const {
        std::vector<std::pair<size_t, std::string>> segs;
        const size_t last = static_cast<size_t>(std::upper_bound(starts.begin(), starts.end(), l0) - starts.begin());
        for (size_t i = 0; i < last && segs.size() < max_reads; ++i) {
            const bam1_t* aln = reads[i]->alignment.get();
            if (bam_endpos(aln) < l1) continue;
            const hts_pos_t q0 = locus_query_index(aln, l0), q1 = locus_query_index(aln, l1 - 1);
            if (q0 < 0 || q1 < q0) continue;
            std::string seq(static_cast<size_t>(q1 - q0 + 1), 'N');
            const uint8_t* packed = bam_get_seq(aln);
            for (hts_pos_t k = 0; k <= q1 - q0; ++k)
                seq[static_cast<size_t>(k)] = seq_nt16_str[bam_seqi(packed, q0 + k)];
            segs.emplace_back(graph_index[i], std::move(seq));
        }
        return segs;
    }
};

// Alleles competing with an indel at its locus: the reference and every other
// called allele whose edit lies inside the same window. At a two-allele
// heterozygote (1/2) the "other" haplotype carries the second allele, not the
// reference; contrasting REF with ALT alone calls its reads ALT.
struct LocusAllele { hts_pos_t p0; hts_pos_t end0; std::string alt; };

std::vector<std::string> competing_haplotypes(const ReferenceView& ref, hts_pos_t w0, hts_pos_t w1,
                                                     const LocusAllele& self,
                                                     const std::vector<LocusAllele>& sorted_alleles);

int closest_allele_call(const std::string& seq, const std::string& hap_alt,
                               const std::vector<std::string>& competitors, int min_margin);

}  // namespace pgphase_collect

#endif  // PGPHASE_UNION_INTERNAL_HPP
