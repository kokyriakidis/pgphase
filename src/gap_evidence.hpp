#ifndef PGPHASE_GAP_EVIDENCE_HPP
#define PGPHASE_GAP_EVIDENCE_HPP

#include "gap_recovery.hpp"

namespace pgphase_collect {

enum class GapEventRole { Anchor, Private, Graph, Boundary, Unsupported };
enum class GapAlleleStatus { Missing, Observed, LowQuality, Conflicting };

struct GapAllele {
    int index = -1;
    GapAlleleStatus status = GapAlleleStatus::Missing;
    int query_index = -1;
    int base_quality = -1;
};

struct GapEvent {
    // Zero-based, half-open replacement interval; allele 0 is reference.
    hts_pos_t beg = 0, end = 0;
    std::vector<std::string> alleles;
    std::vector<int> msa_to_event;
    std::string id;
    GapEventRole role = GapEventRole::Unsupported;
    CandidateVariant prototype;
    hts_pos_t anchor_ps = -1;
    std::array<int, 2> anchor_alleles{-1, -1};
};

struct GapObservation {
    size_t read = 0, event = 0;
    GapAllele bam, graph, msa;
};

/// Immutable evidence for one original graph phase gap. Projection is the only
/// path into mutable k-means state; repeated trials never change this snapshot.
class GapEvidence {
public:
    GapEvidence(PhasingChunk input, const PhaseGap& gap);
    PhasingChunk project(const Options& opts = Options{}, bool bam_only = false) const;
    const std::vector<GapEvent>& events() const;
    const std::vector<GapObservation>& observations() const;
    const PhaseGap& gap() const;
    void write_audit(const std::string& prefix) const;

private:
    PhaseGap gap_;
    RegionChunk region_;
    hts_pos_t ref_beg_ = 0, ref_end_ = 0;
    std::string reference_;
    std::vector<Interval> low_complexity_, noisy_;
    std::vector<ReadRecord> reads_;
    std::vector<GapEvent> events_;
    std::vector<GapObservation> observations_;
};

/// Footprint ownership, excluding the two graph anchor positions.
bool gap_owns_variant(const PhaseGap& gap, const VariantKey& key);

} // namespace pgphase_collect
#endif
