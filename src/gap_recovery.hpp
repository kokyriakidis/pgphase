#ifndef PGPHASE_GAP_RECOVERY_HPP
#define PGPHASE_GAP_RECOVERY_HPP

#include "phasing_types.hpp"

#include <string_view>
#include <unordered_map>

namespace pgphase_collect {

struct PhaseGap {
    int tid;
    hts_pos_t left_ps;
    hts_pos_t right_ps;
    hts_pos_t left_end;
    hts_pos_t right_beg;
    hts_pos_t region_beg;
    hts_pos_t region_end;
};

/// Index read locations and freeze the pre-recovery haplotype/PS assignments.
struct GapReadIndex {
    using Key = std::pair<int, std::string_view>;
    struct Hash { size_t operator()(const Key& key) const; };
    std::unordered_map<Key, std::vector<std::pair<size_t, size_t>>, Hash> reads;
    std::unordered_map<Key, std::pair<int, hts_pos_t>, Hash> assignments;
    explicit GapReadIndex(const std::vector<PhasingChunk>& chunks);
};

struct GapStitchResult {
    bool left_linked = false;
    bool right_linked = false;
    bool joined = false;
    bool right_flip = false;
    int reads_added = 0;
};

struct GapPhaseEdge {
    hts_pos_t left_ps;
    hts_pos_t right_ps;
    bool right_flip;
};

/// Select graph-observed reads intersecting a gap and retain BAM/MSA evidence.
int select_graph_gap_bam_reads(PhasingChunk& proposal, const PhaseGap& gap,
                               const Options& opts);

/// Find internal gaps between non-overlapping, read-supported phase blocks.
std::vector<PhaseGap> find_phase_gaps(const std::vector<PhasingChunk>& chunks);

/// Stitch a local proposal to its flanks using the normal overlap vote rule.
/// Existing blocks are only relabelled uniformly; unresolved reads may be added.
GapStitchResult stitch_gap_proposal(std::vector<PhasingChunk>& chunks,
                                    const PhasingChunk& proposal,
                                    const PhaseGap& gap, const Options& opts,
                                    const GapReadIndex* read_index = nullptr,
                                    bool defer_phase_set_merge = false,
                                    bool orientation_only = false);

/// Resolve accepted gap relationships and relabel all blocks in one pass.
/// Returns the number of parity-conflicting relationships that were rejected.
int apply_gap_phase_edges(std::vector<PhasingChunk>& chunks,
                          const std::vector<GapPhaseEdge>& edges);

/// Add gap-driven MSA windows, preserving detected noisy intervals intact.
void prepare_gap_msa_regions(PhasingChunk& chunk, hts_pos_t beg, hts_pos_t end);

/// Admit SNPs only at tier 2; retain them and admit indels at tier 3.
void run_gap_msa_tier(PhasingChunk& chunk, const Options& opts,
                      hts_pos_t beg, hts_pos_t end, bool snp_only);

} // namespace pgphase_collect
#endif
