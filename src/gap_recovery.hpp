#ifndef PGPHASE_GAP_RECOVERY_HPP
#define PGPHASE_GAP_RECOVERY_HPP

#include "phasing_types.hpp"

#include <set>
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
    // Every phase-set id already carrying a committed pre-recovery read
    // (the set of `assignments` values' second element). Phase-set ids are
    // genome positions, so a gap's own locally re-derived k-means can, by
    // coincidence, reconstruct the same id as a real, pre-existing block
    // (e.g. its local window's k-means anchors on the same first het site as
    // an adjacent established block). emit_independent_gap_block must never
    // emit into one of these -- it has no orientation-vote safety net the
    // way stitch_gap_proposal does, so silently landing gap-only reads in an
    // already-real phase set pollutes that block with unvalidated orientation.
    std::set<hts_pos_t> established_phase_sets;
    explicit GapReadIndex(const std::vector<PhasingChunk>& chunks);
};

struct GapLinkEvidence {
    hts_pos_t proposal_ps;
    std::array<std::array<int, 4>, 2> votes;
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
                                    bool orientation_only = false,
                                    std::vector<GapLinkEvidence>* evidence = nullptr);

/// Resolve accepted gap relationships and relabel all blocks in one pass.
/// Returns the number of parity-conflicting relationships that were rejected.
int apply_gap_phase_edges(std::vector<PhasingChunk>& chunks,
                          const std::vector<GapPhaseEdge>& edges);

/// Emit a brand-new, independent phase block from a gap's own local solve,
/// for reads that carry no original (pre-recovery) haplotype/phase-set
/// assignment at all -- i.e. reads whose only evidence lives inside the gap
/// itself, on neither known flank. `stitch_gap_proposal` only ever extends an
/// EXISTING flank when it can confidently vote for one; when a gap cannot
/// vote for either flank (both sides thin or absent) its proposal is
/// otherwise discarded even if the gap-only reads formed a perfectly
/// coherent block among themselves. This is checked separately, using the
/// gap's own local phase-set id and requires at least `min_reads` gap-only
/// reads sharing it -- a deliberately additive, minimal-risk operation: it
/// can only assign a read that previously had no phase at all, and never
/// touches, overwrites, or reorients a read or candidate the flanks already
/// own. Phase-set ids are genome positions, NOT unique by construction: a
/// gap's local re-solve can coincidentally reproduce the id of an already-
/// established block (`read_index.established_phase_sets`) or of another
/// gap emitted earlier in the same sequential batch (`emitted_this_round`,
/// in/out -- pass the same set across every call in one batch and this
/// function inserts its own `chosen_ps` into it on success). Both are
/// checked; refuses to choose a colliding id either way.
/// Returns the number of reads newly phased this way.
int emit_independent_gap_block(std::vector<PhasingChunk>& chunks,
                               const PhasingChunk& proposal,
                               const PhaseGap& gap, const Options& opts,
                               const GapReadIndex& read_index,
                               int min_reads,
                               hts_pos_t* emitted_ps = nullptr,
                               std::set<hts_pos_t>* emitted_this_round = nullptr);

/// Attach verified heterozygous sites with no read-supported phase set to
/// stitched read blocks, without changing any read assignment.
int anchor_orphan_msa_sites(std::vector<PhasingChunk>& chunks, const Options& opts,
                            const std::vector<PhaseGap>* gaps = nullptr);

/// Add gap-driven MSA windows, preserving detected noisy intervals intact.
void prepare_gap_msa_regions(PhasingChunk& chunk, hts_pos_t beg, hts_pos_t end);

/// Admit SNPs only at tier 2; retain them and admit indels at tier 3.
void run_gap_msa_tier(PhasingChunk& chunk, const Options& opts,
                      hts_pos_t beg, hts_pos_t end, bool snp_only);

} // namespace pgphase_collect
#endif
