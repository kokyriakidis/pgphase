#ifndef PGPHASE_REFINED_ALIGNMENT_HPP
#define PGPHASE_REFINED_ALIGNMENT_HPP

// --refine-aln in collect-graph-variation: each read's alignment inside the noisy regions the MSA resolved, rewritten
// from its alignment to that region's consensus (longcallD's --refine-aln), taken from the union's alignment solve.
// The writer is aligned_read_writer.hpp.

#include "phasing_types.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include <htslib/sam.h>

namespace pgphase_collect {

/// A read's alignment as the MSA rewrote it: kept only where it differs from the input record's (with = and X
/// read as M). Positions are 0-based, as bam1_core_t::pos; the CIGAR in the input record's notation (M, or = and X).
struct RefinedAlignment {
    std::string qname;
    hts_pos_t orig_pos = 0;
    hts_pos_t pos = 0;
    std::vector<uint32_t> cigar;
};

struct RefineCounts {
    std::size_t refined = 0;    // rewritten, kept
    std::size_t unchanged = 0;  // the MSA's alignment is the input's
    std::size_t rejected = 0;   // failed the check against the read's bases, or moved too far: input kept
};

/// The farthest refinement may move a read's start: further is rejected. It bounds how far back the aligned writer
/// holds records to keep its output sorted (an MSA region is at most kUnionMsaMaxRegionLength, 10 kb).
constexpr hts_pos_t kRefineMaxShift = 20000;

/// The alignments the noisy-region MSA rewrote in an alignment chunk (its reads' digars after
/// update_digars_from_aln_str, with Options::refine_msa_digars), each checked against the read's own bases and the
/// chunk's reference; then frees what the rewrite kept (digars, base qualities, noisy intervals).
std::vector<RefinedAlignment> take_refined_alignments(PhasingChunk& chunk, RefineCounts& counts);

}  // namespace pgphase_collect

#endif
