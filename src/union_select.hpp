#ifndef PGPHASE_UNION_SELECT_HPP
#define PGPHASE_UNION_SELECT_HPP

// The injection plan of union gap phasing (union_select.cpp): which alignment candidates enter the graph chunk and
// how the catalog's rows are written beside them. inject_alignment_sites (union_inject.cpp) applies it.

#include "union_internal.hpp"

#include <set>
#include <unordered_map>
#include <vector>

namespace pgphase_collect {

// Which alignment candidates enter the graph chunk, and as what.
struct Selection {
    std::vector<size_t> selected;
    std::set<size_t> island_only;    // selected only as island sites
    std::set<size_t> island_written;  // of them, the depth gate's: written unphased with the MSA's genotype
    std::set<size_t> low_mapq_only;  // clean hets admitted by the exemption
    std::unordered_map<size_t, std::vector<size_t>> complement_rows;  // bam candidate -> graph rows it overlaps
    std::set<size_t> replaced_rows;  // graph rows an MSA call of an allele they do not spell replaces in the VCF
    std::set<size_t> unwritten_rows;  // graph rows whose allele the consensus sequences of their MSA region lack
    std::set<size_t> hom_rows;        // graph rows both consensus sequences of their MSA region carry: written 1/1
    std::set<size_t> hom_records;     // homozygous calls (the MSA's, the alignment's): written 1/1, no vote
    std::unordered_map<size_t, const char*> row_reason;  // why a graph row is unwritten or 1/1 (stage trace)
    std::set<size_t> spelled_hom;     // clean het calls no read spells the reference at (spell_homozygous_calls)
    std::set<size_t> additive;        // calls adding a change to a catalog row's haplotype (find_additive_calls)
    std::vector<const char*> decision;  // per candidate (diagnostic only)
};

/// The plan for one chunk: every check of the catalog's rows against the MSA's consensus sequences and the reads'
/// spellings, then the alignment candidates selected (select_sites). Reads the graph chunk and the alignment solve and
/// changes neither, so a test can plan on a fixture and inspect the plan. `diagnostics` keeps each candidate's
/// decision (Selection::decision) for write_injection_diagnostics.
Selection plan_injection(const GraphChunkBuildResult& graph_chunk, const PhasingChunk& bam, bool diagnostics);

/// Diagnostic, PGPHASE_DIAG=DIR: each alignment candidate's injection decision (inj.TID.BEG.tsv), each read's calls in
/// the alignment solve (alnreads.TID.BEG.tsv) and the graph chunk's rows before the merge (graph.TID.BEG.tsv).
void write_injection_diagnostics(const char* diag_dir, const GraphChunkBuildResult& graph_chunk,
                                 const PhasingChunk& bam, const Selection& sel);

/// An injected MSA call's phase in the alignment solve, from the clusters of its region (union_select.cpp).
void aln_phase_from_msa_clusters(const PhasingChunk& bam, CandidateVariant& cand);

}  // namespace pgphase_collect

#endif  // PGPHASE_UNION_SELECT_HPP
