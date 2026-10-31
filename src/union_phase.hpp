#pragma once

// Union gap phasing: what collect-graph-variation does with --bam.
//
// The catalog's sites are phased together with the alignment's sample-specific
// heterozygotes: inject_alignment_sites adds them to a graph chunk, the clean
// k-means runs as usual, and phase_chunk_with_alignment_sites then solves the
// whole read x site matrix by EM, adding local haplotype windows over the
// alignment's noisy loci. See docs/IMPLEMENTATION.md, "Union gap phasing".

#include "graph_bam_adapter.hpp"

#include <cstddef>
#include <tuple>
#include <cstdint>
#include <unordered_map>
#include <vector>

namespace pgphase_collect {

/// Error of one read's call at one candidate, from how decisively its own
/// sequence realigns to the called allele: (read << 32 | candidate) -> error.
using ObservationErrors = std::unordered_map<uint64_t, double>;
inline uint64_t observation_key(size_t read, size_t candidate) {
    return (static_cast<uint64_t>(read) << 32) | static_cast<uint64_t>(candidate);
}

/// Merge an alignment chunk's sample-specific heterozygotes into a graph chunk.
///
/// Takes the alignment's clean heterozygotes and MSA-verified noisy
/// heterozygotes that no voting graph row overlaps and whose two most observed
/// alleles each have confidently mapped support, one row per left-aligned
/// allele, with the alignment's category bits and per-read alleles. Reads only
/// the alignment holds are added. Every index-parallel array is rebuilt.
/// Returns the number of injected sites.
size_t inject_alignment_sites(GraphChunkBuildResult& graph_chunk, const PhasingChunk& bam,
                              const char* contig_name);

/// A locus-level site for the EM: each spanning read is assigned to one of
/// two local haplotype sequences over a whole window, with its own error.
struct LocusWindowSite {
    hts_pos_t pos = 0;  // 1-based window start
    hts_pos_t end = 0;  // 1-based inclusive window end
    std::vector<std::tuple<size_t, int, float>> observations;  // graph read, side, error
    // An island-only window is never an EM site: it can only phase reads the
    // solve leaves unlabelled, in a block of its own (see island blocks).
    bool island_only = false;
    // An allele window over one catalog repeat-indel row (build_allele_window_sites): the row, its side 0 the row's
    // allele 0 and side 1 its allele 1. -1 for every other window.
    long catalog_row = -1;
    // The reads' median edit distance to each side's sequence (catalog-row windows): a read nearer the ALT than the
    // REF need not carry that ALT (another length of the repeat).
    int median_edits[2] = {0, 0};
    // Written by phase_chunk_by_global_em: which side haplotype 1 carries
    // (-1 when not phased), the block, and the learned error.
    int phase = -1;
    hts_pos_t block = 0;
    double error = 1.0;
};

/// Local haplotype windows over the alignment's noisy-region loci.
///
/// Loci are the alignment's noisy calls merged, padded and extended through
/// tandem repeats. The two haplotype sequences are medoids of spanning reads
/// the graph chunk's phased sites already label (dominant phase set); a read
/// is assigned when it is at least two edits closer to one of them. A window
/// is kept only when the two sequences differ and both sides hold reads.
std::vector<LocusWindowSite> build_locus_window_sites(const PhasingChunk& bam,
                                                      const PhasingChunk& graph);

/// Allele windows for heterozygous indels the EM does not otherwise observe:
/// the alignment's noisy-region and repeat-context calls, and the catalog's
/// non-voting repeat indels (alleles from their VCF-form metadata). Each spanning read is realigned to the window's two exact allele
/// sequences, reference context with REF and with ALT, and joins the closer
/// one. Windows overlapping `taken` windows or voting indel rows are skipped.
/// Alignment candidates read by more than kUnionMaxDepthRatio times their chunk's median depth: a collapsed
/// duplication, whose two "alleles" separate paralogous copies, not haplotypes (union_select.cpp). Injection keeps
/// them out of the EM; the allele windows must too.
std::vector<char> collapsed_duplication_sites(const PhasingChunk& bam);

std::vector<LocusWindowSite> build_allele_window_sites(const PhasingChunk& bam,
                                                       const GraphChunkBuildResult& graph_chunk,
                                                       const std::vector<LocusWindowSite>& taken);

/// Re-call each read's allele at the injected indel sites by realignment:
/// the read's segment over the site's window (indel, tandem-repeat extension,
/// flank) is aligned to the reference context with REF and with ALT and takes
/// the closer allele, or no call on a tie. Only reads covering the whole window
/// with a call slot at the site are re-called. Returns calls changed.
size_t realign_indel_observations(PhasingChunk& chunk, const PhasingChunk& bam,
                                  ObservationErrors* errors = nullptr);

/// The alignment's base call outranks the graph's at catalog SNP rows: a read
/// whose graph allele disagrees with its alignment call takes the alignment's
/// allele (or no call when the row cannot carry that base), and a graph call
/// at a site the alignment calls but not for that read gets a per-read error
/// (kGraphOnlyCallError). Returns calls changed.
size_t confirm_catalog_snp_calls(GraphChunkBuildResult& graph_chunk, const PhasingChunk& bam,
                                 ObservationErrors* errors);

/// Add the alignment's primary reads (MAPQ >= 20) that the graph chunk lacks,
/// with no observations. Returns a per-read flag marking the added reads.
std::vector<char> add_alignment_only_reads(PhasingChunk& chunk, const PhasingChunk& bam);

/// Give the flagged reads allele calls at voting sites by realignment to the
/// two candidate alleles in reference context (catalog rows from their
/// VCF-form metadata, injected rows from their key): one edit decides a SNP,
/// an indel needs two. Existing calls are never changed. Returns calls added.
size_t fill_missing_observations(GraphChunkBuildResult& graph_chunk, const PhasingChunk& bam,
                                 const std::vector<char>& only_reads);

/// Sites that may be phased but must not carry a join between blocks alone:
/// their two alleles differ by one base of length (for injected 1 bp indels,
/// one homopolymer base), the dominant HiFi error. Indexed by candidate.
std::vector<char> bridge_weak_sites(const GraphChunkBuildResult& graph_chunk, const PhasingChunk& bam);

/// Length windows over repeat loci: at every tandem-repeat tract holding an
/// alignment indel candidate (any category) and no voting site, each spanning
/// read's length change across the tract (between unique flanks) is
/// genotyped against a stutter model learned from the chunk's own loci; a
/// confident two-length locus becomes a window contrasting the two lengths.
/// They enter the EM island-only (see island blocks).
std::vector<LocusWindowSite> build_length_window_sites(const PhasingChunk& bam,
                                                       const GraphChunkBuildResult& graph_chunk,
                                                       const std::vector<LocusWindowSite>& taken);

/// Apply the union EM's local island blocks (PhasingChunk::island_*) after the
/// rescue: reads it left unphased take their island label through the gap_*
/// channel, and island sites still unphased take their island phase. Returns
/// reads labelled.
size_t apply_island_labels(PhasingChunk& chunk);

/// Phase a chunk by EM over its whole read x site matrix.
///
/// Every voting site carrying two observed alleles gets a phase and its own
/// error rate; every read gets a haplotype posterior, initialised from the
/// chunk's read haplotypes. Switch moves flip all sites right of a boundary
/// whenever that raises the likelihood. Phase blocks are cut where flipping
/// across a boundary would cost little; sites whose learned error stays high
/// are left unphased; a read is labelled in the block where its evidence is
/// strongest. Returns the number of phase blocks.
size_t phase_chunk_by_global_em(PhasingChunk& chunk,
                                std::vector<LocusWindowSite>* loci = nullptr,
                                const std::vector<char>* bridge_weak = nullptr,
                                const ObservationErrors* errors = nullptr);

/// After the clean stage: EM, local haplotype windows seeded by its read
/// labels, then EM over sites and windows together.
void phase_chunk_with_alignment_sites(GraphChunkBuildResult& graph_chunk, const PhasingChunk* bam);

} // namespace pgphase_collect
