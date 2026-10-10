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
#include <vector>

namespace pgphase_collect {

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
    // Written by phase_chunk_by_global_em: which side haplotype 1 carries
    // (-1 when not phased), the block, and the learned error.
    mutable int phase = -1;
    mutable hts_pos_t block = 0;
    mutable double error = 1.0;
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
std::vector<LocusWindowSite> build_allele_window_sites(const PhasingChunk& bam,
                                                       const GraphChunkBuildResult& graph_chunk,
                                                       const std::vector<LocusWindowSite>& taken);

/// Re-call each read's allele at the injected indel sites by realignment:
/// the read's segment over the site's window (indel, tandem-repeat extension,
/// flank) is aligned to the reference context with REF and with ALT and takes
/// the closer allele, or no call on a tie. Only reads covering the whole window
/// with a call slot at the site are re-called. Returns calls changed.
size_t realign_indel_observations(PhasingChunk& chunk, const PhasingChunk& bam);

/// Add the alignment's primary reads (MAPQ >= 20) that the graph chunk lacks,
/// with no observations. Returns a per-read flag marking the added reads.
std::vector<char> add_alignment_only_reads(PhasingChunk& chunk, const PhasingChunk& bam);

/// Give the flagged reads allele calls at voting sites by realignment to the
/// two candidate alleles in reference context (catalog rows from their
/// VCF-form metadata, injected rows from their key): one edit decides a SNP,
/// an indel needs two. Existing calls are never changed. Returns calls added.
size_t fill_missing_observations(GraphChunkBuildResult& graph_chunk, const PhasingChunk& bam,
                                 const std::vector<char>& only_reads);

/// Last resort for reads the EM left unlabelled, from the labelled reads
/// themselves: inside each phase set, an alignment event (base, insertion or
/// deletion, left-normalized) that one haplotype's covering reads carry and the
/// other's lack marks a phased heterozygote, whether or not any caller found
/// it. An unlabelled read votes at every marker it covers and is labelled when
/// its votes agree. Decisions go to `deferred_read_labels`. Returns reads decided.
size_t label_reads_from_haplotype_consensus(GraphChunkBuildResult& graph_chunk, const PhasingChunk& bam);

/// Last resort for reads the EM left unlabelled: at each MSA-verified indel of
/// the alignment solve, every spanning read is called by realignment to the two
/// exact alleles. Where at least four already labelled reads of one phase set
/// agree with the calls at 90% or more, the site's orientation is taken from
/// them, and an unlabelled read takes the haplotype its allele implies. A read
/// whose verified indels disagree stays unlabelled. Decisions are recorded in
/// `deferred_read_labels` relative to an anchor read of the same phase set and
/// applied by apply_deferred_read_labels after chunk stitching, so they feed
/// back into nothing. Returns reads decided.
size_t label_reads_from_verified_indels(GraphChunkBuildResult& graph_chunk, const PhasingChunk& bam);

/// Apply the deferred last-resort labels: each read takes its anchor's final
/// phase set and the same or the other haplotype. Returns labels applied.
size_t apply_deferred_read_labels(GraphChunkBuildResult& graph_chunk);

/// Sites that may be phased but must not carry a join between blocks alone:
/// their two alleles differ by one base of length (for injected 1 bp indels,
/// one homopolymer base), the dominant HiFi error. Indexed by candidate.
std::vector<char> bridge_weak_sites(const GraphChunkBuildResult& graph_chunk, const PhasingChunk& bam);

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
                                const std::vector<LocusWindowSite>* loci = nullptr,
                                const std::vector<char>* bridge_weak = nullptr);

/// After the clean stage: EM, local haplotype windows seeded by its read
/// labels, then EM over sites and windows together.
void phase_chunk_with_alignment_sites(GraphChunkBuildResult& graph_chunk, const PhasingChunk* bam);

} // namespace pgphase_collect
