#pragma once

/**
 * @file plegma_msa.hpp
 * @brief plegma as the noisy-region MSA (optional; abPOA stays the default).
 *
 * The default noisy-region MSA. Built with plegma unless `make PLEGMA=0` (third_party/plegma's build/libplegma.so,
 * from `pixi run build-capi`; PGPHASE_PLEGMA defined), chosen per run by --msa (abpoa: longcallD's abPOA). With
 * plegma, a homopolymer model decides run lengths: built in for HiFi and ONT R10.4.1 reads, or --msa-model FILE.
 * Where abPOA gives one MSA per cluster, plegma gives each read's own alignment to its cluster's consensus, as two
 * rows of one length (base codes 0-4, 5 a gap), the pair `make_cons_read_aln_str` takes.
 */

#include "align.hpp"

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace pgphase_collect {

/// Whether this pgphase was built with plegma.
bool plegma_built();

/// The run's noisy-region MSA from opts.msa ("": plegma when built with it, else abPOA) and, with plegma, its
/// homopolymer model from opts.msa_model ("": the built-in one for the read technology, none for short reads; "none";
/// or a shasta2 csv). PGPHASE_MSA and PGPHASE_PLEGMA_MODEL override both. Called once, before any worker thread;
/// without it, the first use takes the defaults. Throws std::runtime_error on a bad choice or model.
void configure_msa(const Options& opts);

/// Whether the run's noisy-region MSA is plegma's.
bool use_plegma_msa();

/// The run's MSA as a line for logs: "plegma 0.3.0, homopolymer model ..." or "abPOA".
std::string msa_description();

/// The MSAs this build has, for --version: "plegma 0.3.0 (default), abPOA" or "abPOA (built without plegma)".
std::string msa_backends();

/// Whether the MSA's het calls stand without the reads' phasing (CandidateVariant::msa_genotyped): plegma with a
/// homopolymer model (PGPHASE_PLEGMA_MODEL).
bool msa_genotypes_alone();

/// With PGPHASE_MSA_DUMP_READS set: the inputs of this thread's last plegma call ("clustered" or "phased", then per
/// read letters:anchor:strand:haplotype, tab-separated), cleared; empty otherwise.
std::string take_plegma_inputs();

/// A read's alignment to its cluster's consensus as an MSA would show the two: equal-length rows, codes 0-4, 5 a gap.
struct PlegmaRows {
    std::vector<uint8_t> cons;
    std::vector<uint8_t> read;
};

/**
 * @brief The consensus sequences of a read set by plegma: at most max_n_cons clusters (abPOA's -d), each with its
 * reads and their rows.
 *
 * read_seqs are base codes 0-4; read_covers the reads' noisy-region cover flags (both ends, the left or the right
 * end alone, or neither), as plegma's anchors; read_strands 1 for a read reverse to the reference (for the
 * homopolymer model PGPHASE_PLEGMA_MODEL names, shasta2's csv, when set). cons_seqs[c] is cluster c's consensus (codes); clu_read_ids[c] its
 * reads as read_ids names them, in input order, and rows[c][k] read k's rows. Returns the number of clusters (0 on
 * failure, or 1 or 2), at opts' scoring and min_af.
 */
int plegma_msa_cons(const Options& opts, const std::vector<int>& read_ids,
                    const std::vector<uint8_t*>& read_seqs, const std::vector<int>& read_lens,
                    const std::vector<int>& read_covers, const std::vector<uint8_t>& read_strands,
                    int max_n_cons,
                    std::vector<std::vector<uint8_t>>& cons_seqs,
                    std::array<int, 2>& clu_n_seqs,
                    std::array<std::vector<int>, 2>& clu_read_ids,
                    std::array<std::vector<PlegmaRows>, 2>& rows);

/**
 * @brief The two haplotypes' consensus sequences of reads pgphase has phased, by plegma in one call: read r on
 * haplotype read_haps[r] (0 or 1). Each haplotype's consensus, reads (as read_ids names them, in input order) and
 * rows as plegma_msa_cons gives them; where the two differ only by homopolymer lengths, with a model, plegma decides
 * each such run on both haplotypes' reads. Returns 2, or 0 on failure.
 */
int plegma_phased_msa_cons(const Options& opts, const std::vector<int>& read_ids,
                           const std::vector<uint8_t*>& read_seqs, const std::vector<int>& read_lens,
                           const std::vector<int>& read_covers, const std::vector<uint8_t>& read_strands,
                           const std::vector<int>& read_haps,
                           std::array<std::vector<uint8_t>, 2>& cons_seqs,
                           std::array<int, 2>& clu_n_seqs,
                           std::array<std::vector<int>, 2>& clu_read_ids,
                           std::array<std::vector<PlegmaRows>, 2>& rows);

}  // namespace pgphase_collect
