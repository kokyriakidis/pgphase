#pragma once

/**
 * @file plegma_msa.hpp
 * @brief plegma as the noisy-region MSA (optional; abPOA stays the default).
 *
 * Built only with `make PLEGMA_DIR=<plegma checkout>` (its build/libplegma.so from `pixi run build-capi`), which
 * defines PGPHASE_PLEGMA; chosen at run time with PGPHASE_MSA=plegma. Where abPOA gives one MSA per cluster, plegma
 * gives each read's own alignment to its cluster's consensus, as two rows of one length (base codes 0-4, 5 a gap),
 * the pair `make_cons_read_aln_str` takes.
 */

#include "align.hpp"

#include <array>
#include <cstdint>
#include <vector>

namespace pgphase_collect {

/// Whether the noisy-region MSA is plegma's: built with it and PGPHASE_MSA=plegma.
bool use_plegma_msa();

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
