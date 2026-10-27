/**
 * @file plegma_msa.cpp
 * @brief plegma as the noisy-region MSA (see plegma_msa.hpp).
 */

#include "plegma_msa.hpp"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>

#ifdef PGPHASE_PLEGMA
#include "plegma.hpp"
#endif

namespace pgphase_collect {

bool use_plegma_msa() {
#ifdef PGPHASE_PLEGMA
    static const bool on = [] {
        const char* msa = std::getenv("PGPHASE_MSA");
        return msa != nullptr && std::strcmp(msa, "plegma") == 0;
    }();
    return on;
#else
    return false;
#endif
}

#ifdef PGPHASE_PLEGMA
namespace {

constexpr uint8_t kGap = 5;

uint8_t code_of(char letter) {
    switch (letter) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return 4;
    }
}

plegma::Anchor anchor_of(int cover) {
    if (noisyIsBothCover(cover)) return plegma::Anchor::Both;
    if (noisyIsLeftCover(cover)) return plegma::Anchor::Left;
    if (noisyIsRightCover(cover)) return plegma::Anchor::Right;
    return plegma::Anchor::Inside;
}

/// The two rows of a read's alignment to the consensus: consensus letters outside its span against gaps, read
/// letters outside its span against gaps, the CIGAR between.
PlegmaRows rows_of(const std::vector<uint8_t>& cons, const uint8_t* read, int read_len,
                   const plegma::Placement& placed) {
    PlegmaRows out;
    const std::size_t width = cons.size() + static_cast<std::size_t>(read_len);
    out.cons.reserve(width);
    out.read.reserve(width);
    for (std::size_t j = 0; j < placed.readStart; ++j) {
        out.cons.push_back(kGap);
        out.read.push_back(read[j]);
    }
    for (std::size_t c = 0; c < placed.consensusStart; ++c) {
        out.cons.push_back(cons[c]);
        out.read.push_back(kGap);
    }
    std::size_t i = placed.consensusStart;
    std::size_t j = placed.readStart;
    std::size_t count = 0;
    for (char op : placed.cigar) {
        if (op >= '0' && op <= '9') {
            count = count * 10 + static_cast<std::size_t>(op - '0');
            continue;
        }
        for (std::size_t k = 0; k < count; ++k) {
            if (op == 'I') {
                out.cons.push_back(kGap);
                out.read.push_back(read[j++]);
            } else if (op == 'D') {
                out.cons.push_back(cons[i++]);
                out.read.push_back(kGap);
            } else {
                out.cons.push_back(cons[i++]);
                out.read.push_back(read[j++]);
            }
        }
        count = 0;
    }
    for (; i < cons.size(); ++i) {
        out.cons.push_back(cons[i]);
        out.read.push_back(kGap);
    }
    for (; j < static_cast<std::size_t>(read_len); ++j) {
        out.cons.push_back(kGap);
        out.read.push_back(read[j]);
    }
    return out;
}

/// The homopolymer model PGPHASE_PLEGMA_MODEL names (shasta2's csv), loaded once; none when unset.
const plegma::Model* model() {
    static const std::unique_ptr<plegma::Model> loaded = []() -> std::unique_ptr<plegma::Model> {
        const char* path = std::getenv("PGPHASE_PLEGMA_MODEL");
        if (path == nullptr || *path == 0) return nullptr;
        try {
            return std::make_unique<plegma::Model>(plegma::Model::load(path));
        } catch (const std::runtime_error& e) {
            // A model asked for and not loaded is a broken run, not one without the model.
            std::fprintf(stderr, "pgphase: PGPHASE_PLEGMA_MODEL: %s\n", e.what());
            std::exit(1);
        }
    }();
    return loaded.get();
}

}  // namespace
#endif

#ifdef PGPHASE_PLEGMA
namespace {

struct PlegmaInput {
    std::vector<std::string> letters;
    std::vector<plegma::Read> reads;
};

/// The reads as plegma takes them: letters, anchors from the cover flags, strands.
void input_of(const std::vector<uint8_t*>& read_seqs, const std::vector<int>& read_lens,
              const std::vector<int>& read_covers, const std::vector<uint8_t>& read_strands, PlegmaInput& in) {
    const std::size_t n = read_seqs.size();
    in.letters.assign(n, std::string());
    in.reads.clear();
    in.reads.reserve(n);
    for (std::size_t r = 0; r < n; ++r) {
        in.letters[r].resize(static_cast<std::size_t>(read_lens[r]));
        for (int i = 0; i < read_lens[r]; ++i)
            in.letters[r][static_cast<std::size_t>(i)] = "ACGTN"[std::min<int>(read_seqs[r][i], 4)];
        const plegma::Strand strand = r < read_strands.size()
                                          ? (read_strands[r] ? plegma::Strand::Reverse : plegma::Strand::Forward)
                                          : plegma::Strand::Unknown;
        in.reads.push_back({in.letters[r], anchor_of(read_covers[r]), strand});
    }
}

/// pgphase's scoring and the PGPHASE_PLEGMA_* knobs as plegma's options.
plegma::Options options_of(const Options& opts) {
    plegma::Options options = plegma::Options()
                                  .pinned(false)
                                  .scoring(opts.match, opts.mismatch)
                                  .gaps(opts.gap_open1, opts.gap_ext1, opts.gap_open2, opts.gap_ext2)
                                  // On HG002 HiFi chr12, 18 and 20, 5 dB against plegma's 15: 22 fewer true hets of
                                  // 28,259 for 77 fewer false of 282, the phasing as good.
                                  .twoLengthsDb(5.0f);
    // PGPHASE_PLEGMA_NOISE_RUN=N: clusters apart only by single letters in runs of N or more are one (0: never).
    if (const char* run = std::getenv("PGPHASE_PLEGMA_NOISE_RUN")) options.noiseRun(std::atoi(run));
    // PGPHASE_PLEGMA_TWO_DB=X: dB two homopolymer lengths must gain over one to stay apart.
    if (const char* db = std::getenv("PGPHASE_PLEGMA_TWO_DB")) options.twoLengthsDb(static_cast<float>(std::atof(db)));
    // PGPHASE_PLEGMA_ROUNDS=E,L: polishing rounds by edit distance, then by likelihood (plegma poa -P).
    if (const char* rounds = std::getenv("PGPHASE_PLEGMA_ROUNDS")) {
        int edit = 50, likely = 20;
        if (std::sscanf(rounds, "%d,%d", &edit, &likely) >= 1) options.polishRounds(edit).likelyRounds(likely);
    }
    return options;
}

/// plegma's clusters as pgphase's: consensus in base codes, reads by read_ids, each read's rows. Returns the
/// clusters kept (at most 2).
int outputs_of(const plegma::Clustering& got, const std::vector<int>& read_ids,
               const std::vector<uint8_t*>& read_seqs, const std::vector<int>& read_lens,
               std::vector<std::vector<uint8_t>>& cons_seqs, std::array<int, 2>& clu_n_seqs,
               std::array<std::vector<int>, 2>& clu_read_ids, std::array<std::vector<PlegmaRows>, 2>& rows) {
    const int n_cons = static_cast<int>(std::min<std::size_t>(got.consensus.size(), 2));
    cons_seqs.assign(static_cast<std::size_t>(n_cons), {});
    for (int c = 0; c < 2; ++c) {
        clu_n_seqs[static_cast<std::size_t>(c)] = 0;
        clu_read_ids[static_cast<std::size_t>(c)].clear();
        rows[static_cast<std::size_t>(c)].clear();
    }
    for (int c = 0; c < n_cons; ++c)
        for (char letter : got.consensus[static_cast<std::size_t>(c)])
            cons_seqs[static_cast<std::size_t>(c)].push_back(code_of(letter));
    for (std::size_t r = 0; r < read_seqs.size(); ++r) {
        const plegma::Placement& placed = got.placements[r];
        if (placed.cluster < 0 || placed.cluster >= n_cons) continue;
        const std::size_t c = static_cast<std::size_t>(placed.cluster);
        clu_read_ids[c].push_back(read_ids[r]);
        rows[c].push_back(rows_of(cons_seqs[c], read_seqs[r], read_lens[r], placed));
        ++clu_n_seqs[c];
    }
    return n_cons;
}

}  // namespace
#endif

int plegma_msa_cons(const Options& opts, const std::vector<int>& read_ids,
                    const std::vector<uint8_t*>& read_seqs, const std::vector<int>& read_lens,
                    const std::vector<int>& read_covers, const std::vector<uint8_t>& read_strands,
                    int max_n_cons,
                    std::vector<std::vector<uint8_t>>& cons_seqs,
                    std::array<int, 2>& clu_n_seqs,
                    std::array<std::vector<int>, 2>& clu_read_ids,
                    std::array<std::vector<PlegmaRows>, 2>& rows) {
#ifdef PGPHASE_PLEGMA
    PlegmaInput in;
    input_of(read_seqs, read_lens, read_covers, read_strands, in);
    plegma::Clustering got;
    try {
        got = plegma::clusteredConsensus(in.reads, model(), options_of(opts), std::min(max_n_cons, 2),
                                         static_cast<float>(opts.min_af));
    } catch (const std::runtime_error& e) {
        if (opts.verbose >= 1) std::fprintf(stderr, "plegma MSA: %s\n", e.what());
        return 0;
    }
    return outputs_of(got, read_ids, read_seqs, read_lens, cons_seqs, clu_n_seqs, clu_read_ids, rows);
#else
    (void)opts; (void)read_ids; (void)read_seqs; (void)read_lens; (void)read_covers; (void)read_strands; (void)max_n_cons;
    (void)cons_seqs; (void)clu_n_seqs; (void)clu_read_ids; (void)rows;
    return 0;
#endif
}

int plegma_phased_msa_cons(const Options& opts, const std::vector<int>& read_ids,
                           const std::vector<uint8_t*>& read_seqs, const std::vector<int>& read_lens,
                           const std::vector<int>& read_covers, const std::vector<uint8_t>& read_strands,
                           const std::vector<int>& read_haps,
                           std::array<std::vector<uint8_t>, 2>& cons_seqs,
                           std::array<int, 2>& clu_n_seqs,
                           std::array<std::vector<int>, 2>& clu_read_ids,
                           std::array<std::vector<PlegmaRows>, 2>& rows) {
#ifdef PGPHASE_PLEGMA
    PlegmaInput in;
    input_of(read_seqs, read_lens, read_covers, read_strands, in);
    plegma::Clustering got;
    try {
        got = plegma::phasedConsensus(in.reads, read_haps, model(), options_of(opts));
    } catch (const std::runtime_error& e) {
        if (opts.verbose >= 1) std::fprintf(stderr, "plegma MSA: %s\n", e.what());
        return 0;
    }
    std::vector<std::vector<uint8_t>> cons;
    const int n_cons = outputs_of(got, read_ids, read_seqs, read_lens, cons, clu_n_seqs, clu_read_ids, rows);
    if (n_cons != 2) return 0;
    cons_seqs[0] = std::move(cons[0]);
    cons_seqs[1] = std::move(cons[1]);
    return 2;
#else
    (void)opts; (void)read_ids; (void)read_seqs; (void)read_lens; (void)read_covers; (void)read_strands;
    (void)read_haps; (void)cons_seqs; (void)clu_n_seqs; (void)clu_read_ids; (void)rows;
    return 0;
#endif
}

}  // namespace pgphase_collect
