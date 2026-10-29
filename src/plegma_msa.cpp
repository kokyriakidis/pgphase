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
#include <mutex>
#include <stdexcept>
#include <string>

#ifdef PGPHASE_PLEGMA
#include "plegma.hpp"
#endif

namespace pgphase_collect {

#ifdef PGPHASE_PLEGMA
// The homopolymer models built in (src/msa_models.cpp, generated from models/ by make).
extern const char kHifiHomopolymerModel[];
extern const std::size_t kHifiHomopolymerModelSize;
extern const char kOntHomopolymerModel[];
extern const std::size_t kOntHomopolymerModelSize;
#endif

namespace {

/// The noisy-region MSA of this run, chosen once by configure_msa (or on first use, with the defaults).
struct MsaChoice {
    bool plegma = false;
    std::string model_name;
#ifdef PGPHASE_PLEGMA
    std::unique_ptr<plegma::Model> model;
#endif
};

MsaChoice& msa_choice() {
    static MsaChoice choice;
    return choice;
}

std::once_flag msa_defaults;
bool msa_configured = false;

void choose_msa(const Options& opts) {
    MsaChoice& choice = msa_choice();
    std::string msa = opts.msa;
    // PGPHASE_MSA and PGPHASE_PLEGMA_MODEL override the options (the development loop sets them).
    if (const char* env = std::getenv("PGPHASE_MSA"); env != nullptr && *env != 0) msa = env;
    if (msa.empty()) msa = plegma_built() ? "plegma" : "abpoa";
    if (msa != "plegma" && msa != "abpoa") throw std::runtime_error("--msa " + msa + ": not plegma or abpoa");
    if (msa == "plegma" && !plegma_built())
        throw std::runtime_error("--msa plegma: this pgphase was built without plegma (make PLEGMA=0); rebuild with it, or use --msa abpoa");
    choice.plegma = msa == "plegma";
    choice.model_name.clear();
#ifdef PGPHASE_PLEGMA
    choice.model.reset();
    if (!choice.plegma) return;
    std::string path = opts.msa_model;
    if (const char* env = std::getenv("PGPHASE_PLEGMA_MODEL"); env != nullptr && *env != 0) path = env;
    if (path == "none") return;
    if (!path.empty()) {
        choice.model = std::make_unique<plegma::Model>(plegma::Model::load(path));
        choice.model_name = path;
    } else if (opts.is_ont()) {
        choice.model = std::make_unique<plegma::Model>(
            plegma::Model::parse(std::string_view(kOntHomopolymerModel, kOntHomopolymerModelSize)));
        choice.model_name = "built-in ONT R10.4.1 sup 5.0.0 (shasta2)";
    } else if (!opts.is_short_reads()) {
        choice.model = std::make_unique<plegma::Model>(
            plegma::Model::parse(std::string_view(kHifiHomopolymerModel, kHifiHomopolymerModelSize)));
        choice.model_name = "built-in HiFi (HG002)";
    }
#endif
}

const MsaChoice& current_msa() {
    std::call_once(msa_defaults, [] {
        if (!msa_configured) choose_msa(Options{});
    });
    return msa_choice();
}

}  // namespace

bool plegma_built() {
#ifdef PGPHASE_PLEGMA
    return true;
#else
    return false;
#endif
}

void configure_msa(const Options& opts) {
    msa_configured = true;
    choose_msa(opts);
}

bool use_plegma_msa() { return current_msa().plegma; }

std::string msa_backends() {
#ifdef PGPHASE_PLEGMA
    return "plegma " + plegma::version() + " (default), abPOA";
#else
    return "abPOA (built without plegma)";
#endif
}

std::string msa_description() {
    const MsaChoice& choice = current_msa();
    if (!choice.plegma) return "abPOA";
    std::string out = "plegma";
#ifdef PGPHASE_PLEGMA
    out += " " + plegma::version();
#endif
    out += choice.model_name.empty() ? ", no homopolymer model" : ", homopolymer model " + choice.model_name;
    return out;
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

/// The run's homopolymer model (configure_msa), or none.
const plegma::Model* model() { return current_msa().model.get(); }

}  // namespace
#endif

bool msa_genotypes_alone() {
#ifdef PGPHASE_PLEGMA
    return use_plegma_msa() && model() != nullptr;
#else
    return false;
#endif
}

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

/// PGPHASE_MSA_DUMP_READS: this thread's last call's inputs, for dump_msa_call to write beside its region.
thread_local std::string last_inputs;

/// The call's kind, then per read its letters, anchor, strand and haplotype (-1: none), colon-separated.
void remember_inputs(const char* kind, const PlegmaInput& in, const std::vector<int>* haps) {
    static const bool on = std::getenv("PGPHASE_MSA_DUMP_READS") != nullptr;
    if (!on) return;
    last_inputs = kind;
    for (std::size_t r = 0; r < in.reads.size(); ++r) {
        last_inputs += '\t';
        last_inputs += in.letters[r];
        last_inputs += ':' + std::to_string(static_cast<int>(in.reads[r].anchor)) + ':' +
                       std::to_string(static_cast<int>(in.reads[r].strand)) + ':' +
                       std::to_string(haps != nullptr && r < haps->size() ? (*haps)[r] : -1);
    }
}

/// The no-PS path's `two_lengths_db`: clusters of every spanning read apart only by homopolymer lengths stay apart
/// only above it. HG002 HiFi chr12/18/20 against the PS-aware path's 5 dB: at 8, false het regions 221 -> 205, true
/// 29,450 -> 29,446, arm reads +38 correct, -24 discordant; 12 dB 197 false for 29,441 true, 20 dB 186 for 29,430.
constexpr float kClusteredTwoLengthsDb = 8.0f;

/// pgphase's scoring and the PGPHASE_PLEGMA_* knobs as plegma's options; `two_lengths_db` the PS-aware path's 5 dB
/// unless given (the no-PS path's is kClusteredTwoLengthsDb).
plegma::Options options_of(const Options& opts, float two_lengths_db = 5.0f) {
    plegma::Options options = plegma::Options()
                                  .pinned(false)
                                  .scoring(opts.match, opts.mismatch)
                                  .gaps(opts.gap_open1, opts.gap_ext1, opts.gap_open2, opts.gap_ext2)
                                  // On HG002 HiFi chr12, 18 and 20, 5 dB against plegma's 15: 22 fewer true hets of
                                  // 28,259 for 77 fewer false of 282, the phasing as good.
                                  .twoLengthsDb(two_lengths_db);
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
    remember_inputs("clustered", in, nullptr);
    plegma::Clustering got;
    try {
        got = plegma::clusteredConsensus(in.reads, model(), options_of(opts, kClusteredTwoLengthsDb),
                                         std::min(max_n_cons, 2),
                                         static_cast<float>(opts.min_af));
    } catch (const std::runtime_error& e) {
        if (opts.verbose >= 1) std::fprintf(stderr, "plegma MSA: %s\n", e.what());
        note_msa_failure(std::string("plegma_error: ") + e.what());
        return 0;
    }
    const int n_cons = outputs_of(got, read_ids, read_seqs, read_lens, cons_seqs, clu_n_seqs, clu_read_ids, rows);
    if (n_cons == 0) note_msa_failure("plegma_clusters 0");
    return n_cons;
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
    remember_inputs("phased", in, &read_haps);
    plegma::Clustering got;
    try {
        got = plegma::phasedConsensus(in.reads, read_haps, model(), options_of(opts));
    } catch (const std::runtime_error& e) {
        if (opts.verbose >= 1) std::fprintf(stderr, "plegma MSA: %s\n", e.what());
        note_msa_failure(std::string("plegma_error: ") + e.what());
        return 0;
    }
    std::vector<std::vector<uint8_t>> cons;
    const int n_cons = outputs_of(got, read_ids, read_seqs, read_lens, cons, clu_n_seqs, clu_read_ids, rows);
    if (n_cons != 2) {
        note_msa_failure("plegma_clusters " + std::to_string(n_cons));
        return 0;
    }
    cons_seqs[0] = std::move(cons[0]);
    cons_seqs[1] = std::move(cons[1]);
    return 2;
#else
    (void)opts; (void)read_ids; (void)read_seqs; (void)read_lens; (void)read_covers; (void)read_strands;
    (void)read_haps; (void)cons_seqs; (void)clu_n_seqs; (void)clu_read_ids; (void)rows;
    return 0;
#endif
}

namespace {
thread_local std::string msa_failure;
}

void note_msa_failure(std::string why) { msa_failure = std::move(why); }

std::string take_msa_failure() {
    std::string out;
    out.swap(msa_failure);
    return out;
}

std::string take_plegma_inputs() {
#ifdef PGPHASE_PLEGMA
    std::string out;
    out.swap(last_inputs);
    return out;
#else
    return {};
#endif
}

}  // namespace pgphase_collect
