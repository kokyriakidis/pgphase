#ifndef PGPHASE_ALLELE_GENOTYPE_HPP
#define PGPHASE_ALLELE_GENOTYPE_HPP

#include <array>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace pgphase_collect {

struct MoleculeAlleleCosts {
    std::string molecule;
    std::vector<int> costs;
};

struct DiploidAlleleFit {
    std::optional<std::array<size_t, 2>> alleles;
    std::optional<int64_t> cost;
    std::optional<int64_t> runner_up_cost;
    size_t tied_pairs = 0;
};

struct HeldoutAlleleFit {
    std::string molecule;
    DiploidAlleleFit fit;
    // Full-table allele index; all unselected hypotheses still compete.
    int allele = -1;
    bool same_pair = false;
};

struct DiploidAlleleModel {
    DiploidAlleleFit fit;
    std::vector<HeldoutAlleleFit> heldout;
    std::vector<std::string> conflicting_molecules;
};

/// Minimize sum of per-molecule distance to the nearer allele over all pairs,
/// including homozygous pairs. Ties have no selected genotype. For each molecule,
/// refit with that molecule excluded; neither source HP nor truth enters the fit.
/// Repeated identical molecule costs are idempotent; conflicting costs abstain.
/// Costs and held-out stability are diagnostics, not calibrated confidence.
DiploidAlleleModel fit_diploid_alleles(
    size_t allele_count, const std::vector<MoleculeAlleleCosts>& observations);

} // namespace pgphase_collect

#endif
