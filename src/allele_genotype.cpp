#include "allele_genotype.hpp"

#include <algorithm>
#include <map>
#include <stdexcept>

namespace pgphase_collect {

static void consider_diploid_pair(DiploidAlleleFit& fit,
        const std::array<size_t, 2>& pair, int64_t cost) {
    if (!fit.cost || cost < *fit.cost) {
        fit.runner_up_cost = fit.cost;
        fit.cost = cost;
        fit.alleles = pair;
        fit.tied_pairs = 1;
    } else {
        if (!fit.runner_up_cost || cost < *fit.runner_up_cost) fit.runner_up_cost = cost;
        if (cost == *fit.cost) {
            ++fit.tied_pairs;
            fit.alleles.reset();
        }
    }
}

DiploidAlleleModel fit_diploid_alleles(
        size_t allele_count, const std::vector<MoleculeAlleleCosts>& observations) {
    DiploidAlleleModel model;
    std::map<std::string, std::optional<std::vector<int>>> by_molecule;
    for (const auto& observation : observations) {
        if (observation.costs.size() != allele_count ||
            std::any_of(observation.costs.begin(), observation.costs.end(), [](int c) { return c < 0; }))
            throw std::invalid_argument("diploid allele costs must be nonnegative and cover every hypothesis");
        const auto [entry, inserted] = by_molecule.emplace(observation.molecule, observation.costs);
        if (!inserted && entry->second != std::optional<std::vector<int>>(observation.costs))
            entry->second.reset();
    }
    std::vector<const std::vector<int>*> costs;
    for (const auto& [name, row] : by_molecule) {
        if (!row) model.conflicting_molecules.push_back(name);
        else {
            model.heldout.push_back({name, {}, -1, false});
            costs.push_back(&*row);
        }
    }
    if (costs.empty() || allele_count == 0) return model;

    for (size_t first = 0; first < allele_count; ++first) {
        for (size_t second = first; second < allele_count; ++second) {
            const std::array<size_t, 2> pair{first, second};
            int64_t total = 0;
            for (const auto* row : costs) total += std::min((*row)[first], (*row)[second]);
            consider_diploid_pair(model.fit, pair, total);
            if (costs.size() == 1) continue;
            for (size_t i = 0; i < costs.size(); ++i) {
                const int own_cost = std::min((*costs[i])[first], (*costs[i])[second]);
                consider_diploid_pair(model.heldout[i].fit, pair, total - own_cost);
            }
        }
    }
    for (size_t i = 0; i < costs.size(); ++i) {
        auto& heldout = model.heldout[i];
        heldout.same_pair = model.fit.alleles && heldout.fit.alleles == model.fit.alleles;
        if (!heldout.fit.alleles) continue;
        const auto& row = *costs[i];
        const auto best = std::min_element(row.begin(), row.end());
        if (std::count(row.begin(), row.end(), *best) != 1) continue;
        const size_t allele = static_cast<size_t>(best - row.begin());
        const auto& pair = *heldout.fit.alleles;
        if (allele == pair[0] || allele == pair[1]) heldout.allele = static_cast<int>(allele);
    }
    return model;
}

} // namespace pgphase_collect
