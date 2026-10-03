#include <algorithm>
#include <array>
#include <climits>
#include <functional>
#include <iostream>
#include <map>
#include <vector>
struct RecoveryMecRow {
    int fixed_mismatches = 0;
    int fixed_observations = 0;
    std::map<size_t, std::array<int, 2>> variable_mismatches;
};

struct RecoveryMecEffect {
    size_t row = 0;
    int mismatch0 = 0;
    int mismatch1 = 0;
    int observations = 0;
};

struct RecoveryMecOptimum {
    int score = INT_MAX;
    std::vector<int> bits;
};

// Find the exact minimum-error diploid assignment for one fixed orientation of
// the right block. A read pays the smaller Hamming distance to haplotype 1 or
// its complement. Branch-and-bound is exact: its lower bound lets every
// unassigned observation choose any mismatch count, which can only
// underestimate the attainable cost and therefore cannot prune an optimum.
static RecoveryMecOptimum solve_recovery_mec(
        const std::vector<RecoveryMecRow>& rows,
        size_t variable_count) {
    std::vector<std::vector<RecoveryMecEffect>> effects(variable_count);
    std::vector<int> total_observations(rows.size(), 0);
    std::vector<int> mismatches(rows.size(), 0);
    std::vector<int> assigned_observations(rows.size(), 0);
    for (size_t row_i = 0; row_i < rows.size(); ++row_i) {
        const RecoveryMecRow& row = rows[row_i];
        mismatches[row_i] = row.fixed_mismatches;
        assigned_observations[row_i] = row.fixed_observations;
        total_observations[row_i] = row.fixed_observations;
        for (const auto& [variable, costs] : row.variable_mismatches) {
            const int observations = costs[0] + costs[1];
            if (variable >= variable_count || observations == 0) continue;
            effects[variable].push_back(
                RecoveryMecEffect{row_i, costs[0], costs[1], observations});
            total_observations[row_i] += observations;
        }
    }

    std::vector<size_t> order(variable_count);
    for (size_t i = 0; i < variable_count; ++i) order[i] = i;
    std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        const auto support = [&](size_t variable) {
            int n = 0;
            for (const RecoveryMecEffect& effect : effects[variable])
                n += effect.observations;
            return n;
        };
        return support(a) > support(b);
    });

    const auto lower_bound = [&]() {
        int bound = 0;
        for (size_t row_i = 0; row_i < rows.size(); ++row_i) {
            const int remaining = total_observations[row_i] -
                                  assigned_observations[row_i];
            const int low = mismatches[row_i];
            const int high = low + remaining;
            const int total = total_observations[row_i];
            bound += std::min(std::min(low, total - low),
                              std::min(high, total - high));
        }
        return bound;
    };
    const auto score_assignment = [&](const std::vector<int>& assignment) {
        int score = 0;
        for (const RecoveryMecRow& row : rows) {
            int mismatch = row.fixed_mismatches;
            int observations = row.fixed_observations;
            for (const auto& [variable, costs] : row.variable_mismatches) {
                if (variable >= assignment.size()) continue;
                mismatch += costs[static_cast<size_t>(assignment[variable])];
                observations += costs[0] + costs[1];
            }
            score += std::min(mismatch, observations - mismatch);
        }
        return score;
    };

    // A coordinate-descent seed supplies a tight feasible upper bound before
    // exact search. It affects only search order and pruning; branch-and-bound
    // still proves that no lower score exists.
    std::vector<int> seed(variable_count, 0);
    int seed_score = score_assignment(seed);
    bool improved = true;
    while (improved) {
        improved = false;
        for (size_t variable : order) {
            seed[variable] ^= 1;
            const int flipped_score = score_assignment(seed);
            if (flipped_score < seed_score) {
                seed_score = flipped_score;
                improved = true;
            } else {
                seed[variable] ^= 1;
            }
        }
    }

    RecoveryMecOptimum optimum;
    optimum.score = seed_score;
    optimum.bits = seed;
    std::vector<int> bits(variable_count, 0);
    std::function<void(size_t)> search = [&](size_t depth) {
        const int bound = lower_bound();
        // The incumbent is already feasible. Equal-score assignments cannot
        // change source-to-sink parity, so only a strict improvement matters.
        if (bound >= optimum.score) return;
        if (depth == order.size()) {
            optimum.score = bound;
            optimum.bits = bits;
            return;
        }

        const size_t variable = order[depth];
        for (int branch = 0; branch <= 1; ++branch) {
            const int bit = branch == 0 ? seed[variable]
                                        : 1 - seed[variable];
            bits[variable] = bit;
            for (const RecoveryMecEffect& effect : effects[variable]) {
                mismatches[effect.row] +=
                    bit == 0 ? effect.mismatch0 : effect.mismatch1;
                assigned_observations[effect.row] += effect.observations;
            }
            search(depth + 1);
            for (const RecoveryMecEffect& effect : effects[variable]) {
                mismatches[effect.row] -=
                    bit == 0 ? effect.mismatch0 : effect.mismatch1;
                assigned_observations[effect.row] -= effect.observations;
            }
        }
    };
    search(0);
    return optimum;
}


int main() {
    size_t nrows = 0, nvars = 0;
    if (!(std::cin >> nrows >> nvars)) return 1;
    std::vector<RecoveryMecRow> rows(nrows);
    for (auto& row : rows) {
        size_t n = 0;
        std::cin >> row.fixed_mismatches >> row.fixed_observations >> n;
        for (size_t j = 0; j < n; ++j) {
            size_t vi; int c0, c1;
            std::cin >> vi >> c0 >> c1;
            row.variable_mismatches[vi] = {c0, c1};
        }
    }
    auto optimum = solve_recovery_mec(rows, nvars);
    std::cout << optimum.score;
    for (int bit : optimum.bits) std::cout << ' ' << bit;
    std::cout << '\n';
}
