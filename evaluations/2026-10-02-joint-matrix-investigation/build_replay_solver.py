from pathlib import Path
s=Path('src/collect_phase.cpp').read_text();a=s.index('struct RecoveryMecRow {');b=s.index('// Solve one adjacent recovery edge',a)
headers='''#include <algorithm>
#include <array>
#include <climits>
#include <functional>
#include <iostream>
#include <map>
#include <vector>
'''
main=r'''
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
'''
Path('test_data/tmp_gap_next29/mec_replay.cpp').write_text(headers+s[a:b]+main)
