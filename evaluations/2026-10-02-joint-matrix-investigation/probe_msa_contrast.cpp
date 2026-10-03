#include "collect_phase_noisy.hpp"
#include <cassert>
#include <iostream>
using namespace pgphase_collect;
static AlnStr alignment(const std::string& target, const std::string& query) {
    assert(target.size() == query.size());
    AlnStr result;
    for (char base : target) result.target_aln.push_back(base == '-' ? 5 : base_to_nt4(base));
    for (char base : query) result.query_aln.push_back(base == '-' ? 5 : base_to_nt4(base));
    result.aln_len = static_cast<int>(target.size());
    result.target_end = result.query_end = result.aln_len - 1;
    return result;
}
int main() {
    const AlnStr four = alignment("ACG----TGC", "ACGTTTTTGC");
    const AlnStr eight = alignment("ACG--------TGC", "ACGTTTTTTTTTGC");
    const AlnStr reference = alignment("ACGTGC", "ACGTGC");
    const std::array<AlnStr, 2> consensuses{four, eight};
    VariantKey key;
    key.pos = 103;
    key.type = VariantType::Insertion;
    key.alt = "TTTT";
    const int selected = call_msa_site_allele({four, four}, key, 100, &consensuses);
    const int other = call_msa_site_allele({eight, eight}, key, 100, &consensuses);
    const int literal_ref = call_msa_site_allele({reference, reference}, key, 100, &consensuses);
    std::cout << "selected=" << selected << " other_verified_consensus=" << other
              << " literal_reference=" << literal_ref << '\n';
    assert(selected == 1 && other == -1 && literal_ref == 0);
}
