#ifndef PGPHASE_OUTPUT_RECONCILE_HPP
#define PGPHASE_OUTPUT_RECONCILE_HPP

// Post-phasing record reconciliation ("output" trace stage): merges adjacent deletion+insertion pairs on one
// haplotype, drops contradicted catalog walks, and resolves co-located alleles claiming the same haplotype.

#include "collect_types.hpp"

#include <functional>

namespace pgphase_collect {

/// Drop records that put two DIFFERENT alleles on the SAME haplotype at one
/// position (the POS each is written at, vcf_record_pos). Contained alleles lose to the complete one (a SNP against the
/// insertion containing it); when neither contains the other, the better
/// supported by allele depth wins, ties going to the shorter allele.
/// Make two co-located records claiming one haplotype complementary, one allele
/// per haplotype, the way longcallD's structure does by construction. Returns
/// the number re-oriented.
size_t make_colocated_alleles_complementary(CandidateTable& records, int min_alt_support);

void drop_conflicting_haplotype_alleles(CandidateTable& records);

/// A deletion and an insertion side by side on one haplotype (the insertion right after the deleted bases, or right
/// before them) written as one record: the deletion's, its ALT the inserted bases. Returns the number merged.
size_t merge_adjacent_deletion_insertion(CandidateTable& records);

/// A long deletion or replacement (a catalog walk) a record inside it contradicts on one haplotype of a phase set -- it
/// changes a base the walk deletes or keeps -- is dropped. Returns the number dropped.
size_t drop_contradicted_walks(CandidateTable& records, const std::function<char(int, hts_pos_t)>& ref_base);


} // namespace pgphase_collect

#endif // PGPHASE_OUTPUT_RECONCILE_HPP
