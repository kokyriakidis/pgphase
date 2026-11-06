#pragma once
// Haploid regions of the sample (--haploid-contigs, --par-bed): an XY sample carries one copy of chrX and chrY outside
// their pseudo-autosomal regions, so a site there has no heterozygous genotype. Configured once, before any worker
// starts, and read-only after.

#include <htslib/hts.h>

#include <string>
#include <string_view>
#include <vector>

namespace pgphase_haploid {

/// CONTIGS name the haploid contigs; PAR_BED (may be empty) lists the intervals on them that stay diploid, 0-based
/// half-open. Names match by their part after the last '#' ("CHM13#0#chrX" is "chrX"). Throws std::runtime_error on a
/// BED it cannot read.
void configure(const std::vector<std::string>& contigs, const std::string& par_bed);

bool enabled();

/// True when 1-based POS on CONTIG lies on a haploid contig, outside every PAR interval.
bool at(std::string_view contig, hts_pos_t pos);

}  // namespace pgphase_haploid
