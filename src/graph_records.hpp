#ifndef PGPHASE_GRAPH_RECORDS_HPP
#define PGPHASE_GRAPH_RECORDS_HPP

// The graph command's records (graph_records.cpp), apart from the run that produces the chunks (graph_collect.cpp).

#include "bam_io_types.hpp"
#include "graph_bam_adapter.hpp"

#include <htslib/sam.h>

#include <string>
#include <unordered_map>
#include <vector>

namespace pgphase_collect {

/// A batch's records: each chunk's rows as VCF records for the shared writer -- every surviving allele of a row
/// trimmed to minimal VCF form and keyed as the BAM path keys it, with the row's genotype, phase set and provenance,
/// and place_unlinked_hets' placements applied. Rows written no record (GraphSiteMeta::no_record, an empty REF, an
/// injected site without a het consensus) are skipped, each traced with its reason.
CandidateTable graph_chunks_to_candidate_table(const std::vector<GraphChunkBuildResult>& graph_chunks,
                                              const std::unordered_map<std::string, int>& contig_to_tid,
                                              const Options& opts, ReferenceCache& reference,
                                              const bam_hdr_t* header);

}  // namespace pgphase_collect

#endif  // PGPHASE_GRAPH_RECORDS_HPP
