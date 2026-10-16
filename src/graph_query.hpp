#ifndef PGPHASE_GRAPH_QUERY_HPP
#define PGPHASE_GRAPH_QUERY_HPP

#include "collect_types.hpp"
#include "graph_sites.hpp"

#include <htslib/hts.h>
#include <htslib/tbx.h>
#include <string>
#include <unordered_map>
#include <vector>

namespace pgphase_collect {

struct GraphReadAllele {
    std::string site_id;
    std::string chrom;
    hts_pos_t pos = 0;
    std::string read_name;
    int allele = kGraphAlleleMissing;
    int mapq = 255;
    bool reverse = false; // true when the read traverses the snarl in reverse-complement
};

// RAII handle for a reference-indexed GAF (pgphase index-gaf; bgzipped + tabix).
// Opened once per thread, reused across chunks to avoid repeated file
// open / tabix index load overhead.
struct IndexedGafHandle {
    htsFile* fp  = nullptr;
    tbx_t*   tbx = nullptr;

    IndexedGafHandle() = default;
    explicit IndexedGafHandle(const std::string& path);
    ~IndexedGafHandle();
    IndexedGafHandle(const IndexedGafHandle&) = delete;
    IndexedGafHandle& operator=(const IndexedGafHandle&) = delete;
    IndexedGafHandle(IndexedGafHandle&& o) noexcept;
    IndexedGafHandle& operator=(IndexedGafHandle&& o) noexcept;
    explicit operator bool() const { return fp && tbx; }
};

// The `#max_span <contig> <bases>` header lines of an indexed GAF (written by
// pgphase index-gaf): the most graph bases any alignment on each contig spans.
// Empty when the GAF has none.
std::unordered_map<std::string, hts_pos_t> read_gaf_max_spans(const std::string& path);

// Per-chunk query using a pre-opened IndexedGafHandle.
// The GAF must have index-gaf's three leading coordinate columns (contig, start, end).
std::vector<GraphReadAllele>
scan_indexed_gaf_chunk(IndexedGafHandle& handle,
                       const std::string& contig,
                       hts_pos_t beg,
                       hts_pos_t end,
                       const GraphSiteCatalogView& catalog,
                       int min_mapq);

// Reports how much graph evidence the walk matcher accepted vs discarded.
void graph_query_report_match_stats();

// Validates that the GAF file is bgzip-compressed and has a .tbi index.
void require_indexed_gaf(const std::string& indexed_gaf_file);

} // namespace pgphase_collect

#endif
