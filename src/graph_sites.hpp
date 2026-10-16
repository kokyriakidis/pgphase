#ifndef PGPHASE_GRAPH_SITES_HPP
#define PGPHASE_GRAPH_SITES_HPP

// Snarl site catalog: reading the sites of a catalog (pgcat .pgsites files, or a
// sites-only VCF in the format of `vg deconstruct -a`), validating allele
// traversals, and providing the GraphSiteCatalog used by graph_query for allele
// matching and graph_bam_adapter for phasing.

#include "collect_types.hpp"

#include <htslib/tbx.h>
#include <cstdint>
#include <iosfwd>
#include <map>
#include <memory>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace pgphase_collect {

/// The most nodes an allele walk may visit for reads to be matched to its site
/// (the compact index stores walk lengths in 16 bits). A site with a longer
/// walk is never matched, so loaders skip it.
constexpr size_t kMaxSiteWalkSteps = 65535;

// One step in an allele traversal walk: a graph node ID and its orientation.
struct GraphWalkStep {
    uint64_t node = 0;
    bool reverse = false;

    bool operator==(const GraphWalkStep& other) const {
        return node == other.node && reverse == other.reverse;
    }
};

using GraphWalk = std::vector<GraphWalkStep>;

constexpr int kGraphAlleleMissing = -1;    // read doesn't span the site
constexpr int kGraphAlleleAmbiguous = -2;  // read matches multiple alleles

// A snarl variation site parsed from a VCF record.
struct GraphSite {
    std::string chrom;          // VCF CHROM (pangenome path name)
    std::string ref_contig;     // suffix after last '#' (e.g. "chr20")
    hts_pos_t pos = 0;          // VCF 1-based POS
    std::string id;              // VCF ID field (snarl identifier)
    std::string ref;
    std::vector<std::string> alts;
    std::vector<std::string> allele_traversals;  // raw AT field entries
    std::vector<GraphWalk> allele_walks;          // parsed walks (one per allele)
    int level = -1;              // snarl nesting level (LV info field)
    std::string parent;          // parent snarl site key (PS info field)
    std::string root;            // root snarl (unused, reserved)
    hts_pos_t ref_beg = 0;      // linearized ref interval start
    hts_pos_t ref_end = 0;      // linearized ref interval end
    std::vector<int> conditional_parent_alleles;  // PA field: parent alleles that gate this child
    bool has_spanning_deletion = false;  // true if any ALT is '*'
    bool bam_alt_deletion_no_ref = false;
    bool bam_low_fraction_snp = false;
    bool bam_homozygous_alt = false;  // physical BAM evidence contradicts graph heterozygosity
    bool eligible = true;        // false if validation failed
    std::string skip_reason;     // why the site was marked ineligible

    hts_pos_t order_pos() const { return pos; }
};

struct GraphSiteCatalogView;

/// What a load did, including everything it chose not to keep.
///
/// A loader that drops a record silently is indistinguishable from one that
/// loads a file correctly, so every path that discards something counts it
/// here. `by_skip_reason` covers sites that parsed but failed structural
/// validation; the other counters cover records that never became sites.
struct GraphSiteLoadStats {
    size_t data_lines = 0;         ///< non-header lines seen
    size_t sites_parsed = 0;       ///< lines that became a GraphSite
    size_t short_line = 0;         ///< fewer than the 8 mandatory VCF columns
    size_t bad_position = 0;       ///< POS or END not a positive integer
    size_t beyond_read_span = 0;   ///< every allele longer than any alignment on the contig: skipped
    size_t walk_too_long = 0;      ///< an allele walk over kMaxSiteWalkSteps nodes: never matched, skipped
    size_t unsupported_allele = 0; ///< symbolic <DEL>, or a character outside ACGTN*.
    std::map<std::string, size_t> by_skip_reason;  ///< ineligible sites, by reason

    size_t eligible() const {
        size_t bad = 0;
        for (const auto& kv : by_skip_reason) bad += kv.second;
        return sites_parsed >= bad ? sites_parsed - bad : 0;
    }
    /// One line for a log or a test failure message.
};

struct GraphSiteCatalog {
    std::vector<GraphSite> sites;
    GraphSiteLoadStats stats;

    // Return a view covering all sites (no copy).
    GraphSiteCatalogView view_all() const;
};

// Lightweight non-owning view into a subset of a GraphSiteCatalog's sites.
// Avoids deep-copying GraphSite objects when building per-chunk catalogs.
// The referenced catalog must outlive the view.
struct GraphSiteCatalogView {
    const std::vector<GraphSite>* source = nullptr;
    std::vector<size_t> indices;

    size_t size() const { return indices.size(); }
    bool empty() const { return indices.empty(); }
    const GraphSite& operator[](size_t i) const { return (*source)[indices[i]]; }
};

inline GraphSiteCatalogView GraphSiteCatalog::view_all() const {
    GraphSiteCatalogView v;
    v.source = &sites;
    v.indices.resize(sites.size());
    for (size_t i = 0; i < sites.size(); ++i) v.indices[i] = i;
    return v;
}

// When `filters` is non-empty only sites overlapping at least one filter are loaded.
// For bgzipped VCFs with a .tbi/.csi index, tabix is used to seek directly to each
// region; otherwise the file is streamed and non-overlapping lines are skipped.
GraphSiteCatalog load_graph_site_catalog_from_vcf(
    const std::string& path,
    const std::vector<RegionFilter>& filters = {},
    bool keep_allele_traversal_strings = true);

struct PgSitesCatalog;

// RAII handle for a site catalog: a tabix-indexed sites VCF, a .pgsites file
// (the same VCF lines in independently zstd-compressed position blocks, written
// by pgcat), or a folder of <contig>.pgsites files (pgcat build's output), each
// opened when its contig is first queried. Opened once per thread, reused
// across chunks to avoid repeated file open / index load.
struct SiteCatalogHandle {
    htsFile* fp  = nullptr;
    tbx_t*   tbx = nullptr;
    std::unique_ptr<PgSitesCatalog> pgsites;
    // The most graph bases any alignment on each contig spans (bare contig
    // name), from the GAF's #max_span lines. A site whose every allele is
    // longer cannot be traversed whole by any read, so it is not loaded; for a
    // contig not listed, every site is.
    std::unordered_map<std::string, hts_pos_t> max_span;

    SiteCatalogHandle() = default;
    explicit SiteCatalogHandle(const std::string& path);
    ~SiteCatalogHandle();
    SiteCatalogHandle(const SiteCatalogHandle&) = delete;
    SiteCatalogHandle& operator=(const SiteCatalogHandle&) = delete;
    SiteCatalogHandle(SiteCatalogHandle&& o) noexcept;
    SiteCatalogHandle& operator=(SiteCatalogHandle&& o) noexcept;
    explicit operator bool() const { return (fp && tbx) || pgsites; }
};

/// Load sites overlapping the 1-based inclusive VCF interval [beg, end].
/// The returned catalog is finalized (sorted, validated). Callers querying
/// zero-based APIs such as GAF must convert those coordinates separately.
GraphSiteCatalog load_sites_for_region(SiteCatalogHandle& handle,
                                       const std::string& contig,
                                       hts_pos_t beg, hts_pos_t end);

/// Throws unless `path` is a usable site catalog: a bgzipped VCF with a tabix
/// index, a .pgsites file, or a folder holding .pgsites files.
void require_site_catalog(const std::string& path);

GraphSiteCatalog load_graph_site_catalog_from_gfa_text(const std::string& gfa_text,
                                                       const std::string& reference_sample,
                                                       const std::string& contig);

void write_graph_site_catalog_tsv(std::ostream& out, const GraphSiteCatalog& catalog);


GraphWalk parse_graph_walk(std::string_view walk);

std::string graph_walk_to_string(const GraphWalk& walk);

GraphWalk reverse_graph_walk(const GraphWalk& walk);

int match_graph_allele_exact(const GraphWalk& read_walk,
                             const std::vector<GraphWalk>& allele_walks,
                             bool* reverse_out = nullptr);

std::string graph_site_between_query(const GraphSite& site);

std::string graph_site_validation_skip_reason(const GraphSite& site);

bool graph_site_is_queryable(const GraphSite& site);

// Canonical site key: VCF ID if present, otherwise "chrom:pos:ref".
inline std::string graph_site_key_str(const GraphSite& site) {
    if (!site.id.empty() && site.id != ".") return site.id;
    return site.chrom + ":" + std::to_string(site.pos) + ":" + site.ref;
}

} // namespace pgphase_collect

#endif
