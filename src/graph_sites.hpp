#ifndef PGPHASE_GRAPH_SITES_HPP
#define PGPHASE_GRAPH_SITES_HPP

// Snarl site catalog: reading the sites of dromos's one-file site catalog
// (site_catalog.hpp), validating allele traversals, and providing the
// GraphSiteCatalog used by graph_query for allele matching and
// graph_bam_adapter for phasing.

#include "collect_types.hpp"

#include <cstdint>
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

// A snarl variation site: one record of the site catalog, as its VCF line
// (dromos view) gives it.
struct GraphSite {
    std::string chrom;          // CHROM (pangenome path name)
    std::string ref_contig;     // suffix after last '#' (e.g. "chr20")
    hts_pos_t pos = 0;          // 1-based POS
    std::string id;              // ID (snarl identifier)
    std::string ref;
    std::vector<std::string> alts;
    std::vector<std::string> allele_traversals;  // AT entries (the unit tests' GFA fixtures)
    std::vector<GraphWalk> allele_walks;          // parsed walks (one per allele)
    hts_pos_t ref_beg = 0;      // linearized ref interval start
    hts_pos_t ref_end = 0;      // linearized ref interval end
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
    size_t data_lines = 0;         ///< records read
    size_t sites_parsed = 0;       ///< records that became a GraphSite
    size_t bad_position = 0;       ///< POS not a positive integer
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

class SiteCatalogReader;

// RAII handle for dromos's site catalog (`dromos catalog`: one file for every
// contig of the graph, its header read once per process and shared; the
// handle holds its thread's reader). Opened once per thread, reused across
// chunks.
struct SiteCatalogHandle {
    std::unique_ptr<SiteCatalogReader> catalog;
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
    explicit operator bool() const { return catalog != nullptr; }
};

/// Load sites overlapping the 1-based inclusive interval [beg, end].
/// The returned catalog is finalized (sorted, validated). Callers querying
/// zero-based APIs such as GAF must convert those coordinates separately.
GraphSiteCatalog load_sites_for_region(SiteCatalogHandle& handle,
                                       const std::string& contig,
                                       hts_pos_t beg, hts_pos_t end);

/// Throws unless `path` is dromos's site catalog (its preamble and header
/// checked).
void require_site_catalog(const std::string& path);

/// The sites of a small graph given as GFA text (S and W lines): every bubble of
/// the reference walk against the others. The unit tests' fixtures; pgphase
/// itself reads the site catalog.
GraphSiteCatalog load_graph_site_catalog_from_gfa_text(const std::string& gfa_text,
                                                       const std::string& reference_sample,
                                                       const std::string& contig);

GraphWalk parse_graph_walk(std::string_view walk);

std::string graph_walk_to_string(const GraphWalk& walk);

GraphWalk reverse_graph_walk(const GraphWalk& walk);

int match_graph_allele_exact(const GraphWalk& read_walk,
                             const std::vector<GraphWalk>& allele_walks,
                             bool* reverse_out = nullptr);

std::string graph_site_between_query(const GraphSite& site);

/// The same contig: "chr20" names any sample's ("CHM13#0#chr20"), but two
/// pangenome path names must also name the same sample.
bool graph_chrom_names_match(const std::string& a, const std::string& b);

/// An allele a site may have: sequence of A, C, G, T and N, the spanning
/// deletion '*', or the missing '.'.
bool graph_allele_is_representable(const std::string& allele);

std::string graph_site_validation_skip_reason(const GraphSite& site);

bool graph_site_is_queryable(const GraphSite& site);

// Canonical site key: VCF ID if present, otherwise "chrom:pos:ref".
inline std::string graph_site_key_str(const GraphSite& site) {
    if (!site.id.empty() && site.id != ".") return site.id;
    return site.chrom + ":" + std::to_string(site.pos) + ":" + site.ref;
}

} // namespace pgphase_collect

#endif
