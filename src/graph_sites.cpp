#include "graph_sites.hpp"
#include "site_catalog.hpp"

#include <algorithm>
#include <charconv>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string_view>

namespace pgphase_collect {

namespace {

// Extract the suffix after the last '#' for pangenome path matching.
// "CHM13#0#chr20" → "chr20", "chr20" → "chr20".
std::string chrom_suffix_match_key(const std::string& s) {
    const size_t h = s.rfind('#');
    return h == std::string::npos ? s : s.substr(h + 1);
}

// The sample of a pangenome path name ("CHM13" of "CHM13#0#chr20"); empty for
// a plain contig name.
std::string chrom_sample(const std::string& s) {
    const size_t h = s.find('#');
    return h == std::string::npos ? std::string() : s.substr(0, h);
}

// The same contig: "chr20" matches "CHM13#0#chr20" and "GRCh38#0#chr20", but
// two path names must also name the same sample, as their coordinates differ.
bool chrom_names_match(const std::string& a, const std::string& b) {
    if (a == b) return true;
    if (chrom_suffix_match_key(a) != chrom_suffix_match_key(b)) return false;
    const std::string sa = chrom_sample(a), sb = chrom_sample(b);
    return sa.empty() || sb.empty() || sa == sb;
}

// Sort sites by (chrom, pos, id) and run structural validation on each; a
// site given a skip reason when it was read keeps it.
void finalize_graph_site_catalog_inplace(GraphSiteCatalog& catalog) {
    const auto before = [](const GraphSite& lhs, const GraphSite& rhs) {
        if (lhs.chrom != rhs.chrom) return lhs.chrom < rhs.chrom;
        if (lhs.pos != rhs.pos) return lhs.pos < rhs.pos;
        return lhs.id < rhs.id;
    };
    // A catalog's records come in this order already (POS, then ID).
    if (!std::is_sorted(catalog.sites.begin(), catalog.sites.end(), before))
        std::stable_sort(catalog.sites.begin(), catalog.sites.end(), before);
    for (GraphSite& site : catalog.sites) {
        if (!site.skip_reason.empty()) continue;
        site.skip_reason = graph_site_validation_skip_reason(site);
        site.eligible = site.skip_reason.empty();
    }
    catalog.stats.by_skip_reason.clear();
    for (const GraphSite& site : catalog.sites)
        if (!site.eligible) ++catalog.stats.by_skip_reason[site.skip_reason];
}

/// Alleles this pipeline can turn into a VariantKey: plain sequence, the
/// spanning-deletion marker, or the missing-value dot. A symbolic allele
/// (<DEL>, <INS>) is none of those -- vcf_to_variant_key would treat the angle
/// brackets as sequence and mint a candidate whose ALT is literally "<DEL>".
bool allele_is_representable(const std::string& a) {
    if (a.empty()) return false;
    if (a == "*" || a == ".") return true;
    for (const char ch : a) {
        switch (ch) {
            case 'A': case 'C': case 'G': case 'T': case 'N':
                break;
            default:
                return false;
        }
    }
    return true;
}

} // namespace

// Throws, saying what to use, unless `path` is dromos's site catalog.
static void refuse_other_sites(const std::string& path) {
    if (is_site_catalog_file(path)) return;
    refuse_pgcatalog_sites(path);
    throw std::runtime_error("not dromos's site catalog (dromos catalog: one file, DRMSITES); a sites VCF is not "
                             "read: " + path);
}

void require_site_catalog(const std::string& path) {
    refuse_other_sites(path);
    SiteCatalogFile::open(path);  // its preamble and header checked
}

bool graph_chrom_names_match(const std::string& a, const std::string& b) { return chrom_names_match(a, b); }

bool graph_allele_is_representable(const std::string& allele) { return allele_is_representable(allele); }

GraphWalk parse_graph_walk(std::string_view walk) {
    GraphWalk out;
    size_t i = 0;
    while (i < walk.size()) {
        const char c = walk[i];
        if (c != '>' && c != '<') {
            throw std::runtime_error("invalid graph walk step in: " + std::string(walk));
        }
        const bool reverse = (c == '<');
        ++i;
        const size_t node_beg = i;
        while (i < walk.size() && walk[i] != '>' && walk[i] != '<') ++i;
        if (node_beg == i) {
            throw std::runtime_error("empty node id in graph walk: " + std::string(walk));
        }
        uint64_t node = 0;
        const auto parsed = std::from_chars(walk.data() + node_beg, walk.data() + i, node);
        if (parsed.ec != std::errc() || parsed.ptr != walk.data() + i) {
            throw std::runtime_error("invalid node id in graph walk: " + std::string(walk));
        }
        out.push_back(GraphWalkStep{node, reverse});
    }
    return out;
}

std::string graph_walk_to_string(const GraphWalk& walk) {
    std::string out;
    for (const GraphWalkStep& step : walk) {
        out.push_back(step.reverse ? '<' : '>');
        out += std::to_string(step.node);
    }
    return out;
}

// Structural validation: checks that allele walks have consistent boundaries,
// at least 2 unique alleles, and walks of length ≥2 (left + right boundary).
std::string graph_site_validation_skip_reason(const GraphSite& site) {
    if (!site.skip_reason.empty() && !site.eligible) return site.skip_reason;
    if (site.allele_walks.size() < 2) return "too_few_alleles";

    const GraphWalk* boundary_walk = nullptr;
    bool two_distinct = false;  // two different non-empty walks
    for (const GraphWalk& walk : site.allele_walks) {
        if (walk.empty()) continue;
        if (walk.size() < 2) return "allele_walk_without_boundaries";
        if (boundary_walk != nullptr && !two_distinct && walk != *boundary_walk) two_distinct = true;
        if (boundary_walk == nullptr) {
            boundary_walk = &walk;
            continue;
        }
        if (!(walk.front() == boundary_walk->front()) ||
            !(walk.back() == boundary_walk->back())) {
            return "incompatible_allele_boundaries";
        }
    }
    if (boundary_walk == nullptr) return "no_queryable_allele_walks";
    if (!two_distinct) return "too_few_unique_alleles";
    return "";
}

bool graph_site_is_queryable(const GraphSite& site) {
    return site.eligible && graph_site_validation_skip_reason(site).empty();
}

// --- SiteCatalogHandle RAII implementation ---

SiteCatalogHandle::SiteCatalogHandle(const std::string& path) {
    refuse_other_sites(path);
    catalog = std::make_unique<SiteCatalogReader>(SiteCatalogFile::open(path));
}

SiteCatalogHandle::~SiteCatalogHandle() = default;

GraphSiteCatalog load_sites_for_region(SiteCatalogHandle& handle,
                                       const std::string& contig,
                                       hts_pos_t beg, hts_pos_t end) {
    GraphSiteCatalog catalog;
    if (!handle || contig.empty() || end <= beg) return catalog;

    hts_pos_t max_span = std::numeric_limits<hts_pos_t>::max();
    const auto span = handle.max_span.find(chrom_suffix_match_key(contig));
    if (span != handle.max_span.end()) max_span = span->second;
    handle.catalog->append_sites(contig, beg, end, max_span, catalog);
    finalize_graph_site_catalog_inplace(catalog);
    return catalog;
}

} // namespace pgphase_collect
