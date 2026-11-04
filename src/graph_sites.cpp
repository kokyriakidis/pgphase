#include "graph_sites.hpp"
#include "site_catalog.hpp"

#include <algorithm>
#include <charconv>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <unordered_set>

namespace pgphase_collect {

namespace {

std::vector<std::string> split_char(const std::string& value, char delim) {
    std::vector<std::string> out;
    std::string cur;
    std::istringstream in(value);
    while (std::getline(in, cur, delim)) out.push_back(cur);
    if (!value.empty() && value.back() == delim) out.emplace_back();
    return out;
}

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

struct GfaWalkRecord {
    std::string sample;
    std::string hap;
    std::string contig;
    hts_pos_t start = 0;
    hts_pos_t end = 0;
    GraphWalk walk;
};

std::string step_key(const GraphWalkStep& step) {
    std::string out = std::to_string(step.node);
    out.push_back(step.reverse ? '-' : '+');
    return out;
}

hts_pos_t parse_pos_or_zero(const std::string& value) {
    if (value.empty() || value == "*") return 0;
    return static_cast<hts_pos_t>(std::stoll(value));
}

std::vector<GfaWalkRecord>
parse_gfa_walk_records(const std::string& gfa_text,
                       std::unordered_map<std::string, std::string>& node_seqs) {
    std::vector<GfaWalkRecord> walks;
    std::istringstream in(gfa_text);
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> fields = split_char(line, '\t');
        if (fields.empty()) continue;
        if (fields[0] == "S" && fields.size() >= 3) {
            node_seqs[fields[1]] = fields[2];
        } else if (fields[0] == "W" && fields.size() >= 7) {
            try {
                GfaWalkRecord rec;
                rec.sample = fields[1];
                rec.hap = fields[2];
                rec.contig = fields[3];
                rec.start = parse_pos_or_zero(fields[4]);
                rec.end = parse_pos_or_zero(fields[5]);
                rec.walk = parse_graph_walk(fields[6]);
                if (!rec.walk.empty()) walks.push_back(std::move(rec));
            } catch (const std::exception&) {
                continue;
            }
        }
    }
    return walks;
}

size_t find_step_at_or_after(const GraphWalk& walk, const GraphWalkStep& step, size_t beg) {
    for (size_t i = beg; i < walk.size(); ++i) {
        if (walk[i] == step) return i;
    }
    return std::numeric_limits<size_t>::max();
}

GraphWalk slice_walk(const GraphWalk& walk, size_t beg, size_t end_inclusive) {
    if (beg >= walk.size() || end_inclusive >= walk.size() || end_inclusive < beg) return {};
    return GraphWalk(walk.begin() + static_cast<std::ptrdiff_t>(beg),
                     walk.begin() + static_cast<std::ptrdiff_t>(end_inclusive + 1));
}

GraphWalk orient_walk_to_reference(const GraphWalk& walk, const GraphWalk& reference_walk) {
    if (walk.empty() || reference_walk.empty()) return walk;

    std::unordered_set<std::string> ref_steps;
    ref_steps.reserve(reference_walk.size() * 2);
    for (const GraphWalkStep& step : reference_walk) ref_steps.insert(step_key(step));

    const auto count_ref_steps = [&](const GraphWalk& candidate) {
        size_t n = 0;
        for (const GraphWalkStep& step : candidate) {
            if (ref_steps.count(step_key(step))) ++n;
        }
        return n;
    };

    const GraphWalk rev = reverse_graph_walk(walk);
    return count_ref_steps(rev) > count_ref_steps(walk) ? rev : walk;
}

std::vector<hts_pos_t>
reference_step_positions(const GfaWalkRecord& ref,
                         const std::unordered_map<std::string, std::string>& node_seqs) {
    std::vector<hts_pos_t> positions;
    positions.reserve(ref.walk.size());
    hts_pos_t pos = ref.start;
    for (const GraphWalkStep& step : ref.walk) {
        positions.push_back(pos);
        auto it = node_seqs.find(std::to_string(step.node));
        pos += (it == node_seqs.end() || it->second.empty())
                 ? 1
                 : static_cast<hts_pos_t>(it->second.size());
    }
    return positions;
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

GraphSiteCatalog load_graph_site_catalog_from_gfa_text(const std::string& gfa_text,
                                                       const std::string& reference_sample,
                                                       const std::string& contig) {
    std::unordered_map<std::string, std::string> node_seqs;
    const std::vector<GfaWalkRecord> walks = parse_gfa_walk_records(gfa_text, node_seqs);

    int ref_i = -1;
    for (size_t i = 0; i < walks.size(); ++i) {
        if (walks[i].contig == contig && walks[i].sample == reference_sample) {
            ref_i = static_cast<int>(i);
            break;
        }
    }
    if (ref_i < 0) {
        for (size_t i = 0; i < walks.size(); ++i) {
            if (walks[i].contig == contig) {
                ref_i = static_cast<int>(i);
                break;
            }
        }
    }
    if (ref_i < 0) {
        throw std::runtime_error("GFA text did not contain a W line for contig: " + contig);
    }

    const GfaWalkRecord& ref = walks[static_cast<size_t>(ref_i)];
    const std::vector<hts_pos_t> ref_positions = reference_step_positions(ref, node_seqs);

    GraphSiteCatalog catalog;
    std::unordered_map<std::string, size_t> site_by_key;

    const auto add_event = [&](size_t ref_beg, size_t ref_end,
                               const GraphWalk& ref_sub, const GraphWalk& alt_sub) {
        if (ref_sub.size() < 2 || alt_sub.size() < 2 || ref_sub == alt_sub) return;
        const std::string ref_text = graph_walk_to_string(ref_sub);
        const std::string alt_text = graph_walk_to_string(alt_sub);
        const std::string key = contig + "|" + step_key(ref_sub.front()) + "|" +
                                step_key(ref_sub.back()) + "|" + ref_text;

        size_t site_i = 0;
        auto it = site_by_key.find(key);
        if (it == site_by_key.end()) {
            GraphSite site;
            site.chrom = contig;
            site.ref_contig = contig;
            site.pos = ref_beg < ref_positions.size() ? ref_positions[ref_beg] + 1 : 1;
            site.id = "auto_" + std::to_string(catalog.sites.size() + 1) + "_" +
                      step_key(ref_sub.front()) + "_" + step_key(ref_sub.back());
            site.ref = ".";
            site.ref_beg = site.pos;
            site.ref_end = ref_end < ref_positions.size() ? ref_positions[ref_end] + 1 : site.pos;
            site.allele_traversals.push_back(ref_text);
            site.allele_walks.push_back(ref_sub);
            catalog.sites.push_back(std::move(site));
            site_i = catalog.sites.size() - 1;
            site_by_key.emplace(key, site_i);
        } else {
            site_i = it->second;
        }

        GraphSite& site = catalog.sites[site_i];
        if (std::find(site.allele_traversals.begin(),
                      site.allele_traversals.end(),
                      alt_text) == site.allele_traversals.end()) {
            site.allele_traversals.push_back(alt_text);
            site.allele_walks.push_back(alt_sub);
            site.alts.push_back(".");
        }
        (void)ref_end;
    };

    for (size_t walk_i = 0; walk_i < walks.size(); ++walk_i) {
        if (static_cast<int>(walk_i) == ref_i || walks[walk_i].contig != ref.contig) continue;
        const GraphWalk alt = orient_walk_to_reference(walks[walk_i].walk, ref.walk);
        size_t ref_pos = 0;
        size_t alt_scan = 0;
        while (ref_pos + 1 < ref.walk.size()) {
            const size_t alt_left = find_step_at_or_after(alt, ref.walk[ref_pos], alt_scan);
            if (alt_left == std::numeric_limits<size_t>::max()) {
                ++ref_pos;
                continue;
            }
            if (alt_left + 1 < alt.size() && alt[alt_left + 1] == ref.walk[ref_pos + 1]) {
                alt_scan = alt_left + 1;
                ++ref_pos;
                continue;
            }

            size_t ref_right = std::numeric_limits<size_t>::max();
            size_t alt_right = std::numeric_limits<size_t>::max();
            for (size_t j = ref_pos + 1; j < ref.walk.size(); ++j) {
                const size_t found = find_step_at_or_after(alt, ref.walk[j], alt_left + 1);
                if (found != std::numeric_limits<size_t>::max()) {
                    ref_right = j;
                    alt_right = found;
                    break;
                }
            }
            if (ref_right == std::numeric_limits<size_t>::max()) break;

            add_event(ref_pos, ref_right,
                      slice_walk(ref.walk, ref_pos, ref_right),
                      slice_walk(alt, alt_left, alt_right));
            ref_pos = ref_right;
            alt_scan = alt_right;
        }
    }

    std::stable_sort(catalog.sites.begin(), catalog.sites.end(),
                     [](const GraphSite& lhs, const GraphSite& rhs) {
                         if (lhs.chrom != rhs.chrom) return lhs.chrom < rhs.chrom;
                         if (lhs.pos != rhs.pos) return lhs.pos < rhs.pos;
                         return lhs.id < rhs.id;
                     });
    return catalog;
}

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

GraphWalk reverse_graph_walk(const GraphWalk& walk) {
    GraphWalk out;
    out.reserve(walk.size());
    for (auto it = walk.rbegin(); it != walk.rend(); ++it) {
        out.push_back(GraphWalkStep{it->node, !it->reverse});
    }
    return out;
}

int match_graph_allele_exact(const GraphWalk& read_walk,
                             const std::vector<GraphWalk>& allele_walks,
                             bool* reverse_out) {
    int match = kGraphAlleleMissing;
    bool match_reverse = false;
    const GraphWalk reverse_walk = reverse_graph_walk(read_walk);
    for (size_t allele_i = 0; allele_i < allele_walks.size(); ++allele_i) {
        if (read_walk == allele_walks[allele_i]) {
            if (match != kGraphAlleleMissing) return kGraphAlleleAmbiguous;
            match = static_cast<int>(allele_i);
            match_reverse = false;
        } else if (reverse_walk == allele_walks[allele_i]) {
            if (match != kGraphAlleleMissing) return kGraphAlleleAmbiguous;
            match = static_cast<int>(allele_i);
            match_reverse = true;
        }
    }
    if (match >= 0 && reverse_out) *reverse_out = match_reverse;
    return match;
}

// Structural validation: checks that allele walks have consistent boundaries,
// at least 2 unique alleles, and walks of length ≥2 (left + right boundary).
std::string graph_site_validation_skip_reason(const GraphSite& site) {
    if (!site.skip_reason.empty() && !site.eligible) return site.skip_reason;
    if (!site.allele_traversals.empty() &&
        site.allele_walks.size() != site.allele_traversals.size()) {
        return "allele_walk_traversal_count_mismatch";
    }
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

// A site's boundary nodes as "left+:right+" (its first queryable allele walk's first and last steps).
std::string graph_site_between_query(const GraphSite& site) {
    if (!graph_site_is_queryable(site)) return "";
    const GraphWalk* query_walk = nullptr;
    for (const GraphWalk& candidate : site.allele_walks) {
        if (!candidate.empty()) {
            query_walk = &candidate;
            break;
        }
    }
    if (query_walk == nullptr) return "";
    const GraphWalk& walk = *query_walk;
    const GraphWalkStep& left = walk.front();
    const GraphWalkStep& right = walk.back();
    std::string out;
    out += std::to_string(left.node);
    out.push_back(left.reverse ? '-' : '+');
    out.push_back(':');
    out += std::to_string(right.node);
    out.push_back(right.reverse ? '-' : '+');
    return out;
}


// --- SiteCatalogHandle RAII implementation ---

SiteCatalogHandle::SiteCatalogHandle(const std::string& path) {
    refuse_other_sites(path);
    catalog = std::make_unique<SiteCatalogReader>(SiteCatalogFile::open(path));
}

SiteCatalogHandle::~SiteCatalogHandle() = default;
SiteCatalogHandle::SiteCatalogHandle(SiteCatalogHandle&& o) noexcept = default;
SiteCatalogHandle& SiteCatalogHandle::operator=(SiteCatalogHandle&& o) noexcept = default;

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
