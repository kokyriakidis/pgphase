#include "graph_sites.hpp"

#include <algorithm>
#include <charconv>
#include <cctype>
#include <cerrno>
#include <cstring>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <fstream>
#include <iostream>
#include <mutex>
#include <unordered_set>

#include <dirent.h>
#include <sys/stat.h>

#include <zstd.h>

#include <htslib/hts.h>
#include <htslib/kstring.h>
#include <htslib/tbx.h>

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

// Calls f with each `delim`-separated item of s, as std::getline splitting
// does: nothing for an empty s, and a final empty item after a trailing delim.
template <class F>
void for_each_item(std::string_view s, char delim, F&& f) {
    if (s.empty()) return;
    size_t start = 0;
    while (true) {
        const size_t d = s.find(delim, start);
        if (d == std::string_view::npos) {
            f(s.substr(start));
            return;
        }
        f(s.substr(start, d - start));
        start = d + 1;
    }
}

// A VCF INFO field as key/value views into the line. A key given twice keeps
// its first value; a key without '=' has an empty value.
struct InfoView {
    std::vector<std::pair<std::string_view, std::string_view>> items;

    explicit InfoView(std::string_view field) {
        if (field.empty() || field == ".") return;
        for_each_item(field, ';', [&](std::string_view item) {
            if (item.empty()) return;
            const size_t eq = item.find('=');
            if (eq == std::string_view::npos) items.emplace_back(item, std::string_view());
            else items.emplace_back(item.substr(0, eq), item.substr(eq + 1));
        });
    }

    // The value of the first of `keys` present; empty if none is.
    std::string_view first_of(std::initializer_list<std::string_view> keys) const {
        for (std::string_view key : keys)
            for (const auto& kv : items)
                if (kv.first == key) return kv.second;
        return {};
    }

    bool has(std::string_view key) const {
        for (const auto& kv : items)
            if (kv.first == key) return true;
        return false;
    }
};

int parse_optional_int(const InfoView& info, std::string_view key) {
    const std::string_view v = info.first_of({key});
    if (v.empty()) return -1;
    return std::stoi(std::string(v));
}

std::vector<int> parse_optional_int_list(const InfoView& info, std::initializer_list<std::string_view> keys) {
    const std::string_view value = info.first_of(keys);
    std::vector<int> out;
    if (value.empty() || value == ".") return out;
    for_each_item(value, ',', [&](std::string_view field) {
        if (field.empty() || field == ".") return;
        out.push_back(std::stoi(std::string(field)));
    });
    return out;
}

bool contains_spanning_deletion(const std::vector<std::string>& alts) {
    return std::find(alts.begin(), alts.end(), "*") != alts.end();
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
bool chrom_names_match_vcf(const std::string& a, const std::string& b) {
    if (a == b) return true;
    if (chrom_suffix_match_key(a) != chrom_suffix_match_key(b)) return false;
    const std::string sa = chrom_sample(a), sb = chrom_sample(b);
    return sa.empty() || sb.empty() || sa == sb;
}

// Sort sites by (chrom, pos, id) and run structural validation on each.
// Parse-time skip reasons (e.g. malformed_allele_traversal) are preserved
// because the validator can't detect them from the stored empty walk.
void finalize_graph_site_catalog_inplace(GraphSiteCatalog& catalog) {
    std::stable_sort(catalog.sites.begin(), catalog.sites.end(),
                     [](const GraphSite& lhs, const GraphSite& rhs) {
                         if (lhs.chrom != rhs.chrom) return lhs.chrom < rhs.chrom;
                         if (lhs.pos != rhs.pos) return lhs.pos < rhs.pos;
                         return lhs.id < rhs.id;
                     });
    for (GraphSite& site : catalog.sites) {
        // Preserve parse-time skip reasons (e.g. malformed_allele_traversal)
        // that the structural validator cannot detect because the malformed
        // walk was stored as an empty entry.
        if (!site.skip_reason.empty()) continue;
        site.skip_reason = graph_site_validation_skip_reason(site);
        site.eligible = site.skip_reason.empty();
    }
    catalog.stats.by_skip_reason.clear();
    for (const GraphSite& site : catalog.sites)
        if (!site.eligible) ++catalog.stats.by_skip_reason[site.skip_reason];
}

/// A positive integer, or nothing. std::stoll throws on malformed input, and a
/// throw out of the parser aborts the whole run over one bad line with no
/// indication of which; a catalog is external data and a single corrupt record
/// should cost that record, not the analysis.
bool parse_positive(const std::string& text, hts_pos_t& out) {
    if (text.empty()) return false;
    errno = 0;
    char* end = nullptr;
    const long long v = std::strtoll(text.c_str(), &end, 10);
    if (errno != 0 || end == text.c_str() || *end != '\0' || v <= 0) return false;
    out = static_cast<hts_pos_t>(v);
    return true;
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

void upper_in_place(std::string& s) {
    for (char& ch : s) ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
}

// Parse one VCF data line into a GraphSite and append it to the catalog.
// Extracts CHROM, POS, ID, REF, ALT, and INFO fields (AT, LV, PS, PA, END).
void append_graph_site_from_vcf_data_line(const char* buf,
                                          size_t len,
                                          GraphSiteCatalog& catalog,
                                          bool keep_allele_traversal_strings) {
    if (len == 0 || buf[0] == '#') return;
    ++catalog.stats.data_lines;
    // Fields are views into the line; only what the site keeps is copied.
    std::vector<std::string_view> fields;
    fields.reserve(10);
    for_each_item(std::string_view(buf, len), '\t', [&](std::string_view f) { fields.push_back(f); });
    if (fields.size() < 8) { ++catalog.stats.short_line; return; }

    GraphSite site;
    site.chrom = std::string(fields[0]);
    const auto hash_pos = fields[0].rfind('#');
    site.ref_contig = std::string((hash_pos != std::string_view::npos)
                                  ? fields[0].substr(hash_pos + 1)
                                  : fields[0]);
    if (!parse_positive(std::string(fields[1]), site.pos)) { ++catalog.stats.bad_position; return; }
    site.id = std::string(fields[2]);
    site.ref = std::string(fields[3]);
    if (fields[4] != ".")
        for_each_item(fields[4], ',', [&](std::string_view a) { site.alts.emplace_back(a); });
    // Case is not meaning: a soft-masked reference writes acgt, and an allele
    // that differs from a candidate's only in case would never match it.
    upper_in_place(site.ref);
    for (std::string& a : site.alts) upper_in_place(a);
    if (!allele_is_representable(site.ref)) { ++catalog.stats.unsupported_allele; return; }
    for (const std::string& a : site.alts) {
        if (allele_is_representable(a)) continue;
        ++catalog.stats.unsupported_allele;
        return;
    }
    const InfoView info(fields[7]);
    // AT (or UT, LVAT from older vg) walks, parsed straight from the line.
    // Spanning deletions ('*') and missing entries ('.') are empty walks; an
    // entry that does not parse marks the site ineligible.
    bool malformed_walk = false;
    for_each_item(info.first_of({"AT", "UT", "LVAT"}), ',', [&](std::string_view t) {
        if (keep_allele_traversal_strings) site.allele_traversals.emplace_back(t);
        if (t == "*" || t == ".") {
            site.allele_walks.emplace_back();
            return;
        }
        try {
            site.allele_walks.push_back(parse_graph_walk(t));
        } catch (const std::exception&) {
            malformed_walk = true;
            site.allele_walks.emplace_back();
        }
    });
    site.level = parse_optional_int(info, "LV");
    site.parent = std::string(info.first_of({"PS", "PARENT", "Parent"}));
    site.root = std::string(info.first_of({"RS", "ROOT", "Root"}));
    const std::string rc_value(info.first_of({"RC", "REF_CONTIG", "RefContig"}));
    if (!rc_value.empty()) {
        const auto hp = rc_value.rfind('#');
        site.ref_contig = (hp != std::string::npos) ? rc_value.substr(hp + 1) : rc_value;
    }
    site.ref_beg = site.pos;
    site.ref_end = site.pos;
    const std::string end_value(info.first_of({"END"}));
    if (!end_value.empty() && !parse_positive(end_value, site.ref_end)) {
        ++catalog.stats.bad_position;
        return;
    }
    site.conditional_parent_alleles =
        parse_optional_int_list(info, {"PA", "PARENT_ALLELE", "PARENT_ALLELES"});
    site.has_spanning_deletion = contains_spanning_deletion(site.alts);
    if (malformed_walk) {
        site.eligible = false;
        site.skip_reason = "malformed_allele_traversal";
    } else {
        site.skip_reason = graph_site_validation_skip_reason(site);
        site.eligible = site.skip_reason.empty();
    }
    ++catalog.stats.sites_parsed;
    catalog.sites.push_back(std::move(site));
}

// The fewest bases any allele of a VCF data line spans, the padding base not
// counted (alleles of unequal lengths carry one), as pgcat computes it; 0 when
// the line has no REF/ALT fields.
hts_pos_t shortest_allele(const char* line, size_t len) {
    const char* end = line + len;
    const char* p = line;
    for (int tabs = 0; tabs < 3; ++p) {          // to REF, the 4th field
        if (p >= end) return 0;
        if (*p == '\t') ++tabs;
    }
    hts_pos_t shortest = -1, first = -1;
    bool equal = true;
    int fields = 0;                               // REF, then ALT
    for (const char* a = p; p <= end; ++p) {
        if (p != end && *p != '\t' && *p != ',') continue;
        const hts_pos_t l = p - a;
        if (first < 0) first = l;
        equal = equal && l == first;
        if (shortest < 0 || l < shortest) shortest = l;
        a = p + 1;
        if (p == end || (*p == '\t' && ++fields == 2)) break;
    }
    if (shortest < 0) return 0;
    return equal ? shortest : shortest - 1;
}

// The most nodes any allele walk in the line's AT field visits.
size_t longest_walk(const char* line, size_t len) {
    const char* end = line + len;
    const char* at = nullptr;
    for (const char* p = line; p + 3 <= end; ++p) {
        if ((p == line || p[-1] == '\t' || p[-1] == ';') && p[0] == 'A' && p[1] == 'T' && p[2] == '=') {
            at = p + 3;
            break;
        }
    }
    if (at == nullptr) return 0;
    size_t longest = 0, steps = 0;
    for (const char* p = at; p < end && *p != ';' && *p != '\t'; ++p) {
        if (*p == '>' || *p == '<') ++steps;
        else if (*p == ',') { longest = std::max(longest, steps); steps = 0; }
    }
    return std::max(longest, steps);
}

// Whether no read can be matched to the record's site, and so it need not be
// loaded: every allele spans more graph bases than any alignment, or some
// allele walk is longer than the matcher takes. Counted in the load stats.
bool unmatchable_record(const char* line, size_t len, hts_pos_t max_span, GraphSiteLoadStats& stats) {
    if (shortest_allele(line, len) > max_span) {
        ++stats.beyond_read_span;
        return true;
    }
    // A walk of more than kMaxSiteWalkSteps nodes takes over two bytes a node.
    if (len > 2 * kMaxSiteWalkSteps && longest_walk(line, len) > kMaxSiteWalkSteps) {
        ++stats.walk_too_long;
        return true;
    }
    return false;
}

// Load sites from a tabix-indexed VCF, restricted to the given regions.
// Handles pangenome path name mismatches (e.g. filter asks for "chr20" but
// VCF uses "CHM13#0#chr20") via suffix matching.
void append_graph_sites_tabix_filtered(htsFile* fp,
                                       tbx_t* tbx,
                                       const std::vector<RegionFilter>& filters,
                                       hts_pos_t max_span,
                                       GraphSiteCatalog& catalog,
                                       bool keep_allele_traversal_strings) {
    struct HtsItrDeleter {
        void operator()(hts_itr_t* p) const {
            if (p) hts_itr_destroy(p);
        }
    };

    int n_vcf_seqs = 0;
    const char** vcf_seqs = tbx_seqnames(tbx, &n_vcf_seqs);
    std::unordered_map<std::string, std::string> filter_to_vcf_chrom;
    for (const RegionFilter& filter : filters) {
        if (!filter.enabled || filter_to_vcf_chrom.count(filter.chrom)) continue;
        if (tbx_name2id(tbx, filter.chrom.c_str()) >= 0) {
            filter_to_vcf_chrom[filter.chrom] = filter.chrom;
        } else {
            for (int si = 0; si < n_vcf_seqs; ++si) {
                if (chrom_names_match_vcf(filter.chrom, std::string(vcf_seqs[si]))) {
                    filter_to_vcf_chrom[filter.chrom] = vcf_seqs[si];
                    break;
                }
            }
        }
    }
    free(vcf_seqs);

    kstring_t str = KS_INITIALIZE;
    struct KsGuard {
        kstring_t* s;
        ~KsGuard() {
            if (s) ks_free(s);
        }
    } ks_guard{&str};

    auto make_rstr = [](const std::string& chrom, hts_pos_t beg, hts_pos_t end) {
        std::string s = chrom;
        if (beg > 0) {
            s += ":" + std::to_string(beg);
            if (end > 0) s += "-" + std::to_string(end);
        }
        return s;
    };

    for (const RegionFilter& filter : filters) {
        if (!filter.enabled) continue;
        auto it = filter_to_vcf_chrom.find(filter.chrom);
        if (it == filter_to_vcf_chrom.end()) continue;
        std::unique_ptr<hts_itr_t, HtsItrDeleter> itr(
            tbx_itr_querys(tbx,
                           make_rstr(it->second, filter.beg, filter.end).c_str()));
        if (!itr) continue;
        while (tbx_itr_next(fp, tbx, itr.get(), &str) >= 0) {
            if (str.l == 0) continue;
            if (!unmatchable_record(str.s, str.l, max_span, catalog.stats))
                append_graph_site_from_vcf_data_line(str.s, str.l, catalog,
                                                     keep_allele_traversal_strings);
        }
    }
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

GraphSiteCatalog load_graph_site_catalog_from_vcf(
    const std::string& path,
    const std::vector<RegionFilter>& filters,
    bool keep_allele_traversal_strings)
{
    GraphSiteCatalog catalog;

    bool has_filters = false;
    for (const RegionFilter& f : filters) {
        if (f.enabled) { has_filters = true; break; }
    }

    struct TbxDeleter { void operator()(tbx_t* p) const { tbx_destroy(p); } };
    struct HtsFileDeleter { void operator()(htsFile* p) const { hts_close(p); } };

    if (has_filters) {
        // Region-filtered: require tabix index for efficient random access.
        std::unique_ptr<tbx_t, TbxDeleter> tbx(tbx_index_load(path.c_str()));
        if (!tbx)
            throw std::runtime_error(
                "sites VCF must be bgzip-compressed with a tabix index: " + path);
        std::unique_ptr<htsFile, HtsFileDeleter> fp(hts_open(path.c_str(), "r"));
        if (!fp) throw std::runtime_error("failed to open graph site VCF: " + path);
        append_graph_sites_tabix_filtered(fp.get(), tbx.get(), filters,
                                           std::numeric_limits<hts_pos_t>::max(), catalog,
                                           keep_allele_traversal_strings);
    } else {
        // No region filters: stream the entire file.
        std::unique_ptr<htsFile, HtsFileDeleter> fp(hts_open(path.c_str(), "r"));
        if (!fp) throw std::runtime_error("failed to open graph site VCF: " + path);
        kstring_t str = KS_INITIALIZE;
        while (hts_getline(fp.get(), '\n', &str) >= 0) {
            if (str.l == 0 || str.s[0] == '#') continue;
            append_graph_site_from_vcf_data_line(str.s, str.l, catalog,
                                                 keep_allele_traversal_strings);
        }
        ks_free(&str);
    }

    finalize_graph_site_catalog_inplace(catalog);
    return catalog;
}

// --- .pgsites catalogs ---------------------------------------------------------
//
// "PGSITES2", u64 header length + zstd(header text), zstd blocks of VCF lines,
// then the index (u64 count; per block u64 first_pos, last_pos, max_end, offset,
// compressed length, raw length, records, min_allele, max_walk) and the index's
// u64 offset last. min_allele is the fewest bases any allele of the block's
// records spans, max_walk the most nodes any of their allele walks visits: a
// block whose records no read can be matched to is skipped unread.
// max_end is the largest POS + len(REF) - 1 in the block: the record span tabix
// uses for VCF, so a region query selects the same records either way. One file
// holds one contig, named by the header's ##contig line.

struct PgSitesFile {
    struct Block { uint64_t first_pos, last_pos, max_end, offset, comp_len, raw_len, records, min_allele, max_walk; };
    std::ifstream in;
    std::string chrom;
    std::vector<Block> blocks;
};

// One .pgsites file, or a folder of them opened as their contigs are queried.
struct PgSitesCatalog {
    std::string dir;  // empty for a single file
    std::vector<std::unique_ptr<PgSitesFile>> files;
    std::unordered_set<std::string> absent;  // contigs the folder has no file for
};

namespace {

constexpr char kPgSitesMagic[8] = {'P', 'G', 'S', 'I', 'T', 'E', 'S', '2'};
constexpr const char* kPgSitesExt = ".pgsites";

bool is_directory(const std::string& path) {
    struct stat st;
    return stat(path.c_str(), &st) == 0 && S_ISDIR(st.st_mode);
}

// A .pgsites file of any version: the version is checked when it is opened.
bool is_pgsites_file(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    char magic[8] = {};
    return in.read(magic, 8) && std::memcmp(magic, kPgSitesMagic, 7) == 0;
}

bool has_pgsites_file(const std::string& dir) {
    DIR* d = opendir(dir.c_str());
    if (d == nullptr) return false;
    bool found = false;
    while (const dirent* e = readdir(d)) {
        const std::string name = e->d_name;
        if (name.size() > std::strlen(kPgSitesExt) &&
            name.compare(name.size() - std::strlen(kPgSitesExt), std::string::npos, kPgSitesExt) == 0 &&
            is_pgsites_file(dir + "/" + name)) {
            found = true;
            break;
        }
    }
    closedir(d);
    return found;
}

uint64_t read_u64(std::istream& in) {
    unsigned char b[8];
    if (!in.read(reinterpret_cast<char*>(b), 8)) throw std::runtime_error("truncated .pgsites file");
    uint64_t v = 0;
    for (int i = 7; i >= 0; --i) v = (v << 8) | b[i];
    return v;
}

// Blocks holding a giant record are written with a long window (log 31), above
// zstd's default decoding limit, so the limit is raised for every block.
std::string zstd_decompress(const std::string& comp, size_t raw_len) {
    struct Ctx {
        ZSTD_DCtx* c = ZSTD_createDCtx();
        Ctx() { ZSTD_DCtx_setParameter(c, ZSTD_d_windowLogMax, 31); }
        ~Ctx() { ZSTD_freeDCtx(c); }
    };
    thread_local Ctx ctx;
    std::string out(raw_len, '\0');
    const size_t n = ZSTD_decompressDCtx(ctx.c, out.data(), raw_len, comp.data(), comp.size());
    if (ZSTD_isError(n) || n != raw_len) throw std::runtime_error("corrupt .pgsites block");
    return out;
}

std::unique_ptr<PgSitesFile> open_pgsites(const std::string& path) {
    auto f = std::make_unique<PgSitesFile>();
    f->in.open(path, std::ios::binary);
    if (!f->in) throw std::runtime_error("failed to open site catalog: " + path);
    char magic[8] = {};
    if (!f->in.read(magic, 8) || std::memcmp(magic, kPgSitesMagic, 8) != 0)
        throw std::runtime_error(path + ": an older .pgsites format; rebuild it with pgphase build-snarl-catalog");
    const uint64_t header_len = read_u64(f->in);
    std::string comp(header_len, '\0');
    f->in.read(comp.data(), static_cast<std::streamsize>(header_len));
    const unsigned long long raw = ZSTD_getFrameContentSize(comp.data(), comp.size());
    if (raw == ZSTD_CONTENTSIZE_ERROR || raw == ZSTD_CONTENTSIZE_UNKNOWN)
        throw std::runtime_error("corrupt .pgsites header: " + path);
    const std::string header = zstd_decompress(comp, raw);
    const size_t at = header.find("##contig=<ID=");
    if (at == std::string::npos) throw std::runtime_error(".pgsites header has no ##contig: " + path);
    const size_t b = at + 13;
    f->chrom = header.substr(b, header.find_first_of(",>", b) - b);
    f->in.seekg(-8, std::ios::end);
    const uint64_t index_offset = read_u64(f->in);
    f->in.seekg(static_cast<std::streamoff>(index_offset));
    const uint64_t n = read_u64(f->in);
    f->blocks.resize(n);
    for (auto& blk : f->blocks) {
        blk.first_pos = read_u64(f->in);
        blk.last_pos = read_u64(f->in);
        blk.max_end = read_u64(f->in);
        blk.offset = read_u64(f->in);
        blk.comp_len = read_u64(f->in);
        blk.raw_len = read_u64(f->in);
        blk.records = read_u64(f->in);
        blk.min_allele = read_u64(f->in);
        blk.max_walk = read_u64(f->in);
    }
    return f;
}

// The catalog's file for `contig`: from a folder, <contig>.pgsites, or for a
// pangenome path name (CHM13#0#chr20) <chr20>.pgsites; null when there is none.
PgSitesFile* pgsites_for_contig(PgSitesCatalog& cat, const std::string& contig) {
    for (const auto& f : cat.files)
        if (chrom_names_match_vcf(contig, f->chrom)) return f.get();
    if (cat.dir.empty() || cat.absent.count(contig)) return nullptr;
    std::vector<std::string> names = {contig};
    if (chrom_suffix_match_key(contig) != contig) names.push_back(chrom_suffix_match_key(contig));
    std::string other;  // a file found for the contig on another sample's path
    for (const std::string& name : names) {
        const std::string path = cat.dir + "/" + name + kPgSitesExt;
        if (!is_pgsites_file(path)) continue;
        cat.files.push_back(open_pgsites(path));
        if (chrom_names_match_vcf(contig, cat.files.back()->chrom)) return cat.files.back().get();
        other = path + " holds " + cat.files.back()->chrom;
    }
    cat.absent.insert(contig);
    // Handles are per thread: say it once per run, not once per thread.
    static std::mutex mu;
    static std::unordered_set<std::string> warned;
    std::lock_guard<std::mutex> lock(mu);
    if (warned.insert(contig).second) {
        std::cerr << "warning: no sites on " << contig << ": ";
        if (other.empty()) std::cerr << "no " << chrom_suffix_match_key(contig) << kPgSitesExt << " in " << cat.dir << "\n";
        else std::cerr << other << ", another reference\n";
    }
    return nullptr;
}

// Lines of the blocks that may overlap [beg, end], kept when the record's
// [POS, POS + len(REF) - 1] overlaps it, as a tabix query would.
void append_graph_sites_pgsites(PgSitesCatalog& cat, const std::string& contig, hts_pos_t beg, hts_pos_t end,
                                hts_pos_t max_span, GraphSiteCatalog& catalog) {
    PgSitesFile* f = pgsites_for_contig(cat, contig);
    if (f == nullptr) return;
    std::string comp;
    for (const PgSitesFile::Block& blk : f->blocks) {
        if (static_cast<hts_pos_t>(blk.first_pos) > end || static_cast<hts_pos_t>(blk.max_end) < beg) continue;
        if (static_cast<hts_pos_t>(blk.min_allele) > max_span) {
            catalog.stats.beyond_read_span += blk.records;
            continue;
        }
        if (blk.records == 1 && blk.max_walk > kMaxSiteWalkSteps) {
            ++catalog.stats.walk_too_long;
            continue;
        }
        comp.resize(blk.comp_len);
        f->in.seekg(static_cast<std::streamoff>(blk.offset));
        if (!f->in.read(comp.data(), static_cast<std::streamsize>(blk.comp_len)))
            throw std::runtime_error("truncated .pgsites block");
        const std::string raw = zstd_decompress(comp, blk.raw_len);
        size_t i = 0;
        while (i < raw.size()) {
            size_t nl = raw.find('\n', i);
            if (nl == std::string::npos) nl = raw.size();
            const char* line = raw.data() + i;
            const size_t len = nl - i;
            // POS is field 2, REF field 4.
            const char* t1 = static_cast<const char*>(std::memchr(line, '\t', len));
            const char* t2 = t1 ? static_cast<const char*>(std::memchr(t1 + 1, '\t', len - (t1 + 1 - line))) : nullptr;
            const char* t3 = t2 ? static_cast<const char*>(std::memchr(t2 + 1, '\t', len - (t2 + 1 - line))) : nullptr;
            const char* t4 = t3 ? static_cast<const char*>(std::memchr(t3 + 1, '\t', len - (t3 + 1 - line))) : nullptr;
            if (t4 != nullptr) {
                const hts_pos_t pos = std::strtoll(t1 + 1, nullptr, 10);
                const hts_pos_t rec_end = pos + std::max<hts_pos_t>(1, t4 - (t3 + 1)) - 1;
                if (pos <= end && rec_end >= beg && !unmatchable_record(line, len, max_span, catalog.stats))
                    append_graph_site_from_vcf_data_line(line, len, catalog, false);
            }
            i = nl + 1;
        }
    }
}

bool has_tabix_index(const std::string& path) {
    tbx_t* tbx = tbx_index_load(path.c_str());
    if (!tbx) return false;
    tbx_destroy(tbx);
    return true;
}

}  // namespace


void require_site_catalog(const std::string& path) {
    if (is_directory(path)) {
        if (has_pgsites_file(path)) return;
        throw std::runtime_error("site catalog folder has no .pgsites files: " + path);
    }
    if (is_pgsites_file(path) || has_tabix_index(path)) return;
    throw std::runtime_error(
        "site catalog must be a .pgsites file, a folder of them, or a bgzipped VCF with a tabix "
        "index (.tbi or .csi): " + path);
}

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
        throw std::runtime_error("GBZ interval GFA did not contain a W line for contig: " + contig);
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
            site.level = 0;
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

// Build a "left+:right+" boundary query string for the GBZ FFI interval query.
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


void write_graph_site_catalog_tsv(std::ostream& out, const GraphSiteCatalog& catalog) {
    out << "CHROM\tREF_CONTIG\tPOS\tEND\tID\tLEVEL\tPARENT\tROOT\tCOND_PARENT_ALLELES\tN_ALTS\tN_TRAVERSALS\tHAS_STAR\tELIGIBLE\tSKIP_REASON\tREF\tALT\n";
    for (const GraphSite& site : catalog.sites) {
        const std::string skip_reason = graph_site_validation_skip_reason(site);
        out << site.chrom << '\t'
            << (site.ref_contig.empty() ? site.chrom : site.ref_contig) << '\t'
            << site.pos << '\t'
            << site.ref_end << '\t'
            << site.id << '\t'
            << site.level << '\t'
            << (site.parent.empty() ? "." : site.parent) << '\t'
            << (site.root.empty() ? "." : site.root) << '\t';
        for (size_t i = 0; i < site.conditional_parent_alleles.size(); ++i) {
            if (i > 0) out << ',';
            out << site.conditional_parent_alleles[i];
        }
        if (site.conditional_parent_alleles.empty()) out << '.';
        out << '\t'
            << site.alts.size() << '\t'
            << site.allele_traversals.size() << '\t'
            << (site.has_spanning_deletion ? 1 : 0) << '\t'
            << (skip_reason.empty() && site.eligible ? 1 : 0) << '\t'
            << (skip_reason.empty() ? "." : skip_reason) << '\t'
            << site.ref << '\t';
        for (size_t i = 0; i < site.alts.size(); ++i) {
            if (i > 0) out << ',';
            out << site.alts[i];
        }
        out << '\n';
    }
}

// --- SiteCatalogHandle RAII implementation ---

SiteCatalogHandle::SiteCatalogHandle(const std::string& path) {
    if (is_directory(path)) {
        pgsites = std::make_unique<PgSitesCatalog>();
        pgsites->dir = path;
        return;
    }
    if (is_pgsites_file(path)) {
        pgsites = std::make_unique<PgSitesCatalog>();
        pgsites->files.push_back(open_pgsites(path));
        return;
    }
    fp = hts_open(path.c_str(), "r");
    if (fp == nullptr)
        throw std::runtime_error("failed to open sites VCF: " + path);
    tbx = tbx_index_load(path.c_str());
    if (tbx == nullptr) {
        hts_close(fp);
        fp = nullptr;
        throw std::runtime_error("sites VCF requires a tabix index (.tbi): " + path);
    }
}

SiteCatalogHandle::~SiteCatalogHandle() {
    if (tbx != nullptr) tbx_destroy(tbx);
    if (fp != nullptr) hts_close(fp);
}

SiteCatalogHandle::SiteCatalogHandle(SiteCatalogHandle&& o) noexcept
    : fp(o.fp), tbx(o.tbx), pgsites(std::move(o.pgsites)), max_span(std::move(o.max_span)) {
    o.fp  = nullptr;
    o.tbx = nullptr;
}

SiteCatalogHandle& SiteCatalogHandle::operator=(SiteCatalogHandle&& o) noexcept {
    if (this != &o) {
        if (tbx != nullptr) tbx_destroy(tbx);
        if (fp != nullptr) hts_close(fp);
        fp    = o.fp;
        tbx   = o.tbx;
        pgsites = std::move(o.pgsites);
        max_span = std::move(o.max_span);
        o.fp  = nullptr;
        o.tbx = nullptr;
    }
    return *this;
}

GraphSiteCatalog load_sites_for_region(SiteCatalogHandle& handle,
                                       const std::string& contig,
                                       hts_pos_t beg, hts_pos_t end) {
    GraphSiteCatalog catalog;
    if (!handle || contig.empty() || end <= beg) return catalog;

    hts_pos_t max_span = std::numeric_limits<hts_pos_t>::max();
    const auto span = handle.max_span.find(chrom_suffix_match_key(contig));
    if (span != handle.max_span.end()) max_span = span->second;
    if (handle.pgsites) {
        append_graph_sites_pgsites(*handle.pgsites, contig, beg, end, max_span, catalog);
    } else {
        std::vector<RegionFilter> filters = {
            {true, contig, beg, end}
        };
        append_graph_sites_tabix_filtered(handle.fp, handle.tbx, filters, max_span, catalog, false);
    }
    finalize_graph_site_catalog_inplace(catalog);
    return catalog;
}

} // namespace pgphase_collect
