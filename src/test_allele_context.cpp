#include "allele_context.hpp"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <random>
#include <set>
#include <stdexcept>

using namespace pgphase_collect;

static int failures = 0;
static size_t checks = 0;
static void check(bool ok, const char* message) {
    ++checks;
    if (!ok) { ++failures; std::cerr << "FAIL: " << message << '\n'; }
}

static std::unique_ptr<bam1_t, decltype(&bam_destroy1)> read(
        hts_pos_t pos, const std::string& cigar, const std::string& seq, int flag = 0,
        const std::string& name = "read") {
    std::unique_ptr<sam_hdr_t, decltype(&sam_hdr_destroy)> header(
        sam_hdr_parse(31, "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n"), sam_hdr_destroy);
    std::unique_ptr<bam1_t, decltype(&bam_destroy1)> bam(bam_init1(), bam_destroy1);
    const std::string line = name + "\t" + std::to_string(flag) + "\tchr1\t" + std::to_string(pos) +
        "\t60\t" + cigar + "\t*\t0\t0\t" + seq + "\t" + std::string(seq.size(), 'I');
    kstring_t input{0, 0, nullptr};
    kputs(line.c_str(), &input);
    check(sam_parse1(&input, header.get(), bam.get()) >= 0, "parse original alignment");
    std::free(input.s);
    return bam;
}

// Independent scalar dynamic programming, without edlib or production scoring.
static int distance(const std::string& a, const std::string& b) {
    std::vector<int> row(b.size() + 1);
    std::iota(row.begin(), row.end(), 0);
    for (size_t i = 0; i < a.size(); ++i) {
        int diagonal = row[0];
        row[0] = i + 1;
        for (size_t j = 0; j < b.size(); ++j) {
            const int old = row[j + 1];
            row[j + 1] = std::min({row[j] + 1, row[j + 1] + 1, diagonal + (a[i] != b[j])});
            diagonal = old;
        }
    }
    return row.back();
}

// Enumerate every optimal traceback and intersect its exact matches. This
// oracle does not use production forward/backward edge certification.
static std::vector<int> reference_map_oracle(const std::string& ref, const std::string& alt) {
    std::vector<int> path(ref.size(), -1), common;
    bool first = true;
    const int best = distance(ref, alt);
    std::function<void(size_t, size_t, int)> visit = [&](size_t i, size_t j, int cost) {
        if (cost + distance(ref.substr(i), alt.substr(j)) != best) return;
        if (i == ref.size() && j == alt.size()) {
            if (first) { common = path; first = false; }
            else for (size_t k = 0; k < path.size(); ++k)
                if (common[k] != path[k]) common[k] = -1;
            return;
        }
        if (i < ref.size() && j < alt.size()) {
            path[i] = ref[i] == alt[j] ? static_cast<int>(j) : -1;
            visit(i + 1, j + 1, cost + (ref[i] != alt[j]));
        }
        if (i < ref.size()) { path[i] = -1; visit(i + 1, j, cost + 1); }
        if (j < alt.size()) visit(i, j + 1, cost + 1);
    };
    visit(0, 0, 0);
    return common;
}

static void check_reference_map(const std::string& ref, const std::string& alt) {
    const auto mapping = map_allele_reference(ref, alt);
    check(mapping.status == AlleleMapStatus::Complete, "all-optimal reference map completes");
    check(mapping.distance == distance(ref, alt), "reference map has independent optimal distance");
    std::vector<int> actual(ref.size(), -1);
    size_t previous_ref_end = 0, previous_alt_end = 0;
    for (const auto& span : mapping.matches) {
        check(span.length > 0 && span.ref_beg + span.length <= ref.size() &&
            span.alt_beg + span.length <= alt.size(), "matched subpath has valid bounds");
        check(span.ref_beg >= previous_ref_end && span.alt_beg >= previous_alt_end,
            "matched subpaths remain ordered and disjoint");
        check(span.ref_beg != previous_ref_end || span.alt_beg != previous_alt_end || previous_ref_end == 0,
            "consecutive matched subpaths are coalesced");
        check(ref.substr(span.ref_beg, span.length) == alt.substr(span.alt_beg, span.length),
            "mapped bases are exact matches");
        for (size_t k = 0; k < span.length; ++k) actual[span.ref_beg + k] = span.alt_beg + k;
        previous_ref_end = span.ref_beg + span.length;
        previous_alt_end = span.alt_beg + span.length;
    }
    check(actual == reference_map_oracle(ref, alt), "matched map equals intersection of every optimal traceback");
}

static std::vector<std::string> split(const std::string& line, char separator) {
    std::vector<std::string> fields;
    size_t beg = 0;
    for (;;) {
        const size_t end = line.find(separator, beg);
        fields.push_back(line.substr(beg, end - beg));
        if (end == std::string::npos) return fields;
        beg = end + 1;
    }
}

static int replay(const char* context_path, const char* sequence_path) {
    std::ifstream contexts_file(context_path);
    std::ifstream sequences_file(sequence_path);
    check(contexts_file.good() && sequences_file.good(), "saved sequence evidence is readable");
    std::map<std::string, AlleleSequenceContext> contexts;
    std::string line;
    std::getline(contexts_file, line);
    while (std::getline(contexts_file, line)) {
        const auto f = split(line, '\t');
        check(f.size() == 8, "saved context field count");
        if (f.size() != 8 || f[1] != "valid") continue;
        contexts.emplace(f[0], AlleleSequenceContext{std::stoll(f[2]), std::stoll(f[3]),
            {f[4], f[5]}, f[6].empty() ? std::vector<std::string>{} : split(f[6], ',')});
    }
    size_t scored = 0;
    std::getline(sequences_file, line);
    while (std::getline(sequences_file, line)) {
        const auto f = split(line, '\t');
        check(f.size() == 13, "saved molecule field count");
        if (f.size() != 13 || f[3] != "scored") continue;
        const auto context = contexts.find(f[0]);
        check(context != contexts.end(), "saved query has complete context");
        if (context == contexts.end()) continue;
        const auto score = score_allele_context(context->second, f[7]);
        check(std::stoi(f[5]) - std::stoi(f[4]) == static_cast<int>(f[7].size()) &&
            f[8].size() == f[7].size() * 2, "saved physical query bounds and qualities");
        check(score.distances[0] == std::stoi(f[9]) && score.distances[1] == std::stoi(f[10]) &&
            score.other_distance == std::stoi(f[11]) && score.nearest == std::stoi(f[12]),
            "saved-state production scoring replay");
        ++scored;
    }
    check(scored > 0, "saved replay cannot silently skip every molecule");
    std::cout << scored << " scored molecules, " << checks << " checks, " << failures << " failures\n";
    return failures != 0;
}

static int replay_compositions(const std::string& folder, const std::string& prefix, bool multi_neighbor = false,
        bool reference_maps = false, bool nested = false, bool read_catalogs = false) {
    const auto table = [&](const std::string& suffix) {
        std::ifstream input(folder + "/matrix.chunk0." + suffix + ".tsv");
        if (!input) throw std::runtime_error("cannot read composition state: " + folder + "/" + suffix);
        std::vector<std::vector<std::string>> rows;
        std::string line;
        std::getline(input, line);
        while (std::getline(input, line)) rows.push_back(split(line, '\t'));
        return rows;
    };
    std::map<size_t, std::vector<std::string>> alleles;
    for (const auto& row : table("physical-alleles")) alleles[std::stoull(row[0])].push_back(row[2]);
    std::map<size_t, std::set<size_t>> candidates_by_locus;
    for (const auto& row : table("joint-parents")) candidates_by_locus[std::stoull(row[0])].insert(std::stoull(row[1]));
    std::map<size_t, std::set<size_t>> parents;
    for (const auto& row : table("physical-members")) {
        const auto& candidates = candidates_by_locus.at(std::stoull(row[1]));
        parents[std::stoull(row[0])].insert(candidates.begin(), candidates.end());
    }
    const auto catalog = table("site-catalog");
    std::map<size_t, std::vector<MoleculeAlleleSequence>> observations;
    if (read_catalogs) for (const auto& row : table("physical-cohort"))
        observations[std::stoull(row[0])].push_back({row[1], row[5]});
    std::ofstream sites(prefix + ".composition-sites.tsv"), compositions(prefix + ".compositions.tsv"), hypotheses(prefix + ".composed-alleles.tsv");
    if (!sites || !compositions || !hypotheses) throw std::runtime_error("cannot write composition state: " + prefix);
    sites << "physical\tcandidate\tpos\tref\talt_index\talt\tbam_injected\tmsa_verified\tcategory\n";
    compositions << "physical\tparent_allele\tcandidate\talt_index\tstatus\tsequence\n";
    hypotheses << "physical\tallele\tsequence\n";
    std::ofstream path_status, sequence_paths, path_alleles;
    if (multi_neighbor) {
        path_status.open(prefix + ".path-status.tsv");
        sequence_paths.open(prefix + ".sequence-paths.tsv");
        path_alleles.open(prefix + ".path-alleles.tsv");
        if (!path_status || !sequence_paths || !path_alleles)
            throw std::runtime_error("cannot write path state: " + prefix);
        path_status << "physical\tstatus\tvisited_prefixes\toverlap_prefixes\tunsupported_paths\tpaths\n";
        sequence_paths << "physical\tpath\tparent_allele\tneighbors\tsequence\n";
        path_alleles << "physical\tallele\tsequence\n";
    }
    size_t count = 0;
    std::ofstream parent_maps, matched_subpaths;
    std::ofstream nested_compositions, nested_alleles;
    std::ofstream read_status, read_hypotheses, read_support, read_exclusions, read_alleles;
    if (read_catalogs) {
        read_status.open(prefix + ".read-catalog-status.tsv");
        read_hypotheses.open(prefix + ".read-hypotheses.tsv");
        read_support.open(prefix + ".read-support.tsv");
        read_exclusions.open(prefix + ".read-exclusions.tsv");
        read_alleles.open(prefix + ".read-alleles.tsv");
        if (!read_status || !read_hypotheses || !read_support || !read_exclusions || !read_alleles)
            throw std::runtime_error("cannot write read catalog: " + prefix);
        read_status << "physical\tstatus\tmolecules\tconflicting\tunsupported\thypotheses\n";
        read_hypotheses << "physical\thypothesis\tsequence\tmolecules\n";
        read_support << "physical\thypothesis\tmolecule\n";
        read_exclusions << "physical\tmolecule\treason\n";
        read_alleles << "physical\tallele\tsequence\n";
    }
    if (nested) {
        nested_compositions.open(prefix + ".nested-compositions.tsv");
        nested_alleles.open(prefix + ".nested-alleles.tsv");
        if (!nested_compositions || !nested_alleles) throw std::runtime_error("cannot write nested composition state: " + prefix);
        nested_compositions << "physical\tparent_allele\tcandidate\talt_index\tstatus\tsequence\n";
        nested_alleles << "physical\tallele\tsequence\n";
    }
    if (reference_maps) {
        parent_maps.open(prefix + ".parent-maps.tsv");
        matched_subpaths.open(prefix + ".matched-subpaths.tsv");
        if (!parent_maps || !matched_subpaths) throw std::runtime_error("cannot write reference maps: " + prefix);
        parent_maps << "physical\tparent_allele\tstatus\tdistance\tmatched_bases\tsubpaths\n";
        matched_subpaths << "physical\tparent_allele\tref_beg\talt_beg\tlength\n";
    }
    for (const auto& context : table("composition-contexts")) {
        const size_t id = std::stoull(context[0]);
        const hts_pos_t beg = std::stoll(context[1]);
        const hts_pos_t end = std::stoll(context[2]);
        std::vector<AlleleReferenceMap> parent_mappings;
        if (reference_maps) for (size_t parent = 0; parent < alleles.at(id).size(); ++parent) {
            parent_mappings.push_back(map_allele_reference(context[3], alleles.at(id)[parent]));
            const auto& mapping = parent_mappings.back();
            const char* status = mapping.status == AlleleMapStatus::Complete ? "complete"
                : mapping.status == AlleleMapStatus::Limited ? "limited" : "unsupported";
            size_t matched = 0;
            for (const auto& span : mapping.matches) {
                matched += span.length;
                matched_subpaths << id << '\t' << parent << '\t' << span.ref_beg << '\t'
                    << span.alt_beg << '\t' << span.length << '\n';
            }
            parent_maps << id << '\t' << parent << '\t' << status << '\t' << mapping.distance
                << '\t' << matched << '\t' << mapping.matches.size() << '\n';
        }
        const auto base = [&](hts_pos_t pos) { return pos >= beg && pos <= end ? context[3][pos - beg] : 'N'; };
        std::set<std::string> complete(alleles.at(id).begin(), alleles.at(id).end());
        std::set<std::string> nested_sequences;
        std::map<size_t, AlleleNeighborSite> scoped;
        for (const auto& row : catalog) {
            const size_t candidate = std::stoull(row[0]);
            const hts_pos_t pos = std::stoll(row[1]);
            if (parents[id].count(candidate) || row[2].empty() || pos <= beg ||
                pos + static_cast<hts_pos_t>(row[2].size()) - 1 >= end) continue;
            sites << id;
            for (const auto& field : row) sites << '\t' << field;
            sites << '\n';
            auto [site, inserted] = scoped.emplace(candidate, AlleleNeighborSite{candidate, pos, row[2], {}});
            (void)inserted;
            site->second.alts.push_back(row[4]);
            for (size_t parent = 0; parent < alleles[id].size(); ++parent) {
                const auto result = compose_allele_sequence(beg, end,
                    {{beg, context[3], alleles[id][parent]}, {pos, row[2], row[4]}}, base);
                const char* status = result.status == AlleleCompositionStatus::Valid ? "valid"
                    : result.status == AlleleCompositionStatus::Overlap ? "overlap" : "unsupported";
                compositions << id << '\t' << parent << '\t' << candidate << '\t' << row[3] << '\t'
                    << status << '\t' << (result.sequence.empty() ? "." : result.sequence) << '\n';
                if (result.status == AlleleCompositionStatus::Valid) complete.insert(result.sequence);
                if (nested) {
                    const auto result = compose_allele_on_subpaths(beg, parent_mappings[parent], {{pos, row[2], row[4]}});
                    const char* status = result.status == AlleleCompositionStatus::Valid ? "valid"
                        : result.status == AlleleCompositionStatus::Overlap ? "overlap" : "unsupported";
                    nested_compositions << id << '\t' << parent << '\t' << candidate << '\t' << row[3] << '\t'
                        << status << '\t' << (result.sequence.empty() ? "." : result.sequence) << '\n';
                    if (result.status == AlleleCompositionStatus::Valid) nested_sequences.insert(result.sequence);
                }
                ++count;
            }
        }
        size_t allele = 0;
        for (const auto& sequence : complete) hypotheses << id << '\t' << allele++ << '\t' << sequence << '\n';
        if (multi_neighbor) {
            std::vector<AlleleNeighborSite> neighbors;
            for (const auto& [candidate, site] : scoped) { (void)candidate; neighbors.push_back(site); }
            const auto catalog = build_allele_sequence_paths(beg, end, alleles.at(id), neighbors, base);
            const char* status = catalog.status == AllelePathStatus::Complete ? "complete"
                : catalog.status == AllelePathStatus::Limited ? "limited" : "unsupported";
            path_status << id << '\t' << status << '\t' << catalog.visited_prefixes << '\t'
                << catalog.overlap_prefixes << '\t' << catalog.unsupported_paths << '\t' << catalog.paths.size() << '\n';
            for (size_t path = 0; path < catalog.paths.size(); ++path) {
                const auto& entry = catalog.paths[path];
                sequence_paths << id << '\t' << path << '\t' << entry.parent_allele << '\t';
                if (entry.neighbors.empty()) sequence_paths << '.';
                for (size_t i = 0; i < entry.neighbors.size(); ++i) {
                    if (i) sequence_paths << ',';
                    sequence_paths << entry.neighbors[i].first << ':' << entry.neighbors[i].second;
                }
                sequence_paths << '\t' << entry.sequence << '\n';
                complete.insert(entry.sequence);
            }
            allele = 0;
            for (const auto& sequence : complete) path_alleles << id << '\t' << allele++ << '\t' << sequence << '\n';
        }
        if (nested) {
            complete.insert(nested_sequences.begin(), nested_sequences.end());
            allele = 0;
            for (const auto& sequence : complete) nested_alleles << id << '\t' << allele++ << '\t' << sequence << '\n';
        }
        if (read_catalogs) {
            const auto catalog = build_read_allele_catalog(observations[id]);
            read_status << id << '\t' << (catalog.status == ReadAlleleCatalogStatus::Complete ? "complete" : "limited")
                << '\t' << catalog.eligible_molecules << '\t' << catalog.conflicting_molecules.size()
                << '\t' << catalog.unsupported_molecules.size() << '\t' << catalog.hypotheses.size() << '\n';
            for (const auto& name : catalog.conflicting_molecules) read_exclusions << id << '\t' << name << "\tconflicting\n";
            for (const auto& name : catalog.unsupported_molecules) read_exclusions << id << '\t' << name << "\tunsupported\n";
            for (size_t i = 0; i < catalog.hypotheses.size(); ++i) {
                const auto& hypothesis = catalog.hypotheses[i];
                read_hypotheses << id << '\t' << i << '\t' << hypothesis.sequence << '\t' << hypothesis.molecules.size() << '\n';
                for (const auto& name : hypothesis.molecules) read_support << id << '\t' << i << '\t' << name << '\n';
                complete.insert(hypothesis.sequence);
            }
            allele = 0;
            for (const auto& sequence : complete) read_alleles << id << '\t' << allele++ << '\t' << sequence << '\n';
        }
    }
    std::cout << count << " parent/neighbor compositions\n";
    return 0;
}

int main(int argc, char** argv) {
    if (argc == 4 && std::string(argv[1]) == "--compose") return replay_compositions(argv[2], argv[3]);
    if (argc == 4 && std::string(argv[1]) == "--paths") return replay_compositions(argv[2], argv[3], true);
    if (argc == 4 && std::string(argv[1]) == "--maps") return replay_compositions(argv[2], argv[3], true, true);
    if (argc == 4 && std::string(argv[1]) == "--nested") return replay_compositions(argv[2], argv[3], true, true, true);
    if (argc == 4 && std::string(argv[1]) == "--read-catalogs") return replay_compositions(argv[2], argv[3], true, true, true, true);
    if (argc == 4 && std::string(argv[1]) == "--replay") return replay(argv[2], argv[3]);
    if (argc != 1) {
        std::cerr << "Usage: test_allele_context [--replay CONTEXTS.tsv SEQUENCES.tsv | --compose INPUT_FOLDER OUTPUT_PREFIX | --paths INPUT_FOLDER OUTPUT_PREFIX | --maps INPUT_FOLDER OUTPUT_PREFIX | --nested INPUT_FOLDER OUTPUT_PREFIX | --read-catalogs INPUT_FOLDER OUTPUT_PREFIX]\n";
        return 1;
    }
    std::mt19937 random(416);
    const std::vector<MoleculeAlleleSequence> discovery{
        {"a", "ACGT"}, {"a", "acgt"}, {"b", "ACGT"}, {"c", "TGCA"}, {"d", "TGCA"}, {"e", "TGCA"},
        {"f", "AAAA"}, {"x", "CCCC"}, {"x", "GGGG"}, {"y", "N"}, {"y", "ACGT"}};
    const auto discovered = build_read_allele_catalog(discovery);
    check(discovered.status == ReadAlleleCatalogStatus::Complete && discovered.eligible_molecules == 6 &&
        discovered.hypotheses.size() == 2, "read catalog counts distinct eligible molecules without counting aliases");
    check(discovered.hypotheses[0].sequence == "ACGT" && discovered.hypotheses[0].molecules == std::vector<std::string>{"a", "b"},
        "two exact physical supporters establish a sequence hypothesis");
    check(discovered.hypotheses[1].sequence == "TGCA" && discovered.hypotheses[1].molecules == std::vector<std::string>{"c", "d", "e"},
        "hypotheses and supporters have canonical order");
    check(discovered.conflicting_molecules == std::vector<std::string>{"x"} &&
        discovered.unsupported_molecules == std::vector<std::string>{"y"}, "conflicting or invalid names abstain completely");
    check(build_read_allele_catalog({{"a", "ACGT"}, {"a", "ACGT"}}).hypotheses.empty(),
        "duplicated descriptions of one molecule cannot discover an allele");
    const auto heldout_two = build_read_allele_catalog(discovery, "a");
    check(heldout_two.eligible_molecules == 5 && heldout_two.hypotheses.size() == 1 &&
        heldout_two.hypotheses[0].sequence == "TGCA", "excluding a molecule removes every alias before discovery");
    const auto heldout_three = build_read_allele_catalog(discovery, "c");
    check(heldout_three.hypotheses.size() == 2 && heldout_three.hypotheses[1].molecules == std::vector<std::string>{"d", "e"},
        "a three-molecule hypothesis retains independent support after exclusion");
    check(build_read_allele_catalog(discovery, "x").conflicting_molecules.empty(), "excluded conflicts cannot poison discovery");
    check(build_read_allele_catalog(discovery, "y").unsupported_molecules.empty(), "excluded invalid names cannot poison discovery");
    for (const size_t cap : {size_t(0), size_t(1)}) {
        const auto limited = build_read_allele_catalog(discovery, {}, cap);
        check(limited.status == ReadAlleleCatalogStatus::Limited && limited.hypotheses.empty() &&
            limited.eligible_molecules == 6, "limited read catalogs clear every partial hypothesis but retain admission counts");
    }
    check(build_read_allele_catalog(discovery, {}, 2).status == ReadAlleleCatalogStatus::Complete,
        "exact read hypothesis budget completes");
    check(build_read_allele_catalog({}, {}, 0).status == ReadAlleleCatalogStatus::Complete,
        "empty read discovery is a complete empty catalog");
    for (const std::string& invalid : {std::string(), std::string("N"), std::string("*"), std::string(4097, 'A')}) {
        const auto rejected = build_read_allele_catalog({{"a", invalid}, {"b", invalid}});
        check(rejected.eligible_molecules == 0 && rejected.hypotheses.empty() && rejected.unsupported_molecules.size() == 2,
            "unsupported read DNA and length cannot discover sequences");
    }
    check(build_read_allele_catalog({{"", "ACGT"}, {"b", "ACGT"}}).hypotheses.empty(), "unnamed reads cannot contribute support");
    const auto long_reads = build_read_allele_catalog({{"a", std::string(312, 'A')}, {"b", std::string(312, 'A')}, {"c", std::string(55, 'A')}});
    check(long_reads.hypotheses.size() == 1 && long_reads.hypotheses[0].sequence.size() == 312,
        "read hypotheses retain complete insertion sequence beyond the existing reference length");
    const auto maximum_reads = build_read_allele_catalog({{"a", std::string(4096, 'A')}, {"b", std::string(4096, 'A')}});
    check(maximum_reads.hypotheses.size() == 1, "maximum supported context sequence length remains admissible");
    for (size_t trial = 0; trial < 200; ++trial) {
        auto permuted = discovery;
        std::shuffle(permuted.begin(), permuted.end(), random);
        const auto reordered = build_read_allele_catalog(permuted);
        check(reordered.eligible_molecules == discovered.eligible_molecules && reordered.conflicting_molecules == discovered.conflicting_molecules &&
            reordered.unsupported_molecules == discovered.unsupported_molecules && reordered.hypotheses.size() == 2 &&
            reordered.hypotheses[0].sequence == discovered.hypotheses[0].sequence && reordered.hypotheses[0].molecules == discovered.hypotheses[0].molecules &&
            reordered.hypotheses[1].sequence == discovered.hypotheses[1].sequence && reordered.hypotheses[1].molecules == discovered.hypotheses[1].molecules,
            "read source order and conflict order do not alter physical discovery");
    }
    random.seed(416);
    const auto compound_parent = map_allele_reference("ACGTA", "TCGTT");
    const auto nested_check = [&](const std::vector<AlleleSequenceEdit>& edits, const std::string& expected) {
        const auto result = compose_allele_on_subpaths(100, compound_parent, edits);
        check(result.status == AlleleCompositionStatus::Valid && result.sequence == expected,
            "nested edit retains parent changes and maps into its matched island");
    };
    nested_check({{102, "G", "A"}}, "TCATT");
    nested_check({{101, "C", "A"}}, "TAGTT");
    nested_check({{103, "T", "A"}}, "TCGAT");
    nested_check({{102, "G", "GA"}}, "TCGATT");
    nested_check({{102, "G", "GA"}, {103, "T", "C"}}, "TCGACT");
    nested_check({{103, "T", "C"}, {102, "G", "GA"}}, "TCGACT");
    nested_check({{102, "G", ""}}, "TCTT");
    nested_check({{102, "G", "A"}, {102, "G", "A"}}, "TCATT");
    nested_check({{101, "CG", "CA"}, {102, "G", "A"}}, "TCATT");
    nested_check({{101, "C", "A"}, {103, "T", "A"}}, "TAGAT");
    nested_check({{103, "T", "A"}, {101, "C", "A"}}, "TAGAT");
    nested_check({{100, "A", "A"}, {102, "G", "A"}}, "TCATT");
    nested_check({{102, "g", "a"}}, "TCATT");
    for (const auto& edits : std::vector<std::vector<AlleleSequenceEdit>>{
            {{100, "A", "C"}}, {{101, "C", "CA"}}, {{103, "T", "TA"}},
            {{102, "G", "A"}, {102, "G", "C"}}, {{101, "CG", ""}, {102, "G", "A"}}}) {
        const auto result = compose_allele_on_subpaths(100, compound_parent, edits);
        check(result.status == AlleleCompositionStatus::Overlap && result.sequence.empty(),
            "unmapped bases, competing edits and unguarded indels remain unresolved");
    }
    for (const auto& edits : std::vector<std::vector<AlleleSequenceEdit>>{
            {{102, "A", "C"}}, {{99, "A", "C"}}, {{105, "A", "C"}}, {{102, "G", "N"}}}) {
        const auto result = compose_allele_on_subpaths(100, compound_parent, edits);
        check(result.status == AlleleCompositionStatus::Unsupported && result.sequence.empty(),
            "nested composition validates raw REF, bounds and DNA");
    }
    const auto limited_parent = map_allele_reference("ACGTA", "TCGTT", 0);
    check(compose_allele_on_subpaths(100, limited_parent, {{102, "G", "A"}}).status ==
        AlleleCompositionStatus::Unsupported, "limited mapping cannot authorize nested composition");
    check(compose_allele_on_subpaths(100, limited_parent, {}).sequence == "TCGTT",
        "mapping limit does not discard an already supported anchored parent");
    check(compose_allele_on_subpaths(0, compound_parent, {}).status == AlleleCompositionStatus::Unsupported,
        "nested contexts require one-based coordinates");
    check(compose_allele_on_subpaths(1, map_allele_reference("", ""), {}).status ==
        AlleleCompositionStatus::Unsupported, "empty mapping is not a physical reference window");
    const auto boundary_parent = map_allele_reference("AC", "TAC");
    const auto shifted_parent = map_allele_reference("ACGTA", "TTTCGTT");
    const auto shifted_nested = compose_allele_on_subpaths(100, shifted_parent, {{102, "G", "A"}});
    check(shifted_nested.status == AlleleCompositionStatus::Valid && shifted_nested.sequence == "TTTCATT",
        "nested edits use parent ALT offsets after a length-changing parent edit");
    check(compose_allele_on_subpaths(1, boundary_parent, {{1, "AC", "GAC"}}).status ==
        AlleleCompositionStatus::Overlap, "matched bases do not order competing boundary insertions");
    const auto repeat_parent = map_allele_reference("CAAAG", "CAAAAG");
    check(compose_allele_on_subpaths(1, repeat_parent, {{3, "A", "C"}}).status ==
        AlleleCompositionStatus::Overlap, "ambiguous repeat mapping cannot place a neighbor SNP");
    for (size_t trial = 0; trial < 500; ++trial) {
        std::string ref(8, 'A');
        for (char& b : ref) b = "ACGT"[random() % 4];
        std::string parent = ref;
        parent.front() = ref.front() == 'A' ? 'C' : 'A';
        parent.back() = ref.back() == 'A' ? 'C' : 'A';
        const auto oracle = reference_map_oracle(ref, parent);
        const auto mapping = map_allele_reference(ref, parent);
        const size_t offset = 1 + random() % 6;
        const std::string alt(1, ref[offset] == 'A' ? 'C' : 'A');
        const auto result = compose_allele_on_subpaths(100, mapping,
            {{static_cast<hts_pos_t>(100 + offset), ref.substr(offset, 1), alt}});
        if (oracle[offset] < 0)
            check(result.status == AlleleCompositionStatus::Overlap && result.sequence.empty(),
                "random nested substitutions abstain on ambiguous or altered parent bases");
        else {
            parent[oracle[offset]] = alt[0];
            check(result.status == AlleleCompositionStatus::Valid && result.sequence == parent,
                "random nested substitutions match independent all-optimal mapped coordinates");
        }
    }
    random.seed(416);
    for (const auto& pair : std::vector<std::pair<std::string, std::string>>{
            {"ACGTA", "TCGTT"}, {"CAAAG", "CAAAAG"}, {"CAAAAG", "CAAAG"},
            {"CAG", "CTAG"}, {"CTAG", "CAG"}, {"AG", "GA"}, {"", ""}, {"", "AC"}, {"AC", ""}})
        check_reference_map(pair.first, pair.second);
    const auto internal = map_allele_reference("ACGTA", "TCGTT");
    check(internal.matches.size() == 1 && internal.matches[0].ref_beg == 1 &&
        internal.matches[0].alt_beg == 1 && internal.matches[0].length == 3,
        "compound parent exposes its verified internal matched island");
    const auto repeated = map_allele_reference("CAAAG", "CAAAAG");
    check(repeated.matches.size() == 2 && repeated.matches[0].length == 1 &&
        repeated.matches[1].ref_beg == 4 && repeated.matches[1].alt_beg == 5,
        "repeat placements abstain while unique outer anchors survive");
    check(map_allele_reference("AG", "GA").matches.empty(), "equally optimal substitutions and gaps stay unresolved");
    const auto lowercase = map_allele_reference("acgt", "ACGT");
    check(lowercase.status == AlleleMapStatus::Complete && lowercase.matches.size() == 1 &&
        lowercase.matches[0].length == 4, "mapping normalizes DNA case");
    for (const auto& pair : std::vector<std::pair<std::string, std::string>>{
            {"AN", "AA"}, {"AA", "N"}, {"*", "A"}, {std::string(4097, 'A'), "A"}}) {
        const auto invalid = map_allele_reference(pair.first, pair.second);
        check(invalid.status == AlleleMapStatus::Unsupported && invalid.distance == -1 && invalid.matches.empty(),
            "unsupported mapping never publishes an alignment");
    }
    check(map_allele_reference("AC", "AC", 9).status == AlleleMapStatus::Complete,
        "exact cell budget completes");
    for (const size_t budget : {size_t(0), size_t(8)}) {
        const auto limited = map_allele_reference("AC", "AC", budget);
        check(limited.status == AlleleMapStatus::Limited && limited.distance == -1 && limited.matches.empty(),
            "insufficient work budget leaves no partial matches or cost");
    }
    check(map_allele_reference(std::string(4096, 'A'), std::string(4096, 'A')).status == AlleleMapStatus::Limited,
        "maximum diagnostic DNA length remains bounded by cell budget");
    for (size_t trial = 0; trial < 500; ++trial) {
        std::string ref(random() % 7, 'A'), alt(random() % 7, 'A');
        for (char& b : ref) b = "ACGT"[random() % 4];
        for (char& b : alt) b = "ACGT"[random() % 4];
        check_reference_map(ref, alt);
    }
    // Keep all pre-existing random fixtures on their original seed/state.
    random.seed(416);
    std::string reference(500, 'A');
    for (char& base : reference) base = "ACGT"[random() % 4];
    const auto base = [&](hts_pos_t pos) { return pos > 0 && pos <= 500 ? reference[pos - 1] : 'N'; };
    const auto substitute = [](char b) { return std::string(1, b == 'A' ? 'C' : 'A'); };
    const std::vector<AlleleSequenceEdit> neighboring{
        {80, reference.substr(79, 1), substitute(base(80))},
        {120, reference.substr(119, 1), reference.substr(119, 1) + "GG"},
        {160, reference.substr(159, 3), reference.substr(159, 1)}};
    std::string expected_composed = reference.substr(49, 201);
    for (auto i = neighboring.rbegin(); i != neighboring.rend(); ++i)
        expected_composed.replace(i->pos - 50, i->ref.size(), i->alt);
    const auto composed = compose_allele_sequence(50, 250, neighboring, base);
    check(composed.status == AlleleCompositionStatus::Valid && composed.sequence == expected_composed,
        "compose SNP, insertion and deletion in original reference coordinates");
    check(compose_allele_sequence(50, 250, {}, base).sequence == reference.substr(49, 201),
        "empty edit selection retains complete reference sequence");
    const auto conflict = compose_allele_sequence(50, 250,
        {{80, reference.substr(79, 1), substitute(base(80))},
         {80, reference.substr(79, 1), base(80) == 'G' ? "T" : "G"}}, base);
    check(conflict.status == AlleleCompositionStatus::Overlap && conflict.sequence.empty(),
        "contradictory substitutions do not overwrite one another");
    check(compose_allele_sequence(50, 250,
        {{80, reference.substr(79, 5), ""}, {82, reference.substr(81, 1), substitute(base(82))}}, base).status ==
        AlleleCompositionStatus::Overlap, "deletion masking another edit remains unresolved");
    check(compose_allele_sequence(50, 250,
        {{80, reference.substr(79, 5), ""}, {82, "", "CG"}}, base).status ==
        AlleleCompositionStatus::Overlap, "insertion inside a consumed interval remains unresolved");
    check(compose_allele_sequence(50, 250,
        {{80, "", "C"}, {80, "", "G"}}, base).status == AlleleCompositionStatus::Overlap,
        "different insertions at one boundary cannot acquire an invented order");
    const auto boundary = compose_allele_sequence(50, 250,
        {{80, "", std::string(1, base(79) == 'C' ? 'G' : 'C')},
         {80, reference.substr(79, 1), substitute(base(80))}}, base);
    std::string boundary_expected = reference.substr(49, 201);
    boundary_expected.replace(30, 1, std::string(1, base(79) == 'C' ? 'G' : 'C') + substitute(base(80)));
    check(boundary.status == AlleleCompositionStatus::Valid && boundary.sequence == boundary_expected,
        "boundary insertion precedes a substitution without shifting its reference position");
    check(compose_allele_sequence(50, 250, {{80, "N", "A"}}, base).status == AlleleCompositionStatus::Unsupported,
        "unknown REF cannot certify composition");
    check(compose_allele_sequence(50, 250, {{80, substitute(base(80)), substitute(base(80))}}, base).status ==
        AlleleCompositionStatus::Unsupported, "reference mismatch rejects even a no-op hypothesis");
    check(compose_allele_sequence(50, 250, {{80, reference.substr(79, 1), "N"}}, base).status ==
        AlleleCompositionStatus::Unsupported, "unknown ALT rejects composition");
    check(compose_allele_sequence(50, 250, {{49, reference.substr(48, 1), "G"}}, base).status ==
        AlleleCompositionStatus::Unsupported, "partial edit coverage cannot certify composition");
    check(compose_allele_sequence(50, 250, {{80, reference.substr(79, 1), std::string(4097, 'A')}}, base).status ==
        AlleleCompositionStatus::Unsupported, "large compositions abstain");
    const auto repeat_base = [](hts_pos_t pos) { return pos > 0 && pos <= 500 ? 'A' : 'N'; };
    check(compose_allele_sequence(1, 20, {{8, "A", "AA"}, {10, "AA", "AAA"}}, repeat_base).status ==
        AlleleCompositionStatus::Overlap, "shifted equivalent descriptions do not invent a second event");
    check(compose_allele_sequence(1, 20, {{8, "A", "AA"}, {8, "AA", "AAA"}}, repeat_base).status ==
        AlleleCompositionStatus::Overlap, "differently padded ambiguous edits cannot certify one event");
    check(compose_allele_sequence(1, 20, {{8, "A", "AA"}, {8, "A", "AA"}}, repeat_base).sequence ==
        std::string(21, 'A'), "exact duplicate anchored repeat edits apply once");
    check(compose_allele_sequence(1, 20, {{8, "A", "C"}, {7, "AA", "AC"}}, repeat_base).sequence ==
        std::string(7, 'A') + 'C' + std::string(12, 'A'), "unambiguous padded duplicates apply once");
    const std::string repeat_context = "CTTTTTTTCTTTCCTTCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGTCTG";
    const auto repeat_context_base = [&](hts_pos_t pos) { return repeat_context.at(pos - 1); };
    check(compose_allele_sequence(1, repeat_context.size(),
        {{1, repeat_context, "CTTTTTTTCTTTCCTTCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGTCTG"},
         {17, "CTTTTTTTTTTTTTTTTTT", "CTTTTTTTTTTTTTTTTTTTT"}}, repeat_context_base).status ==
        AlleleCompositionStatus::Overlap, "parent and nested repeat insertions cannot silently collapse");
    check(compose_allele_sequence(5, 20, {{8, "A", "AA"}}, repeat_base).sequence ==
        std::string(17, 'A'), "composition needs only the complete anchored window");
    check(compose_allele_sequence(1, 20, {{8, "A", ""}, {8, "A", "C"}}, repeat_base).status ==
        AlleleCompositionStatus::Overlap, "left-aligning a deletion cannot unmask a conflicting SNP");
    check(compose_allele_sequence(1, 20, {{8, "A", "AA"}, {8, "A", "C"}}, repeat_base).status ==
        AlleleCompositionStatus::Overlap, "ambiguous anchor trimming cannot choose an insertion/SNP order");
    check(compose_allele_sequence(1, 20, {{1, std::string(20, 'A'), std::string(21, 'A')},
        {10, "A", "C"}}, repeat_base).status == AlleleCompositionStatus::Overlap,
        "a full parent repeat allele cannot invent a placement beside another edit");
    for (int iteration = 0; iteration < 1000; ++iteration) {
        std::vector<AlleleSequenceEdit> edits;
        for (int i = 0; i < 2 + iteration % 4; ++i) {
            const hts_pos_t pos = 50 + i * 30;
            const std::string raw_ref = reference.substr(pos - 1, 1 + random() % 4);
            std::string raw_alt;
            for (size_t j = random() % 6; j > 0; --j) raw_alt += "ACGT"[random() % 4];
            edits.push_back({pos, raw_ref, raw_alt});
        }
        std::string oracle = reference.substr(19, 261);
        for (auto i = edits.rbegin(); i != edits.rend(); ++i)
            oracle.replace(i->pos - 20, i->ref.size(), i->alt);
        const auto result = compose_allele_sequence(20, 280, edits, base);
        check(result.status == AlleleCompositionStatus::Valid && result.sequence == oracle,
            "independent raw-edit oracle across random mixed neighboring alleles");
        edits.push_back(edits.front());
        std::shuffle(edits.begin(), edits.end(), random);
        const auto permuted = compose_allele_sequence(20, 280, edits, base);
        check(permuted.status == result.status && permuted.sequence == oracle,
            "composition is independent of input ordering and duplicate descriptions");
    }
    const std::string path_reference = reference.substr(49, 201);
    const std::vector<AlleleNeighborSite> path_sites = {
        {20, 140, reference.substr(139, 1), {substitute(base(140))}},
        {10, 80, reference.substr(79, 1), {substitute(base(80)), reference.substr(79, 1) + "GG"}}};
    const auto paths = build_allele_sequence_paths(50, 250, {path_reference}, path_sites, base);
    check(paths.status == AllelePathStatus::Complete && paths.paths.size() == 6,
        "all one-ALT-per-site and omitted-site combinations are represented");
    std::set<std::string> path_sequences;
    for (const auto& path : paths.paths) {
        std::string oracle = path_reference;
        for (auto selected = path.neighbors.rbegin(); selected != path.neighbors.rend(); ++selected) {
            const auto site = std::find_if(path_sites.begin(), path_sites.end(),
                [&](const auto& s) { return s.candidate == selected->first; });
            oracle.replace(site->pos - 50, site->ref.size(), site->alts.at(selected->second - 1));
        }
        check(path.sequence == oracle && path.parent_allele == 0, "multi-neighbor raw edit oracle");
        path_sequences.insert(path.sequence);
    }
    check(path_sequences.size() == 6 && paths.overlap_prefixes == 0 && paths.unsupported_paths == 0,
        "disjoint edits preserve every distinct multi-neighbor sequence");
    const auto limited = build_allele_sequence_paths(50, 250, {path_reference}, path_sites, base,
        paths.visited_prefixes - 1);
    check(limited.status == AllelePathStatus::Limited && limited.paths.empty(),
        "work limits cannot expose a truncated catalog as complete");
    const auto exact_budget = build_allele_sequence_paths(50, 250, {path_reference}, path_sites, base,
        paths.visited_prefixes);
    check(exact_budget.status == AllelePathStatus::Complete && exact_budget.paths.size() == paths.paths.size(),
        "exactly enough prefix budget completes the search");
    check(build_allele_sequence_paths(50, 250, {path_reference}, path_sites, base, 0).status ==
        AllelePathStatus::Limited, "zero work budget explicitly abstains");
    auto duplicate_sites = path_sites;
    duplicate_sites.push_back(path_sites.front());
    check(build_allele_sequence_paths(50, 250, {path_reference}, duplicate_sites, base).status ==
        AllelePathStatus::Unsupported, "duplicate candidate IDs cannot select two ALTs of one site");
    check(build_allele_sequence_paths(0, 250, {path_reference}, {}, base).status ==
        AllelePathStatus::Unsupported, "invalid reference window cannot yield paths");
    const auto unsupported_parent = build_allele_sequence_paths(50, 250, {path_reference, "N"}, {}, base);
    check(unsupported_parent.status == AllelePathStatus::Complete && unsupported_parent.paths.size() == 1 &&
        unsupported_parent.unsupported_paths == 1, "unsupported parent paths retain explicit uncertainty");
    const auto independent_parents = build_allele_sequence_paths(50, 250,
        {path_reference, path_reference}, path_sites, base);
    check(independent_parents.paths.size() == 12 && independent_parents.paths[6].parent_allele == 1,
        "identical sequences retain independent parent provenance without changing their edits");
    std::vector<AlleleNeighborSite> deep_sites;
    for (size_t i = 0; i < 65; ++i) deep_sites.push_back({i, 80 + static_cast<hts_pos_t>(i), reference.substr(79 + i, 1), {}});
    check(build_allele_sequence_paths(50, 250, {path_reference}, deep_sites, base).status ==
        AllelePathStatus::Limited, "diagnostic recursion depth is bounded and marked incomplete");
    const auto conflicting_paths = build_allele_sequence_paths(1, 20, {std::string(20, 'A')},
        {{1, 8, "A", {"C"}}, {2, 8, "A", {"G"}}}, repeat_base);
    check(conflicting_paths.status == AllelePathStatus::Complete && conflicting_paths.paths.size() == 3 &&
        conflicting_paths.overlap_prefixes == 1, "conflicting ALT selections prune only their descendants");
    const auto duplicate_paths = build_allele_sequence_paths(1, 20, {std::string(20, 'A')},
        {{1, 8, "A", {"C"}}, {2, 8, "A", {"C"}}}, repeat_base);
    check(duplicate_paths.paths.size() == 2 && duplicate_paths.paths.back().neighbors.size() == 2 &&
        duplicate_paths.paths.back().sequence ==
        std::string(7, 'A') + 'C' + std::string(12, 'A'),
        "two source descriptions retain provenance while applying one edit");
    const auto atomic_aliases = build_allele_sequence_paths(1, 20, {std::string(20, 'A')},
        {{1, 8, "AA", {"CA", "AC"}}, {2, 8, "AA", {"AC", "CA"}}}, repeat_base);
    check(atomic_aliases.paths.size() == 3 && std::none_of(atomic_aliases.paths.begin(), atomic_aliases.paths.end(),
        [](const auto& path) { return path.sequence.find("CC") != std::string::npos; }),
        "duplicate complete sites cannot combine mutually exclusive physical ALTs");
    check(atomic_aliases.paths.back().neighbors == std::vector<std::pair<size_t, size_t>>{{1, 1}, {2, 2}},
        "physical alias selection maps back to each original ALT table order");
    const auto different_tables = build_allele_sequence_paths(1, 20, {std::string(20, 'A')},
        {{1, 8, "AA", {"CA"}}, {2, 8, "AA", {"AC"}}}, repeat_base);
    check(different_tables.paths.size() == 4 && different_tables.paths.back().sequence.find("CC") != std::string::npos,
        "different complete tables are not silently merged by position alone");
    const auto repaired_length = build_allele_sequence_paths(1, 4096, {std::string(4096, 'A')},
        {{1, 100, "", {"C"}}, {2, 200, "A", {""}}}, [](hts_pos_t) { return 'A'; });
    check(repaired_length.status == AllelePathStatus::Complete && repaired_length.paths.size() == 3 &&
        repaired_length.unsupported_paths == 1 && repaired_length.paths.back().neighbors.size() == 2 &&
        repaired_length.paths.back().sequence.size() == 4096,
        "overlong intermediate insertion is not pruned before a compensating deletion");
    for (int iteration = 0; iteration < 100; ++iteration) {
        std::vector<AlleleNeighborSite> sites;
        for (size_t i = 0; i < 4; ++i) {
            const hts_pos_t pos = 70 + i * 30 + random() % 5;
            const std::string site_ref(1, base(pos));
            const std::string inserted{ "ACGT"[random() % 4], "ACGT"[random() % 4] };
            sites.push_back({i, pos, site_ref, {substitute(base(pos)), site_ref + inserted}});
        }
        const auto all = build_allele_sequence_paths(50, 250, {path_reference}, sites, base);
        check(all.status == AllelePathStatus::Complete && all.paths.size() == 81,
            "four three-way choices enumerate every complete assignment");
        std::set<std::string> oracle_sequences;
        for (size_t assignment = 0; assignment < 81; ++assignment) {
            size_t value = assignment;
            std::string oracle = path_reference;
            std::vector<size_t> choices;
            for (size_t i = 0; i < sites.size(); ++i) { choices.push_back(value % 3); value /= 3; }
            for (size_t i = sites.size(); i-- > 0;)
                if (choices[i]) oracle.replace(sites[i].pos - 50, sites[i].ref.size(),
                    sites[i].alts[choices[i] - 1]);
            oracle_sequences.insert(oracle);
        }
        std::set<std::string> observed;
        for (const auto& path : all.paths) observed.insert(path.sequence);
        check(observed == oracle_sequences, "exhaustive independent raw multi-site sequence oracle");
        std::shuffle(sites.begin(), sites.end(), random);
        const auto shuffled = build_allele_sequence_paths(50, 250, {path_reference}, sites, base);
        check(shuffled.paths.size() == all.paths.size() && shuffled.visited_prefixes == all.visited_prefixes,
            "neighbor order cannot change paths or deterministic work limits");
        for (size_t i = 0; i < all.paths.size(); ++i)
            check(shuffled.paths[i].neighbors == all.paths[i].neighbors &&
                shuffled.paths[i].sequence == all.paths[i].sequence,
                "path ordering and provenance survive neighbor permutations");
    }
    const std::string ref(1, base(100));
    const std::string alt(1, base(100) == 'C' ? 'G' : 'C');
    const std::string other(1, base(100) != 'A' && alt != "A" ? 'A' : 'T');
    const auto pair = normalize_allele_contrast(100, ref, ref, alt, base);
    check(pair.has_value(), "normalize fixture");
    const auto context = build_allele_sequence_context(pair->key, {{100, ref, {alt, other}}}, base);
    check(context && context->beg == 84 && context->end == 116, "complete context and outside flanks");
    if (!context) return 1;
    check(context->other_alleles.size() == 1, "retain third parent allele");
    for (int allele = 0; allele < 2; ++allele) {
        const auto score = score_allele_context(*context, context->alleles[allele]);
        check(score.nearest == allele && score.distances[allele] == 0, "rank selected full allele");
        for (int flag : {0, BAM_FREVERSE}) {
            auto bam = read(context->beg, "33M", context->alleles[allele], flag);
            const auto slice = extract_allele_read_slice(bam.get(), *context);
            check(slice && slice->sequence == context->alleles[allele] && slice->query_beg == 0 &&
                slice->query_end == 33 && slice->qualities == std::vector<uint8_t>(33, 40) &&
                slice->mapq == 60, "original sequence, quality and reverse alignment orientation");
        }
    }
    std::string other_query = reference.substr(83, 33);
    other_query[16] = other[0];
    const auto other_score = score_allele_context(*context, other_query);
    check(other_score.nearest == -1 && other_score.other_distance == 0,
        "short REF/ALT cannot assign a different full parent allele");
    const auto alt_pair = normalize_allele_contrast(100, ref, alt, other, base);
    const auto alt_context = build_allele_sequence_context(alt_pair->key, {{100, ref, {other, alt}}}, base);
    check(alt_context && alt_context->other_alleles.size() == 1 &&
        score_allele_context(*alt_context, reference.substr(83, 33)).nearest == -1,
        "literal REF remains an unselected hypothesis at ALT/ALT");
    const auto reversed_pair = normalize_allele_contrast(100, ref, alt, ref, base);
    const auto reversed = build_allele_sequence_context(reversed_pair->key, {{100, ref, {alt, other}}}, base);
    check(reversed && reversed->alleles == context->alleles && reversed->other_alleles == context->other_alleles,
        "source allele order cannot alter physical hypotheses");
    auto reordered_context = *alt_context;
    std::swap(reordered_context.alleles[0], reordered_context.alleles[1]);
    reordered_context.other_alleles.push_back(reordered_context.other_alleles.front());
    auto shifted_context = *context;
    ++shifted_context.beg;
    ++shifted_context.end;
    auto longer_context = *context;
    ++longer_context.end;
    auto different_hypotheses = *context;
    different_hypotheses.other_alleles.push_back("GGG");
    const std::vector<std::optional<AlleleSequenceContext>> descriptions{
        context, alt_context, reversed, reordered_context, std::nullopt,
        shifted_context, longer_context, different_hypotheses};
    const auto physical_groups = group_allele_contexts(descriptions);
    check(physical_groups.size() == 4, "physical identity includes both bounds and complete allele table");
    check(physical_groups.at({context->beg, context->end, full_allele_sequences(*context)}) ==
        std::vector<size_t>({0, 1, 2, 3}),
        "REF/ALT and ALT/ALT descriptions, allele reversals and duplicate hypotheses share one physical group");
    check(group_allele_contexts({std::nullopt}).empty(), "unsupported contexts cannot share an identity");
    for (int iteration = 0; iteration < 100; ++iteration) {
        std::vector<size_t> order(descriptions.size());
        std::iota(order.begin(), order.end(), 0);
        std::shuffle(order.begin(), order.end(), random);
        std::vector<std::optional<AlleleSequenceContext>> permuted;
        for (const size_t i : order) permuted.push_back(descriptions[i]);
        const auto groups = group_allele_contexts(permuted);
        check(groups.size() == physical_groups.size(), "physical count is independent of source enumeration");
        auto expected = physical_groups.begin();
        for (const auto& [key, members] : groups) {
            std::vector<size_t> original_members;
            for (const size_t i : members) original_members.push_back(order[i]);
            std::sort(original_members.begin(), original_members.end());
            check(expected != physical_groups.end() && !(key < expected->first) &&
                !(expected->first < key) && original_members == expected->second,
                "physical identities, ordering and memberships survive source permutations");
            if (expected != physical_groups.end()) ++expected;
        }
    }
    const auto padded = build_allele_sequence_context(pair->key,
        {{98, reference.substr(97, 5), {reference.substr(97, 2) + alt + reference.substr(100, 2)}},
         {100, ref, {alt, other}}}, base);
    check(padded && padded->beg == 82 && padded->end == 118 && padded->other_alleles.size() == 1,
        "padded descriptions extend full coverage without duplicate hypotheses");
    const std::string parent_ref = reference.substr(97, 5);
    std::string parent_other = parent_ref;
    parent_other[4] = parent_ref[4] == 'A' ? 'C' : 'A';
    const auto parent = build_allele_sequence_context(pair->key,
        {{98, parent_ref, {parent_ref.substr(0, 2) + alt + parent_ref.substr(3), parent_other}},
         {100, ref, {alt}}}, base);
    check(parent.has_value(), "build full parent around primitive marker");
    if (parent) {
        std::string query = reference.substr(parent->beg - 1, parent->end - parent->beg + 1);
        query.replace(98 - parent->beg, parent_ref.size(), parent_other);
        const auto score = score_allele_context(*parent, query);
        check(query[100 - parent->beg] == ref[0] && score.distances[0] > 0 &&
            score.other_distance == 0 && score.nearest == -1,
            "primitive REF cannot certify full parent REF when another branch changes elsewhere");
    }
    check(!build_allele_sequence_context(pair->key, {{100, "N", {alt}}}, base), "unknown REF rejects context");
    check(!build_allele_sequence_context(pair->key, {{100, ref, {"N"}}}, base), "unknown parent ALT rejects context");
    check(!build_allele_sequence_context(pair->key, {{100, ref, {""}}}, base), "empty VCF ALT rejects context");
    check(!build_allele_sequence_context(pair->key, {{1, std::string(1, base(1)), {"C"}}}, base), "missing flank rejects context");
    check(!build_allele_sequence_context(pair->key, {{100, ref, {std::string(4097, 'A')}}}, base), "large hypotheses abstain");
    auto bam = read(84, "33M", context->alleles[0]);
    bam_get_qual(bam.get())[8] = 255;
    bam_get_qual(bam.get())[9] = 1;
    const auto quality_slice = extract_allele_read_slice(bam.get(), *context);
    check(quality_slice && quality_slice->qualities[8] == 255 && quality_slice->qualities[9] == 1,
        "unknown and low qualities retained rather than fabricated confidence");
    bam = read(85, "32M", context->alleles[0].substr(1));
    check(!extract_allele_read_slice(bam.get(), *context), "partial left flank abstains");
    bam = read(84, "32M", context->alleles[0].substr(0, 32));
    check(!extract_allele_read_slice(bam.get(), *context), "partial right flank abstains");
    bam = read(84, "16M1N16M", context->alleles[0].substr(0, 16) + context->alleles[0].substr(17));
    check(!extract_allele_read_slice(bam.get(), *context), "reference skip cannot certify parent");
    bam = read(84, "1D32M", context->alleles[0].substr(1));
    check(!extract_allele_read_slice(bam.get(), *context), "deleted outer anchor abstains");
    bam = read(84, "2S33M2S", "GG" + context->alleles[0] + "CC");
    const auto clipped = extract_allele_read_slice(bam.get(), *context);
    check(clipped && clipped->query_beg == 2 && clipped->query_end == 35 &&
        clipped->sequence == context->alleles[0], "soft clips outside complete context excluded");
    bam = read(84, "33M", std::string(33, 'N'));
    check(!extract_allele_read_slice(bam.get(), *context), "unknown query bases abstain");
    check(!extract_allele_read_slice(nullptr, *context), "missing alignment abstains");
    auto selected = read(84, "33M", context->alleles[0], 0, "selected");
    auto unobserved = read(84, "33M", other_query, BAM_FREVERSE, "unobserved");
    auto collision = read(85, "32M", context->alleles[0].substr(1), 0, "collision");
    auto collision_covering = read(84, "33M", context->alleles[1], 0, "collision");
    auto secondary = read(84, "33M", context->alleles[1], BAM_FSECONDARY, "selected");
    auto supplementary = read(84, "33M", context->alleles[1], BAM_FSUPPLEMENTARY, "selected");
    auto unmapped = read(84, "33M", context->alleles[1], BAM_FUNMAP, "unmapped");
    auto no_contig = read(84, "33M", context->alleles[1], 0, "no_contig");
    no_contig->core.tid = -1;
    auto duplicate = read(84, "33M", context->alleles[1], BAM_FDUP, "duplicate");
    auto qcfail = read(84, "33M", context->alleles[1], BAM_FQCFAIL, "qcfail");
    auto low_mapq = read(84, "33M", context->alleles[1], 0, "low_mapq");
    low_mapq->core.qual = 0;
    auto skip = read(84, "16M1N16M", other_query.substr(0, 16) + other_query.substr(17), 0, "skip");
    auto deleted = read(84, "1D32M", other_query.substr(1), 0, "deleted");
    auto right_partial = read(84, "32M", other_query.substr(0, 32), 0, "right_partial");
    bam_get_qual(unobserved.get())[4] = 255;
    std::vector<const bam1_t*> alignments{selected.get(), unobserved.get(), collision.get(),
        collision_covering.get(), secondary.get(), supplementary.get(), unmapped.get(),
        duplicate.get(), qcfail.get(), low_mapq.get(), skip.get(), deleted.get(), right_partial.get(),
        no_contig.get(), nullptr};
    const auto cohort = collect_allele_read_slices(alignments, *context, 1, false);
    check(cohort.size() == 2 && cohort.count("selected") && cohort.count("unobserved"),
        "complete primary cohort includes other alleles without source observations");
    check(cohort.at("unobserved").sequence == other_query &&
        cohort.at("unobserved").qualities[4] == 255 && cohort.at("unobserved").mapq == 60,
        "cohort preserves original reverse sequence, unknown quality and MAPQ");
    check(!cohort.count("collision"), "duplicate primary ambiguity precedes physical coverage");
    check(collect_allele_read_slices(alignments, *context, 0, false).count("low_mapq") == 1,
        "cohort honors configured recovery MAPQ floor");
    const auto filtered = collect_allele_read_slices(alignments, *context, 1, true);
    check(filtered.size() == 4 && filtered.count("duplicate") && filtered.count("qcfail"),
        "explicit filtered-read admission matches BAM loading policy");
    check(collect_allele_read_slices(alignments, *context, 61, true).empty(),
        "MAPQ above all alignments produces empty cohort");
    for (int iteration = 0; iteration < 100; ++iteration) {
        std::shuffle(alignments.begin(), alignments.end(), random);
        const auto reordered = collect_allele_read_slices(alignments, *context, 1, false);
        check(reordered.size() == cohort.size() && reordered.count("selected") &&
            reordered.count("unobserved") && reordered.at("unobserved").sequence == other_query,
            "cohort is invariant to primary, secondary and conflicting alignment order");
    }
    alignments.push_back(collision_covering.get());
    check(!collect_allele_read_slices(alignments, *context, 1, false).count("collision"),
        "third primary cannot revive ambiguous molecule");
    check(collect_allele_read_slices({}, *context, 1, false).empty(), "empty BAM cohort abstains");
    const auto insertion_pair = normalize_allele_contrast(100, ref, ref, ref + "GG", base);
    const auto insertion = build_allele_sequence_context(insertion_pair->key, {{100, ref, {ref + "GG"}}}, base);
    check(insertion.has_value(), "build insertion context");
    if (insertion) {
        const auto edit = *insertion_pair->key.alleles[1];
        const int left = edit.pos + 1 - insertion->beg;
        const int right = insertion->end - edit.pos;
        bam = read(insertion->beg, std::to_string(left) + "M2I" + std::to_string(right) + "M", insertion->alleles[1]);
        const auto slice = extract_allele_read_slice(bam.get(), *insertion);
        check(slice && slice->sequence == insertion->alleles[1] &&
            score_allele_context(*insertion, slice->sequence).nearest == 1, "internal insertion and motif shift retained");
    }
    const std::string deletion_ref = reference.substr(99, 3);
    const auto deletion_pair = normalize_allele_contrast(100, deletion_ref, deletion_ref, ref, base);
    const auto deletion = build_allele_sequence_context(deletion_pair->key, {{100, deletion_ref, {ref}}}, base);
    check(deletion.has_value(), "build deletion context");
    if (deletion) {
        const auto edit = *deletion_pair->key.alleles[1];
        const int left = edit.pos + 1 - deletion->beg;
        const int right = deletion->end - edit.pos - 2;
        bam = read(deletion->beg, std::to_string(left) + "M2D" + std::to_string(right) + "M", deletion->alleles[1]);
        const auto slice = extract_allele_read_slice(bam.get(), *deletion);
        check(slice && slice->sequence == deletion->alleles[1] &&
            score_allele_context(*deletion, slice->sequence).nearest == 1, "internal deletion retained");
    }
    std::vector<AlleleSequenceContext> contexts{*context, *alt_context, *padded};
    if (insertion) contexts.push_back(*insertion);
    if (deletion) contexts.push_back(*deletion);
    for (const auto& c : contexts) {
        for (int iteration = 0; iteration < 100; ++iteration) {
            std::string query = c.alleles[random() % 2];
            for (int j = 0; j < iteration % 5; ++j) query[random() % query.size()] = "ACGT"[random() % 4];
            if (iteration % 3 == 0) query.insert(random() % query.size(), 1, 'C');
            if (iteration % 7 == 0) query.erase(random() % query.size(), 1);
            const auto score = score_allele_context(c, query);
            const std::array<int, 2> expected{distance(query, c.alleles[0]), distance(query, c.alleles[1])};
            int other_distance = -1;
            for (const auto& other_allele : c.other_alleles) {
                const int d = distance(query, other_allele);
                if (other_distance < 0 || d < other_distance) other_distance = d;
            }
            const int winner = expected[0] < expected[1] ? 0 : 1;
            const int nearest = expected[winner] < expected[1 - winner] &&
                (other_distance < 0 || expected[winner] < other_distance) ? winner : -1;
            check(score.distances == expected && score.other_distance == other_distance &&
                score.nearest == nearest, "independent edit-distance and decoy oracle");
        }
    }
    std::cout << checks << " checks, " << failures << " failures\n";
    return failures != 0;
}
