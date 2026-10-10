#include "allele_genotype.hpp"
#include "allele_context.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <tuple>

using namespace pgphase_collect;

static int failures = 0;
static size_t checks = 0;
static void check(bool ok, const char* message) {
    ++checks;
    if (!ok) { ++failures; std::cerr << "FAIL: " << message << '\n'; }
}

static DiploidAlleleFit oracle(size_t n, const std::vector<MoleculeAlleleCosts>& rows) {
    if (rows.empty() || n == 0) return {};
    std::vector<std::pair<int64_t, std::array<size_t, 2>>> pairs;
    for (size_t a = 0; a < n; ++a) for (size_t b = a; b < n; ++b) {
        int64_t score = 0;
        for (const auto& row : rows) score += std::min(row.costs[a], row.costs[b]);
        pairs.push_back({score, {a, b}});
    }
    std::sort(pairs.begin(), pairs.end());
    DiploidAlleleFit result;
    result.cost = pairs[0].first;
    if (pairs.size() > 1) result.runner_up_cost = pairs[1].first;
    result.tied_pairs = std::count_if(pairs.begin(), pairs.end(), [&](const auto& pair) { return pair.first == *result.cost; });
    if (result.tied_pairs == 1) result.alleles = pairs[0].second;
    return result;
}

static bool equal(const DiploidAlleleFit& a, const DiploidAlleleFit& b) {
    return a.alleles == b.alleles && a.cost == b.cost && a.runner_up_cost == b.runner_up_cost && a.tied_pairs == b.tied_pairs;
}

static std::vector<std::string> split(const std::string& line, char delimiter) {
    std::vector<std::string> result;
    size_t beg = 0;
    for (;;) {
        const size_t end = line.find(delimiter, beg);
        result.push_back(line.substr(beg, end - beg));
        if (end == std::string::npos) return result;
        beg = end + 1;
    }
}

static void write_fit(std::ostream& out, const DiploidAlleleFit& fit) {
    if (fit.cost) out << *fit.cost; else out << '.';
    out << '\t';
    if (fit.runner_up_cost) out << *fit.runner_up_cost; else out << '.';
    out << '\t' << fit.tied_pairs << '\t';
    if (fit.alleles) out << (*fit.alleles)[0] << '\t' << (*fit.alleles)[1];
    else out << ".\t.";
}

// Rescore one saved cohort per physical context, then project compatibility views.
static int replay(const std::string& folder, const std::string& prefix) {
    std::ifstream contexts_file(folder + "/matrix.chunk0.joint-contexts.tsv");
    std::ifstream sequences_file(folder + "/matrix.chunk0.physical-cohort.tsv");
    if (!contexts_file || !sequences_file) { std::cerr << "cannot read saved physical context in " << folder << '\n'; return 1; }
    std::ofstream alleles(prefix + ".alleles.tsv"), costs(prefix + ".costs.tsv"), genotypes(prefix + ".genotypes.tsv"), heldout(prefix + ".heldout.tsv");
    std::ofstream physical_contexts(prefix + ".physical-contexts.tsv"), physical_members(prefix + ".physical-members.tsv");
    std::ofstream physical_alleles(prefix + ".physical-alleles.tsv"), physical_costs(prefix + ".physical-costs.tsv");
    std::ofstream physical_genotypes(prefix + ".physical-genotypes.tsv"), physical_heldout(prefix + ".physical-heldout.tsv");
    if (!alleles || !costs || !genotypes || !heldout || !physical_contexts || !physical_members ||
        !physical_alleles || !physical_costs || !physical_genotypes || !physical_heldout) {
        std::cerr << "cannot write saved physical state: " << prefix << '\n'; return 1;
    }
    for (auto* stream : {&alleles, &physical_alleles})
        *stream << (stream == &alleles ? "locus" : "physical") << "\tallele\tsequence\n";
    for (auto* stream : {&costs, &physical_costs})
        *stream << (stream == &costs ? "locus" : "physical") << "\tread\tcosts\n";
    for (auto* stream : {&genotypes, &physical_genotypes})
        *stream << (stream == &genotypes ? "locus" : "physical") << "\tmolecules\tcost\trunner_up\ttied_pairs\tallele0\tallele1\tconflicting_molecules\n";
    for (auto* stream : {&heldout, &physical_heldout})
        *stream << (stream == &heldout ? "locus" : "physical") << "\tread\tcost\trunner_up\ttied_pairs\tallele0\tallele1\tallele\tsame_pair\n";
    physical_contexts << "physical\tbeg\tend\n";
    physical_members << "physical\tlocus\tallele0\tallele1\n";
    std::vector<std::optional<AlleleSequenceContext>> contexts;
    std::string line;
    std::getline(contexts_file, line);
    while (std::getline(contexts_file, line)) {
        const auto f = split(line, '\t');
        if (f.size() != 8) return 1;
        const size_t locus = std::stoull(f[0]);
        if (contexts.size() <= locus) contexts.resize(locus + 1);
        if (f[1] == "valid") contexts[locus] = AlleleSequenceContext{std::stoll(f[2]), std::stoll(f[3]),
            {f[4], f[5]}, f[6].empty() ? std::vector<std::string>{} : split(f[6], ',')};
    }
    std::vector<std::vector<std::string>> hypotheses;
    std::vector<std::optional<size_t>> physical_by_locus(contexts.size());
    for (const auto& [key, members] : group_allele_contexts(contexts)) {
        const size_t id = hypotheses.size();
        hypotheses.push_back(key.alleles);
        physical_contexts << id << '\t' << key.beg << '\t' << key.end << '\n';
        for (const size_t member : members) {
            physical_by_locus[member] = id;
            physical_members << id << '\t' << member;
            for (const auto& allele : contexts[member]->alleles)
                physical_members << '\t' << std::distance(key.alleles.begin(),
                    std::lower_bound(key.alleles.begin(), key.alleles.end(), allele));
            physical_members << '\n';
        }
    }
    std::vector<std::vector<MoleculeAlleleCosts>> observations(hypotheses.size());
    std::getline(sequences_file, line);
    size_t count = 0;
    while (std::getline(sequences_file, line)) {
        const auto f = split(line, '\t');
        if (f.size() != 7) return 1;
        const size_t id = std::stoull(f[0]);
        observations.at(id).push_back({f[1], score_allele_sequences(hypotheses.at(id), f[5])});
        ++count;
    }
    std::vector<DiploidAlleleModel> models;
    const auto write_state = [&](size_t id, const std::vector<std::string>& sequences,
            const std::vector<MoleculeAlleleCosts>& rows, const DiploidAlleleModel& model,
            std::ostream& allele_out, std::ostream& cost_out, std::ostream& genotype_out, std::ostream& heldout_out) {
        for (size_t i = 0; i < sequences.size(); ++i) allele_out << id << '\t' << i << '\t' << sequences[i] << '\n';
        for (const auto& row : rows) {
            cost_out << id << '\t' << row.molecule << '\t';
            for (size_t i = 0; i < row.costs.size(); ++i) { if (i > 0) cost_out << ','; cost_out << row.costs[i]; }
            cost_out << '\n';
        }
        genotype_out << id << '\t' << model.heldout.size() << '\t';
        write_fit(genotype_out, model.fit);
        genotype_out << '\t' << model.conflicting_molecules.size() << '\n';
        for (const auto& read : model.heldout) {
            heldout_out << id << '\t' << read.molecule << '\t';
            write_fit(heldout_out, read.fit);
            heldout_out << '\t' << read.allele << '\t' << read.same_pair << '\n';
        }
    };
    for (size_t id = 0; id < hypotheses.size(); ++id) {
        models.push_back(fit_diploid_alleles(hypotheses[id].size(), observations[id]));
        check(models.back().heldout.size() == observations[id].size(), "physical state has one row per eligible molecule");
        write_state(id, hypotheses[id], observations[id], models.back(),
            physical_alleles, physical_costs, physical_genotypes, physical_heldout);
    }
    for (size_t locus = 0; locus < contexts.size(); ++locus) {
        if (!physical_by_locus[locus]) {
            write_state(locus, {}, {}, {}, alleles, costs, genotypes, heldout);
            continue;
        }
        const size_t id = *physical_by_locus[locus];
        write_state(locus, hypotheses[id], observations[id], models[id], alleles, costs, genotypes, heldout);
    }
    check(count > 0, "saved-state replay cannot silently skip every molecule");
    std::cout << hypotheses.size() << " physical contexts, " << count << " scored molecules\n";
    return failures != 0;
}

int main(int argc, char** argv) {
    if (argc == 4 && std::string(argv[1]) == "--state") return replay(argv[2], argv[3]);
    if (argc != 1) { std::cerr << "Usage: test_allele_genotype [--state INPUT_FOLDER OUTPUT_PREFIX]\n"; return 1; }
    const AlleleSequenceContext c{1, 3, {"CCC", "AAA"}, {"GGG", "CCC", "GGG"}};
    const auto hypotheses = full_allele_sequences(c);
    check(hypotheses == std::vector<std::string>({"AAA", "CCC", "GGG"}), "all hypotheses are unique and ordered independently of selection");
    check(score_allele_sequences(hypotheses, "GGG") == std::vector<int>({3, 3, 0}), "retain individual unselected allele costs");
    const auto alternatives = full_allele_sequences({1, 3, {"AAA", "CCC"}, {"TTT", "GGG"}});
    std::vector<MoleculeAlleleCosts> missing_pair;
    for (const auto& read : std::vector<std::pair<std::string, std::string>>{
            {"a", "GGG"}, {"b", "GGG"}, {"c", "TTT"}, {"d", "TTT"}})
        missing_pair.push_back({read.first, score_allele_sequences(alternatives, read.second)});
    const auto recovered = fit_diploid_alleles(alternatives.size(), missing_pair);
    check(recovered.fit.alleles == std::optional<std::array<size_t, 2>>({{2, 3}}) &&
        std::all_of(recovered.heldout.begin(), recovered.heldout.end(), [](const auto& r) { return r.same_pair; }),
        "complete parent sequences can support a held-out ALT/ALT pair absent from the original contrast");
    std::vector<MoleculeAlleleCosts> balanced{{"a", {0, 5, 9}}, {"b", {0, 5, 9}}, {"c", {5, 0, 9}}, {"d", {5, 0, 9}}};
    const auto model = fit_diploid_alleles(3, balanced);
    check(model.fit.alleles == std::optional<std::array<size_t, 2>>({{0, 1}}) && model.fit.cost == 0 && model.fit.runner_up_cost == 10, "infer diploid pair from all molecules");
    check(std::all_of(model.heldout.begin(), model.heldout.end(), [](const auto& r) { return r.same_pair && r.allele >= 0; }), "balanced pair supported independently of each held-out molecule");
    balanced.push_back({"outside", {1, 2, 0}});
    const auto outside = fit_diploid_alleles(3, balanced);
    const auto& out = outside.heldout.back();
    check(out.molecule == "outside" && out.same_pair && out.allele == -1, "a closer unselected hypothesis cannot become selected-pair evidence");
    std::vector<MoleculeAlleleCosts> self;
    for (int i = 0; i < 8; ++i) self.push_back({std::to_string(i), {0, 5, 9}});
    self.push_back({"own", {5, 0, 9}});
    const auto singleton = fit_diploid_alleles(3, self);
    check(singleton.fit.alleles.has_value() && !singleton.heldout.back().fit.alleles &&
        !singleton.heldout.back().same_pair && singleton.heldout.back().allele == -1,
        "one molecule cannot certify its own second allele");
    const auto homo = fit_diploid_alleles(1, {{"a", {0}}, {"b", {0}}});
    check(homo.fit.alleles == std::optional<std::array<size_t, 2>>({{0, 0}}), "homozygous hypotheses are included");
    const auto tie = fit_diploid_alleles(3, {{"a", {0, 0, 0}}, {"b", {0, 0, 0}}});
    check(!tie.fit.alleles && tie.fit.tied_pairs == 6 && tie.fit.runner_up_cost == 0, "ties cannot borrow first ALT or invent a genotype");
    const auto single = fit_diploid_alleles(3, {{"a", {0, 1, 2}}});
    check(!single.heldout[0].fit.cost && single.heldout[0].allele == -1, "empty held-out training data abstains");
    auto duplicates = balanced;
    duplicates.insert(duplicates.end(), balanced.begin(), balanced.end());
    const auto repeated = fit_diploid_alleles(3, duplicates);
    check(equal(repeated.fit, outside.fit) && repeated.heldout.size() == outside.heldout.size(), "duplicate source descriptions do not change costs or sample size");
    duplicates.push_back({"outside", {0, 2, 2}});
    duplicates.push_back({"outside", {1, 2, 0}});
    const auto conflict = fit_diploid_alleles(3, duplicates);
    check(conflict.conflicting_molecules == std::vector<std::string>{"outside"} && conflict.heldout.size() == 4 &&
        equal(conflict.fit, model.fit), "conflicting molecule costs remain sticky unknown rather than extra votes");
    const auto huge = fit_diploid_alleles(2, {{"a", {0, 2000000000}}, {"b", {0, 2000000000}},
        {"c", {2000000000, 0}}, {"d", {2000000000, 0}}});
    check(huge.fit.cost == 0 && huge.fit.runner_up_cost == 4000000000LL, "cohort costs retain 64-bit sums");
    for (const auto& invalid : {std::vector<int>{0}, std::vector<int>{0, -1, 2}}) {
        bool rejected = false;
        try { (void)fit_diploid_alleles(3, {{"invalid", invalid}}); }
        catch (const std::invalid_argument&) { rejected = true; }
        check(rejected, "incomplete or negative cost tables cannot be padded into evidence");
    }
    std::mt19937 random(417);
    for (int iteration = 0; iteration < 1000; ++iteration) {
        const size_t n = 1 + random() % 6;
        std::vector<MoleculeAlleleCosts> rows;
        for (size_t i = 0, count = random() % 9; i < count; ++i) {
            std::vector<int> costs(n);
            for (int& cost : costs) cost = random() % 8;
            rows.push_back({std::to_string(i), costs});
        }
        const auto fit = fit_diploid_alleles(n, rows);
        check(equal(fit.fit, oracle(n, rows)), "exhaustive independently sorted full-cohort oracle");
        check(fit.heldout.size() == rows.size(), "one held-out fit per original molecule");
        for (size_t i = 0; i < rows.size(); ++i) {
            auto training = rows;
            training.erase(training.begin() + i);
            const auto expected = oracle(n, training);
            const auto& actual = fit.heldout[i];
            check(equal(actual.fit, expected), "held-out fit matches fresh brute-force training solve");
            const auto& costs = rows[i].costs;
            const int best = *std::min_element(costs.begin(), costs.end());
            const size_t allele = std::find(costs.begin(), costs.end(), best) - costs.begin();
            const int call = expected.alleles && std::count(costs.begin(), costs.end(), best) == 1 &&
                (allele == (*expected.alleles)[0] || allele == (*expected.alleles)[1]) ? static_cast<int>(allele) : -1;
            check(actual.allele == call && actual.same_pair == (fit.fit.alleles && expected.alleles == fit.fit.alleles), "held-out membership retains competition and pair uncertainty");
        }
        std::reverse(rows.begin(), rows.end());
        const auto reordered = fit_diploid_alleles(n, rows);
        check(equal(fit.fit, reordered.fit), "read order cannot select a different genotype");
        std::vector<size_t> permutation(n);
        for (size_t i = 0; i < n; ++i) permutation[i] = i;
        std::shuffle(permutation.begin(), permutation.end(), random);
        for (auto& row : rows) {
            const auto previous = row.costs;
            for (size_t i = 0; i < n; ++i) row.costs[i] = previous[permutation[i]];
        }
        auto permuted = fit_diploid_alleles(n, rows).fit;
        if (permuted.alleles) {
            auto& pair = *permuted.alleles;
            pair = {permutation[pair[0]], permutation[pair[1]]};
            std::sort(pair.begin(), pair.end());
        }
        check(equal(fit.fit, permuted), "allele-index permutation changes only the coordinate gauge");
    }
    std::cout << checks << " checks, " << failures << " failures\n";
    return failures != 0;
}
