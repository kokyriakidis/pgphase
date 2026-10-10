// Evaluation-only read-hypothesis residuals; original catalog costs stay fixed.
#include "allele_context.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

static std::vector<std::string> fields(const std::string& line) {
    std::vector<std::string> result;
    std::istringstream input(line);
    for (std::string field; std::getline(input, field, '\t');) result.push_back(field);
    return result;
}

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: score_read_catalog INPUT_FOLDER\n";
        return 1;
    }
    const std::string prefix = std::string(argv[1]) + "/matrix.chunk0.";
    std::map<int, std::set<std::string>> old;
    std::map<int, std::vector<pgphase_collect::MoleculeAlleleSequence>> observations;
    std::ifstream alleles(prefix + "nested-alleles.tsv"), cohort(prefix + "physical-cohort.tsv");
    if (!alleles || !cohort) throw std::runtime_error("Cannot open raw catalogs/cohort in " + prefix);
    std::string line;
    std::getline(alleles, line);
    while (std::getline(alleles, line)) {
        const auto row = fields(line);
        old[std::stoi(row.at(0))].insert(row.at(2));
    }
    std::getline(cohort, line);
    while (std::getline(cohort, line)) {
        const auto row = fields(line);
        observations[std::stoi(row.at(0))].push_back({row.at(1), row.at(5)});
    }
    std::cout << "physical\tread\tsequence\tdistance\theldout_admitted\n";
    for (const auto& [physical, reads] : observations) {
        const auto full = pgphase_collect::build_read_allele_catalog(reads);
        if (full.status != pgphase_collect::ReadAlleleCatalogStatus::Complete)
            throw std::runtime_error("Full read catalog limited");
        std::vector<std::string> novel;
        for (const auto& hypothesis : full.hypotheses)
            if (!old.at(physical).count(hypothesis.sequence)) novel.push_back(hypothesis.sequence);
        for (const auto& read : reads) {
            const auto heldout = pgphase_collect::build_read_allele_catalog(reads, read.molecule);
            if (heldout.status != pgphase_collect::ReadAlleleCatalogStatus::Complete)
                throw std::runtime_error("Held-out read catalog limited");
            std::set<std::string> admitted;
            for (const auto& hypothesis : heldout.hypotheses) admitted.insert(hypothesis.sequence);
            const auto distances = pgphase_collect::score_allele_sequences(novel, read.sequence);
            for (size_t i = 0; i < novel.size(); ++i)
                std::cout << physical << '\t' << read.molecule << '\t' << novel[i] << '\t'
                          << distances[i] << '\t' << admitted.count(novel[i]) << '\n';
        }
    }
}
