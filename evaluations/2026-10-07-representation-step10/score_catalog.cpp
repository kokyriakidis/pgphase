// Evaluation-only residual lower bounds, with no phase/truth inputs or fitting.
#include "allele_context.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

static std::vector<std::string> fields(const std::string& line) {
    std::vector<std::string> result;
    std::istringstream input(line);
    for (std::string field; std::getline(input, field, '\t');) result.push_back(field);
    return result;
}

static std::map<int, std::vector<std::string>> alleles(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Cannot read " + path);
    std::map<int, std::vector<std::string>> result;
    std::string line;
    std::getline(input, line);
    while (std::getline(input, line)) {
        const auto row = fields(line);
        auto& sequences = result[std::stoi(row.at(0))];
        if (std::stoul(row.at(1)) != sequences.size())
            throw std::runtime_error("Unordered alleles in " + path);
        sequences.push_back(row.at(2));
    }
    return result;
}

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: score_catalog INPUT_FOLDER\n";
        return 1;
    }
    const std::string prefix = std::string(argv[1]) + "/matrix.chunk0.";
    const auto original = alleles(prefix + "physical-alleles.tsv");
    const auto composed = alleles(prefix + "composed-alleles.tsv");
    std::ifstream cohort(prefix + "physical-cohort.tsv");
    if (!cohort) throw std::runtime_error("Cannot read " + prefix + "physical-cohort.tsv");
    std::cout << "physical\tread\told_costs\told_min\tcomposed_min\n";
    std::string line;
    std::getline(cohort, line);
    while (std::getline(cohort, line)) {
        const auto row = fields(line);
        const int physical = std::stoi(row.at(0));
        const auto old_costs = pgphase_collect::score_allele_sequences(original.at(physical), row.at(5));
        const auto new_costs = pgphase_collect::score_allele_sequences(composed.at(physical), row.at(5));
        std::cout << physical << '\t' << row.at(1) << '\t';
        for (size_t i = 0; i < old_costs.size(); ++i) {
            if (i) std::cout << ',';
            std::cout << old_costs[i];
        }
        std::cout << '\t' << *std::min_element(old_costs.begin(), old_costs.end())
                  << '\t' << *std::min_element(new_costs.begin(), new_costs.end()) << '\n';
    }
}
