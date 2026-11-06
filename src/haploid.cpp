#include "haploid.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <utility>

namespace pgphase_haploid {
namespace {

using Intervals = std::vector<std::pair<hts_pos_t, hts_pos_t>>;  // PAR, 0-based half-open

std::unordered_map<std::string, Intervals>& contigs_map() {
    static std::unordered_map<std::string, Intervals> m;
    return m;
}

std::string_view short_name(std::string_view name) {
    const size_t hash = name.rfind('#');
    return hash == std::string_view::npos ? name : name.substr(hash + 1);
}

}  // namespace

void configure(const std::vector<std::string>& contigs, const std::string& par_bed) {
    auto& m = contigs_map();
    m.clear();
    for (const std::string& c : contigs)
        if (!c.empty()) m.emplace(std::string(short_name(c)), Intervals{});
    if (par_bed.empty()) return;
    std::ifstream in(par_bed);
    if (!in) throw std::runtime_error("cannot open --par-bed " + par_bed);
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#' || line.rfind("track", 0) == 0 || line.rfind("browser", 0) == 0) continue;
        std::istringstream f(line);
        std::string chrom;
        hts_pos_t beg = 0, end = 0;
        if (!(f >> chrom >> beg >> end) || beg < 0 || end < beg) throw std::runtime_error("bad --par-bed line: " + line);
        const auto it = m.find(std::string(short_name(chrom)));
        if (it != m.end()) it->second.emplace_back(beg, end);
    }
}

bool enabled() { return !contigs_map().empty(); }

bool at(std::string_view contig, hts_pos_t pos) {
    const auto& m = contigs_map();
    if (m.empty()) return false;
    const auto it = m.find(std::string(short_name(contig)));
    if (it == m.end()) return false;
    for (const auto& [beg, end] : it->second)
        if (pos > beg && pos <= end) return false;
    return true;
}

}  // namespace pgphase_haploid
