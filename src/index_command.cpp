#include "pgcatalog_command.hpp"

#include <cerrno>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>

// pgphase index: the reads' graph alignments indexed by reference position,
// as collect-graph-variation --gaf reads them. dromos (github.com/kokyriakidis/
// dromos) writes the index, placing each alignment by the graph's reference map,
// which pgcatalog ref-map makes from the graph once and which is kept for the
// next run.

namespace {

const char* const kUsage =
    "pgphase index: index the reads' graph alignments by reference position, for collect-graph-variation --gaf\n"
    "\n"
    "  pgphase index --graph GRAPH.gbz --reads READS [options]\n"
    "\n"
    "  --graph FILE              the GBZ graph the reads were aligned to\n"
    "  --reads FILE              their alignments: a GAM, or a GAF (plain, gzip or\n"
    "                            bgzip), told apart by their first bytes; - reads stdin\n"
    "  --out FILE                [READS.coord.gaf.gz, in the current folder], with\n"
    "                            FILE.tbi beside it\n"
    "  --reference-sample NAME   the sample whose paths give the coordinates [the\n"
    "                            graph's declared reference sample, if it declares\n"
    "                            exactly one]\n"
    "  --ref-map FILE            the graph's reference map (pgcatalog ref-map) [made\n"
    "                            once and kept beside the graph as\n"
    "                            GRAPH.gbz[.SAMPLE].refmap, or beside the output if\n"
    "                            the graph's folder is not writable]\n"
    "  -t, --threads N           threads [all cores]\n"
    "  --level N                 BGZF compression level, 0-12 [1]; 6 for a smaller\n"
    "                            file\n"
    "\n"
    "Each alignment is placed on the reference by the walk's nodes on the reference\n"
    "sample's paths and written as `contig  start  end  <its 12 GAF columns>`\n"
    "(1-based, inclusive; tags rq, pu and rs kept), sorted, bgzipped and\n"
    "tabix-indexed. Unmapped alignments, and those without a unique reference\n"
    "interval, are left out and counted by reason. pgcatalog makes the reference\n"
    "map and dromos the index: PGCATALOG and DROMOS name them [beside pgphase, then\n"
    "on PATH].\n";

int fail(const std::string& message) {
    std::cerr << "error: " << message << "\nsee pgphase index --help\n";
    return 2;
}

bool is_number(const std::string& s) {
    return !s.empty() && s.size() < 10 && s.find_first_not_of("0123456789") == std::string::npos;
}

std::string folder_of(const std::string& path) {
    const size_t slash = path.rfind('/');
    return slash == std::string::npos ? "." : path.substr(0, slash == 0 ? 1 : slash);
}

std::string name_of(const std::string& path) {
    const size_t slash = path.rfind('/');
    return slash == std::string::npos ? path : path.substr(slash + 1);
}

bool ends_with(const std::string& s, const std::string& end) {
    return s.size() >= end.size() && s.compare(s.size() - end.size(), end.size(), end) == 0;
}

// pgcat pgphase index's default output: the input's name less .gz, .bgz and .gaf (or
// .gam), with .coord.gaf.gz, in the current folder.
std::string default_out(const std::string& input) {
    std::string stem = name_of(input);
    for (const char* end : {".gz", ".bgz", ".gaf", ".gam"}) {
        while (ends_with(stem, end)) stem.resize(stem.size() - std::strlen(end));
    }
    return stem + ".coord.gaf.gz";
}

// Runs `args` to its end; its exit status (127 if it cannot run).
int run_and_wait(const std::vector<std::string>& args) {
    std::vector<char*> argv;
    for (const std::string& a : args) argv.push_back(const_cast<char*>(a.c_str()));
    argv.push_back(nullptr);
    const pid_t pid = fork();
    if (pid < 0) return 127;
    if (pid == 0) {
        execvp(argv[0], argv.data());
        std::cerr << "error: cannot run " << args[0] << " (" << std::strerror(errno) << ")\n";
        _exit(127);
    }
    int status = 0;
    while (waitpid(pid, &status, 0) < 0) {
        if (errno != EINTR) return 127;
    }
    return WIFEXITED(status) ? WEXITSTATUS(status) : 128 + WTERMSIG(status);
}

// The kept reference map is current when it is newer than the graph.
bool current(const std::string& map, const std::string& graph) {
    struct stat m{}, g{};
    if (stat(map.c_str(), &m) != 0 || stat(graph.c_str(), &g) != 0) return false;
    return m.st_mtime >= g.st_mtime;
}

}  // namespace

int run_index(int argc, char* argv[], const char* argv0) {
    std::string graph, reads, out, sample, ref_map, threads, level;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "-h" || a == "--help") {
            std::cout << kUsage;
            return 0;
        }
        auto value = [&](std::string& into) -> bool {
            if (i + 1 >= argc) return false;
            into = argv[++i];
            return true;
        };
        bool ok = true;
        if (a == "--graph") ok = value(graph);
        else if (a == "--reads") ok = value(reads);
        else if (a == "--gaf" || a == "--gam")
            return fail(a + " is gone: give the alignments with --reads (dromos tells a GAM from a GAF by its bytes)");
        else if (a == "--out") ok = value(out);
        else if (a == "--reference-sample") ok = value(sample);
        else if (a == "--ref-map") ok = value(ref_map);
        else if (a == "-t" || a == "--threads") ok = value(threads);
        else if (a == "--level") ok = value(level);
        else if (a == "--memory" || a == "--tmp-dir")
            return fail(a + " is gone: the index is sorted in memory (about 16 GB for a 30x human genome)");
        else if (a == "--unplaced")
            return fail("--unplaced is gone: the alignments without a reference interval are counted, not kept");
        else if (a == "--slim")
            return fail("--slim is gone: the index always holds each alignment's 12 GAF columns (and rq, pu, rs)");
        else return fail("unknown option " + a);
        if (!ok) return fail(a + " needs a value");
    }
    if (reads.empty()) return fail("give the alignments with --reads");
    if (graph.empty() && ref_map.empty()) return fail("--graph is needed (or --ref-map)");
    if (!threads.empty() && (!is_number(threads) || std::stoi(threads) < 1))
        return fail("--threads needs a whole number of at least 1, not " + threads);
    if (!level.empty() && (!is_number(level) || std::stoi(level) > 12))
        return fail("--level needs a whole number from 0 to 12, not " + level);
    if (out.empty()) out = default_out(reads == "-" ? std::string("reads") : reads);

    if (ref_map.empty()) {
        const std::string name = name_of(graph) + (sample.empty() ? "" : "." + sample) + ".refmap";
        const std::string beside_graph = folder_of(graph) + "/" + name;
        ref_map = current(beside_graph, graph) || access(folder_of(graph).c_str(), W_OK) == 0
                      ? beside_graph
                      : folder_of(out) + "/" + name;
        if (!current(ref_map, graph)) {
            std::cerr << "[pgphase] the graph's reference map, kept for the next run: " << ref_map << "\n";
            const std::string part = ref_map + ".part" + std::to_string(getpid());
            std::vector<std::string> make = {find_tool("pgcatalog", "PGCATALOG", argv0), "ref-map", "--graph", graph,
                                             "--out", part};
            if (!sample.empty()) {
                make.push_back("--reference-sample");
                make.push_back(sample);
            }
            const int status = run_and_wait(make);
            if (status != 0) {
                std::remove(part.c_str());
                std::cerr << "error: pgcatalog ref-map failed (exit " << status << "); build it with `make pgcatalog`\n"
                          << "(needs Rust), or point PGCATALOG at it\n";
                return status == 127 ? 127 : 1;
            }
            if (std::rename(part.c_str(), ref_map.c_str()) != 0) {
                std::remove(part.c_str());
                std::cerr << "error: cannot write " << ref_map << " (" << std::strerror(errno) << ")\n";
                return 1;
            }
        }
    }

    const std::string dromos = find_tool("dromos", "DROMOS", argv0);
    std::vector<std::string> args = {dromos, "index"};
    if (!threads.empty()) args.insert(args.end(), {"-t", threads});
    if (!level.empty()) args.insert(args.end(), {"-l", level});
    args.insert(args.end(), {ref_map, reads, out});
    std::vector<char*> exec_args;
    for (const std::string& a : args) exec_args.push_back(const_cast<char*>(a.c_str()));
    exec_args.push_back(nullptr);
    execvp(dromos.c_str(), exec_args.data());
    std::cerr << "error: cannot run " << dromos << " (" << std::strerror(errno) << ").\n"
              << "pgphase index runs dromos index: build it with `make dromos` (needs pixi),\n"
              << "or point DROMOS at the binary.\n";
    return 127;
}
