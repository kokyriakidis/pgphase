#include <exception>
#include <iostream>
#include <string>

#include "collect_pipeline.hpp"
#include "graph_collect.hpp"
#include "pgcat_command.hpp"
#include "plegma_msa.hpp"
#include "version.hpp"

namespace {

void print_main_help() {
    std::cerr << "pgphase " PGPHASE_VERSION "\n"
              << "\n"
              << "Usage: pgphase <command> [options]\n"
              << "Commands:\n"
              << "  collect-graph-variation  Phase reads from graph alignments (GAF), optionally with a BAM\n"
              << "  collect-bam-variation    Phase reads from linear alignments (BAM/CRAM) only\n"
              << "  build-snarl-catalog      Build the site catalogs from a GBZ graph (runs pgcat build)\n"
              << "  index-gaf                Index graph alignments (GAF) by reference position (runs pgcat index-gaf)\n";
}

} // namespace

int main(int argc, char* argv[]) {
    try {
        if (argc > 1 && (std::string(argv[1]) == "--version" || std::string(argv[1]) == "version")) {
            std::cout << "pgphase " PGPHASE_VERSION "\n"
                      << "noisy-region MSA: " << pgphase_collect::msa_backends() << "\n";
            return 0;
        }
        if (argc > 1 && std::string(argv[1]) == "collect-bam-variation") {
            return collect_bam_variation(argc - 1, argv + 1);
        }
        if (argc > 1 && std::string(argv[1]) == "collect-graph-variation") {
            return collect_graph_variation(argc - 1, argv + 1);
        }
        if (argc > 1 && std::string(argv[1]) == "build-snarl-catalog") {
            return run_pgcat("build", argc - 1, argv + 1, argv[0]);
        }
        if (argc > 1 && std::string(argv[1]) == "index-gaf") {
            return run_pgcat("index-gaf", argc - 1, argv + 1, argv[0]);
        }
    } catch (const std::exception& e) {
        // Turn any uncaught error (bad CLI value, missing file, etc.) into a
        // clean diagnostic and non-zero exit instead of an abort/core dump.
        std::cerr << "error: " << e.what() << '\n';
        return 1;
    }

    print_main_help();
    return 1;
}
