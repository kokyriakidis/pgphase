#include <exception>
#include <iostream>
#include <string>

#include "build_catalog.hpp"
#include "collect_pipeline.hpp"
#include "graph_collect.hpp"
#include "version.hpp"

namespace {

void print_main_help() {
    std::cerr << "pgphase " PGPHASE_VERSION "\n"
              << "\n"
              << "Usage: pgphase <command> [options]\n"
              << "Commands:\n"
              << "  collect-graph-variation  Phase reads from graph alignments (GAF), optionally with a BAM\n"
              << "  collect-bam-variation    Phase reads from linear alignments (BAM/CRAM) only\n"
              << "  build-snarl-catalog      Build the phasing site catalog from a GBZ graph\n";
}

} // namespace

int main(int argc, char* argv[]) {
    try {
        if (argc > 1 && (std::string(argv[1]) == "--version" || std::string(argv[1]) == "version")) {
            std::cout << "pgphase " PGPHASE_VERSION "\n";
            return 0;
        }
        if (argc > 1 && std::string(argv[1]) == "collect-bam-variation") {
            return collect_bam_variation(argc - 1, argv + 1);
        }
        if (argc > 1 && std::string(argv[1]) == "collect-graph-variation") {
            return collect_graph_variation(argc - 1, argv + 1);
        }
        if (argc > 1 && std::string(argv[1]) == "build-snarl-catalog") {
            return build_snarl_catalog(argc - 1, argv + 1);
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
