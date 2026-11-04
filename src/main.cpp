#include <exception>
#include <iostream>
#include <string>

#include "collect_pipeline.hpp"
#include "graph_collect.hpp"
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
              << "\n"
              << "The site catalog and the indexed GAF are made by dromos (github.com/kokyriakidis/dromos):\n"
              << "  dromos catalog --graph graph.gbz --reference-sample CHM13\n"
              << "  dromos index --graph graph.gbz --reference-sample CHM13 reads.gam reads.coord.gaf.gz\n";
}

// Commands pgphase no longer has: the dromos command that replaces each.
int removed_command(const std::string& command, const char* dromos_command) {
    std::cerr << "error: pgphase " << command << " is gone: dromos makes pgphase's inputs"
              << " (github.com/kokyriakidis/dromos).\n"
              << "  " << dromos_command << "\n";
    return 1;
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
            return removed_command(argv[1], "dromos catalog --graph graph.gbz --reference-sample CHM13");
        }
        if (argc > 1 && std::string(argv[1]) == "index") {
            return removed_command(
                argv[1], "dromos index --graph graph.gbz --reference-sample CHM13 reads.gam reads.coord.gaf.gz");
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
