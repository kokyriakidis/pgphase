#include "build_catalog.hpp"

#include "arg_parse.hpp"

#include <cstdio>
#include <cstdlib>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/wait.h>
#include <unistd.h>

namespace {

std::string bc_shell_quote(const std::string& val) {
    std::string out = "'";
    for (char c : val) {
        if (c == '\'') out += "'\\''";
        else out.push_back(c);
    }
    out.push_back('\'');
    return out;
}

void print_help() {
    std::cerr
        << "Usage: pgphase build-snarl-catalog [options] <graph.gbz>\n"
        << "\n"
        << "Run 'vg deconstruct -a' on a GBZ pangenome graph and write a bgzipped,\n"
        << "tabix-indexed sites-only VCF (no sample genotypes) for use as a\n"
        << "phasing site catalog with 'collect-graph-variation'.\n"
        << "\n"
        << "  <graph.gbz> is the pangenome GBZ; the catalog covers the reference path\n"
        << "  <ref-sample>#0#<contig> (all <ref-sample>#0# paths without --contig).\n"
        << "\n"
        << "Options:\n"
        << "  --ref-sample STR    Reference sample name, e.g. CHM13 [required]\n"
        << "  --contig STR        Process one contig, e.g. chr20 (recommended: per-chr)\n"
        << "  -t, --threads N     Threads for vg deconstruct [4]\n"
        << "  --vg-bin PATH       Path to vg binary [vg]\n"
        << "  --snarls FILE       Pre-built snarls file from 'vg snarls' (skips recomputation)\n"
        << "  -o, --output FILE   Output bgzipped VCF file (tabix index created automatically) [required]\n"
        << "  -h, --help          Show this help\n"
        << "\n"
        << "Example:\n"
        << "  vg snarls -t 16 graph.gbz > graph.snarls      # once per graph\n"
        << "  for chr in chr{1..22} chrX; do\n"
        << "    pgphase build-snarl-catalog --ref-sample CHM13 --contig $chr \\\n"
        << "        --snarls graph.snarls -t 16 -o ${chr}.sites.vcf.gz graph.gbz\n"
        << "  done\n"
        << "\n"
        << "  # Run phasing\n"
        << "  pgphase collect-graph-variation \\\n"
        << "      --ref ref.fa \\\n"
        << "      --sites chr20.sites.vcf.gz \\\n"
        << "      --gaf reads.coord.gaf.gz \\\n"
        << "      --bam reads.bam \\\n"
        << "      --phased-vcf-out phased.vcf \\\n"
        << "      --phased-bam-out phased.bam \\\n"
        << "      -t 8\n";
}

} // namespace

int build_snarl_catalog(int argc, char* argv[]) {
    std::string ref_sample;
    std::string contig;
    int threads = 4;
    std::string vg_bin = "vg";
    std::string snarls_file;
    std::string output_file;

    static const option long_opts[] = {
        {"ref-sample", required_argument, nullptr, 'R'},
        {"contig",     required_argument, nullptr, 'N'},
        {"threads",    required_argument, nullptr, 't'},
        {"vg-bin",     required_argument, nullptr, 'V'},
        {"snarls",     required_argument, nullptr, 'S'},
        {"output",     required_argument, nullptr, 'o'},
        {"help",       no_argument,       nullptr, 'h'},
        {nullptr, 0, nullptr, 0}
    };

    int ch;
    while ((ch = getopt_long(argc, argv, "t:o:h", long_opts, nullptr)) != -1) {
        switch (ch) {
        case 'R': ref_sample  = optarg; break;
        case 'N': contig      = optarg; break;
        case 't': threads     = pgphase_collect::parse_int_arg(optarg, "--threads"); break;
        case 'V': vg_bin      = optarg; break;
        case 'S': snarls_file = optarg; break;
        case 'o': output_file = optarg; break;
        case 'h': print_help(); return 0;
        default:  print_help(); return 1;
        }
    }

    if (optind >= argc) {
        std::cerr << "error: missing <graph.gbz>\n";
        print_help();
        return 1;
    }
    if (ref_sample.empty()) {
        std::cerr << "error: --ref-sample is required (e.g. --ref-sample CHM13)\n";
        return 1;
    }
    if (output_file.empty()) {
        std::cerr << "error: -o/--output is required\n";
        print_help();
        return 1;
    }

    const std::string gbz_file = argv[optind];

    // Build the vg deconstruct command.
    std::ostringstream deconstruct_cmd;
    deconstruct_cmd << bc_shell_quote(vg_bin) << " deconstruct";
    if (!contig.empty())
        deconstruct_cmd << " -p " << bc_shell_quote(ref_sample + "#0#" + contig);
    else
        deconstruct_cmd << " -P " << bc_shell_quote(ref_sample + "#0#");
    deconstruct_cmd << " -a -t " << threads;
    if (!snarls_file.empty())
        deconstruct_cmd << " -r " << bc_shell_quote(snarls_file);
    deconstruct_cmd << " " << bc_shell_quote(gbz_file);

    // Keep CHROM..INFO (fields 1-8), the only columns pgphase reads; this also
    // drops FORMAT and the sample names from the #CHROM line. cut and a
    // multithreaded bgzip keep pace with vg (awk ran at ~11 MB/s, bgzip at one
    // core), and give byte-identical output.
    const std::string strip_samples = "cut -f1-8";
    const std::string pipeline =
        deconstruct_cmd.str() + " | " + strip_samples + " | bgzip -@ " + std::to_string(threads) + " -c > " + bc_shell_quote(output_file);
    std::cerr << "[build-snarl-catalog] " << pipeline << "\n";
    // pipefail: a vg failure part way (e.g. out of memory) must fail the
    // command, not leave a truncated catalog behind bgzip's exit status.
    const std::string full_cmd = "bash -o pipefail -c " + bc_shell_quote(pipeline);
    int ret = std::system(full_cmd.c_str());

    if (ret == 0) {
        const std::string tabix_cmd = "tabix -f -p vcf " + bc_shell_quote(output_file);
        std::cerr << "[build-snarl-catalog] " << tabix_cmd << "\n";
        ret = std::system(tabix_cmd.c_str());
    }

    if (ret != 0)
        std::cerr << "[build-snarl-catalog] error: command failed; " << output_file
                  << " is incomplete\n";
    return ret == 0 ? 0 : 1;
}
