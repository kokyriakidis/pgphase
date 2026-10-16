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

// Run one shell pipeline; pipefail so a vg failure part way (e.g. out of
// memory) fails the step instead of hiding behind the last command's status.
int run_step(const std::string& pipeline) {
    std::cerr << "[build-snarl-catalog] " << pipeline << "\n";
    return std::system(("bash -o pipefail -c " + bc_shell_quote(pipeline)).c_str());
}

std::string parent_dir(const std::string& path) {
    const size_t slash = path.rfind('/');
    if (slash == std::string::npos) return ".";
    return slash == 0 ? "/" : path.substr(0, slash);
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
        << "  With --contig, the contig's own snarls are computed first, from its\n"
        << "  connected component alone (vg chunk -C, vg snarls), so vg deconstruct\n"
        << "  visits only that contig's snarls. Given whole-genome snarls, every\n"
        << "  per-contig run would enumerate the haplotypes of every snarl in the\n"
        << "  genome (chr20: ~25 min this way against an estimated ~14 h).\n"
        << "\n"
        << "Options:\n"
        << "  --ref-sample STR    Reference sample name, e.g. CHM13 [required]\n"
        << "  --contig STR        Process one contig, e.g. chr20 (recommended: per-chr)\n"
        << "  -t, --threads N     Threads for vg deconstruct [4]\n"
        << "  --vg-bin PATH       Path to vg binary [vg]\n"
        << "  --snarls FILE       Snarls file from 'vg snarls'; with --contig it should\n"
        << "                      hold that contig's snarls only [computed per contig]\n"
        << "  --tmp-dir DIR       Where the per-contig graph and snarls are built\n"
        << "                      [the output file's directory]\n"
        << "  -o, --output FILE   Output bgzipped VCF file (tabix index created automatically) [required]\n"
        << "  -h, --help          Show this help\n"
        << "\n"
        << "Example:\n"
        << "  for chr in chr{1..22} chrX; do\n"
        << "    pgphase build-snarl-catalog --ref-sample CHM13 --contig $chr \\\n"
        << "        -t 16 -o ${chr}.sites.vcf.gz graph.gbz\n"
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
    std::string tmp_dir;

    static const option long_opts[] = {
        {"ref-sample", required_argument, nullptr, 'R'},
        {"contig",     required_argument, nullptr, 'N'},
        {"threads",    required_argument, nullptr, 't'},
        {"vg-bin",     required_argument, nullptr, 'V'},
        {"snarls",     required_argument, nullptr, 'S'},
        {"output",     required_argument, nullptr, 'o'},
        {"tmp-dir",    required_argument, nullptr, 'T'},
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
        case 'T': tmp_dir     = optarg; break;
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
    const std::string vg = bc_shell_quote(vg_bin);

    // Per-contig snarls: the contig's connected component alone, so the
    // snarls (and every deconstruct traversal) stay within the contig.
    std::string work_dir;
    if (!contig.empty() && snarls_file.empty()) {
        std::string tmpl = (tmp_dir.empty() ? parent_dir(output_file) : tmp_dir) + "/pgphase_catalog_XXXXXX";
        if (mkdtemp(tmpl.data()) == nullptr) {
            std::cerr << "error: cannot create a temporary directory like " << tmpl << "\n";
            return 1;
        }
        work_dir = tmpl;
        const std::string pg = work_dir + "/" + contig + ".pg";
        snarls_file = work_dir + "/" + contig + ".snarls";
        if (run_step(vg + " chunk -x " + bc_shell_quote(gbz_file) + " -C -p " +
                     bc_shell_quote(ref_sample + "#0#" + contig) + " -O pg -t " +
                     std::to_string(threads) + " > " + bc_shell_quote(pg)) != 0 ||
            run_step(vg + " snarls -t " + std::to_string(threads) + " " + bc_shell_quote(pg) +
                     " > " + bc_shell_quote(snarls_file)) != 0) {
            std::cerr << "[build-snarl-catalog] error: computing the snarls of " << contig
                      << " failed (intermediate files kept in " << work_dir << ")\n";
            return 1;
        }
        std::remove(pg.c_str());
    }

    // Build the vg deconstruct command.
    std::ostringstream deconstruct_cmd;
    deconstruct_cmd << vg << " deconstruct";
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
    int ret = run_step(deconstruct_cmd.str() + " | " + strip_samples + " | bgzip -@ " +
                       std::to_string(threads) + " -c > " + bc_shell_quote(output_file));
    if (ret == 0) ret = run_step("tabix -f -p vcf " + bc_shell_quote(output_file));
    if (!work_dir.empty()) {
        std::remove(snarls_file.c_str());
        rmdir(work_dir.c_str());
    }

    if (ret != 0)
        std::cerr << "[build-snarl-catalog] error: command failed; " << output_file
                  << " is incomplete\n";
    return ret == 0 ? 0 : 1;
}
