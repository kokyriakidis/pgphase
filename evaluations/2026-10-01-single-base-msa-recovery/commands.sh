#!/usr/bin/env bash
set -euo pipefail

usage() {
    echo "Usage: commands.sh BEFORE_BINARY AFTER_BINARY OUTPUT_DIRECTORY"
}

die() {
    echo "error: $*" >&2
    exit 1
}

if [[ "${1:-}" == "--help" ]]; then
    usage
    exit 0
fi
[[ $# -eq 3 ]] || { usage; die "expected two binaries and an output directory"; }
readonly before_binary="$1"
readonly after_binary="$2"
readonly output_directory="$3"
readonly repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${repo_root}"

run_chr20() {
    local -r binary="$1"
    local -r directory="$2"
    [[ -x "${binary}" ]] || die "binary is not executable: ${binary}"
    mkdir -p "${directory}"
    /usr/bin/time -v -o "${directory}/resources.txt" \
        "${binary}" collect-graph-variation \
        --ref test_data/chm13v2.0.chr20.renamed.fa \
        --bam test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam \
        --sites test_data/chr20.sites.striped.vcf.gz \
        --gaf test_data/HG002.chr20.annotated.coord.gaf.gz \
        -t 8 -r 'CHM13#0#chr20' \
        -o "${directory}/candidates.tsv" \
        --phased-vcf-out "${directory}/phased.vcf" \
        --phased-bam-out "${directory}/phased.bam" \
        > "${directory}/stdout.log" 2> "${directory}/stderr.log"
}

run_chr20 "${before_binary}" "${output_directory}/before"
run_chr20 "${after_binary}" "${output_directory}/after"
