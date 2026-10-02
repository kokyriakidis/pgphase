#!/usr/bin/env bash
# Repeat the matched chr20 comparison with binaries built before and after the fix.
set -euo pipefail
usage() {
    echo "Usage: ${0} BEFORE_BINARY AFTER_BINARY OUTPUT_DIRECTORY"
    echo "Run both graph+recovery binaries on the committed chr20 inputs."
}
die() { echo "${*}" >&2; exit 1; }
if [[ "${1:-}" == "--help" ]]; then usage; exit 0; fi
[[ "${#}" -eq 3 ]] || { usage >&2; exit 1; }
readonly before_binary="$(realpath "${1}")"
readonly after_binary="$(realpath "${2}")"
readonly output_directory="$(realpath -m "${3}")"
readonly repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
[[ -x "${before_binary}" ]] || die "Missing executable: ${before_binary}"
[[ -x "${after_binary}" ]] || die "Missing executable: ${after_binary}"
cd "${repo}"
run_one() {
    local label="${1}" binary="${2}"
    local directory="${output_directory}/${label}"
    mkdir -p "${directory}"
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
run_one before "${before_binary}"
run_one after "${after_binary}"
