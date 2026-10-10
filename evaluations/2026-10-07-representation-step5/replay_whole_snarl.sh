#!/usr/bin/env bash
set -euo pipefail

if [[ "${1:-}" == "--help" || "$#" -lt 1 || "$#" -gt 3 ]]; then
    echo "Usage: replay_matrix.sh OUTPUT [BINARY] [REGION]"
    echo "Replay an owning chunk with whole multi-allelic snarls and complete-contrast diagnostics."
    exit 0
fi
readonly repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${repo_root}"
readonly output_dir="${1}"
readonly binary="${2:-./pgphase}"
readonly region="${3:-CHM13#0#chr20:4000001-5000000}"
mkdir -p "${output_dir}"
"${binary}" collect-graph-variation \
    --ref test_data/chm13v2.0.chr20.renamed.fa \
    --bam test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam \
    --sites test_data/chr20.sites.striped.vcf.gz \
    --gaf test_data/HG002.chr20.annotated.coord.gaf.gz \
    --snarl-allele-phasing --snarl-keep-whole \
    -t 1 -r "${region}" \
    -o "${output_dir}/candidates.tsv" \
    --phased-vcf-out "${output_dir}/phased.vcf" \
    --phased-bam-out "${output_dir}/phased.bam" \
    --phase-matrix-dump "${output_dir}/matrix" \
    > "${output_dir}/stdout.log" 2> "${output_dir}/stderr.log"
