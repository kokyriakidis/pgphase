#!/usr/bin/env bash
set -euo pipefail

if [[ "${1:-}" == "--help" ]]; then
    echo "Usage: replay.sh [output-directory]"
    echo "Run the full chr20 graph pipeline with four threads."
    exit 0
fi
readonly repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${repo_root}"
readonly output_dir="${1:-test_data/tmp_gap_fix74/frozen_final/0}"
mkdir -p "${output_dir}"
LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu ./pgphase collect-graph-variation \
    --ref test_data/chm13v2.0.chr20.renamed.fa \
    --bam test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam \
    --sites test_data/chr20.sites.striped.vcf.gz \
    --gaf test_data/HG002.chr20.annotated.coord.gaf.gz \
    -t 4 -r CHM13#0#chr20:1-66210255 \
    -o "${output_dir}/candidates.tsv" \
    --phased-vcf-out "${output_dir}/phased.vcf" \
    --phased-bam-out "${output_dir}/phased.bam" \
    > "${output_dir}/stdout.log" 2> "${output_dir}/stderr.log"
