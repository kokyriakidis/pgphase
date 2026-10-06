#!/usr/bin/env bash
set -euo pipefail

if [[ "${1:-}" == "--help" ]]; then
    printf '%s\n' 'Usage: replay.sh' 'Replay the terminal chr20 owner at default, MAPQ1, and permissive thresholds.'
    exit 0
fi
if [[ $# -ne 0 ]]; then
    printf '%s\n' 'Usage: replay.sh [--help]' >&2
    exit 1
fi
readonly repo_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)"
cd -- "${repo_root}"
readonly scratch='test_data/tmp_gap_fix70'
readonly contig_region='CHM13#0#chr20:66000001-66210255'
readonly -a inputs=(--ref test_data/chm13v2.0.chr20.renamed.fa
    --bam test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam)
export LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
for stage in baseline mapq1 permissive; do
    extra=()
    if [[ "${stage}" == mapq1 ]]; then extra=(--min-mapq 1); fi
    if [[ "${stage}" == permissive ]]; then extra=(--min-mapq 0 --min-depth 1 --min-alt-depth 1); fi
    output="${scratch}/${stage}/66"
    mkdir -p -- "${output}"
    PGPHASE_DEBUG_SITE=66206480 ./pgphase collect-graph-variation "${inputs[@]}" \
        --sites test_data/chr20.sites.striped.vcf.gz \
        --gaf test_data/HG002.chr20.annotated.coord.gaf.gz -t 4 -r "${contig_region}" \
        -o "${output}/candidates.tsv" --phased-vcf-out "${output}/phased.vcf" \
        --phased-bam-out "${output}/phased.bam" --filtered-sites-out "${output}/filtered.tsv" \
        "${extra[@]}" > "${output}/stdout.log" 2> "${output}/stderr.log"
done
for depth in 5 3; do
    output="${scratch}/bam_depth${depth}"
    mkdir -p -- "${output}"
    ./pgphase collect-bam-variation "${inputs[@]}" -q 1 -D "${depth}" -V 1 -t 4 \
        -r "${contig_region}" -o "${output}/candidates.tsv" \
        --phased-vcf-out "${output}/phased.vcf" \
        > "${output}/stdout.log" 2> "${output}/stderr.log"
done
