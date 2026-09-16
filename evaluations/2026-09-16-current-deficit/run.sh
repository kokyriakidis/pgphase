#!/usr/bin/env bash
# Whole-chr20 run of the CURRENT full pipeline (graph sites + BAM channel, -q 1,
# gap recovery on) so the deficit can be measured against what the pipeline
# actually leaves open. The earlier deficit list was built from the graph-only
# pass-1 gap inventory, which overstates the remaining work: the hybrid path
# already spans many of those gaps.
set -euo pipefail
readonly REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
readonly OUT="${OUT:-/tmp/hybrid-chr20-current}"
readonly THREADS="${THREADS:-16}"
mkdir -p "${OUT}"
"${REPO}/pgphase" collect-hybrid-variation \
    --ref "${REPO}/test_data/chm13v2.0.chr20.renamed.fa" \
    --bam "${REPO}/test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam" \
    --graph-sites "${REPO}/test_data/chr20.sites.striped.vcf.gz" \
    --gaf "${REPO}/test_data/HG002.chr20.annotated.coord.gaf.gz" \
    -r 'CHM13#0#chr20' -t "${THREADS}" -q 1 \
    --link-by-alleles --block-link-window 8 --min-read-margin 2 \
    --recover-gaps --gap-recovery-report "${OUT}/tiers.tsv" \
    -o "${OUT}/candidates.tsv" --phased-vcf-out "${OUT}/native.vcf" \
    -b "${OUT}/phased.bam" > "${OUT}/stdout.log" 2> "${OUT}/stderr.log"
