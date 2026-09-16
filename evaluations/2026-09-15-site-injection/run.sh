#!/usr/bin/env bash
# Does injecting BAM-derived sites into the graph's site set let one solve span
# the gap? Two arms per window: the site union alone (hybrid without gap
# recovery) and the union plus the gap-recovery tiers, against the graph-only
# baseline that left the gap open.
set -euo pipefail
readonly REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
readonly OUT="${OUT:-/tmp/site-injection}"
readonly WINDOWS="${WINDOWS:-/tmp/windows.tsv}"
readonly FLANK="${FLANK:-50000}"
readonly MIN_MAPQ="${MIN_MAPQ:-1}"
readonly THREADS="${THREADS:-8}"
mkdir -p "${OUT}"
tail -n +2 "${WINDOWS}" | tr -d '\r' | while IFS=$'\t' read -r name beg end kind; do
    [[ "${kind}" == gap ]] || continue
    lo=$(( beg > FLANK ? beg - FLANK : 1 )); hi=$(( end + FLANK ))
    for arm in union union_recover; do
        d="${OUT}/${name}/${arm}"; mkdir -p "${d}"
        extra=()
        [[ "${arm}" == union_recover ]] && extra=(--recover-gaps --gap-recovery-report "${d}/tiers.tsv")
        echo "[$(date +%H:%M:%S)] ${name} ${arm}"
        "${REPO}/pgphase" collect-hybrid-variation \
            --ref "${REPO}/test_data/chm13v2.0.chr20.renamed.fa" \
            --bam "${REPO}/test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam" \
            --graph-sites "${REPO}/test_data/chr20.sites.striped.vcf.gz" \
            --gaf "${REPO}/test_data/HG002.chr20.annotated.coord.gaf.gz" \
            -r "CHM13#0#chr20:${lo}-${hi}" -t "${THREADS}" -q "${MIN_MAPQ}" \
            --link-by-alleles --block-link-window 8 --min-read-margin 2 \
            "${extra[@]}" \
            -o "${d}/candidates.tsv" --phased-vcf-out "${d}/native.vcf" \
            -b "${d}/phased.bam" > "${d}/stdout.log" 2> "${d}/stderr.log" || {
                echo "  FAILED: $(tail -2 "${d}/stderr.log")"; continue; }
    done
done
echo "[$(date +%H:%M:%S)] done"
