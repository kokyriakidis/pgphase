#!/usr/bin/env bash
# Phase each graph gap independently with the BAM channel, over the gap plus a
# flank wide enough to share reads with both adjacent graph blocks. The flank is
# what makes the transfer possible: reads that carry both a BAM haplotype and a
# graph phase set define the orientation mapping between the two gauges.
set -euo pipefail
readonly REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
readonly OUT="${OUT:-/tmp/bam-phase-transfer}"
readonly WINDOWS="${WINDOWS:-/tmp/windows.tsv}"
readonly FLANK="${FLANK:-50000}"
readonly MIN_MAPQ="${MIN_MAPQ:-1}"
readonly THREADS="${THREADS:-8}"
mkdir -p "${OUT}"
tail -n +2 "${WINDOWS}" | tr -d '\r' | while IFS=$'\t' read -r name beg end kind; do
    [[ "${kind}" == gap ]] || continue
    lo=$(( beg > FLANK ? beg - FLANK : 1 )); hi=$(( end + FLANK ))
    d="${OUT}/${name}"; mkdir -p "${d}"
    echo "[$(date +%H:%M:%S)] ${name}  CHM13#0#chr20:${lo}-${hi}"
    "${REPO}/pgphase" collect-bam-variation \
        --ref "${REPO}/test_data/chm13v2.0.chr20.renamed.fa" \
        --bam "${REPO}/test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam" \
        -r "CHM13#0#chr20:${lo}-${hi}" -t "${THREADS}" -q "${MIN_MAPQ}" \
        -o "${d}/candidates.tsv" --phased-vcf-out "${d}/native.vcf" \
        -b "${d}/phased.bam" > "${d}/stdout.log" 2> "${d}/stderr.log" || {
            echo "  FAILED: $(tail -2 "${d}/stderr.log")"; continue; }
    printf '%s\t%s\t%s\t%s\n' "${name}" "${lo}" "${hi}" "${FLANK}" > "${d}/window.txt"
done
echo "[$(date +%H:%M:%S)] done"
