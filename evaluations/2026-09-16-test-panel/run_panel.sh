#!/usr/bin/env bash
# Run one arm over every window in the test panel.
#
# The panel is the fixed set of windows the pipeline leaves unphased and a
# competitor spans at 100% read-level accuracy, so each iteration is one command
# and the target never has to be re-derived. Pass the arm's extra flags in FLAGS;
# with FLAGS empty this is the stock pipeline.
#
#   OUT=/tmp/panel-stock ./run_panel.sh
#   OUT=/tmp/panel-retry FLAGS="--retry-unphased-with-bam" ./run_panel.sh
set -euo pipefail
readonly REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
readonly OUT="${OUT:?set OUT}"
readonly PANEL="${PANEL:-${REPO}/evaluations/2026-09-16-test-panel/panel.tsv}"
readonly FLANK="${FLANK:-50000}"
readonly THREADS="${THREADS:-8}"
readonly FLAGS="${FLAGS:-}"
mkdir -p "${OUT}"
# tr -d '\r' because csv.writer once defaulted to CRLF and a stray \r made the
# last field unparseable, which silently matched zero windows.
tail -n +2 "${PANEL}" | tr -d '\r' | while IFS=$'\t' read -r gl gr bp reads comp acc scored; do
    [[ -z "${gl}" ]] && continue
    d="${OUT}/w${gl}"; mkdir -p "${d}"
    lo=$(( gl - FLANK )); hi=$(( gr + FLANK ))
    # shellcheck disable=SC2086
    "${REPO}/pgphase" collect-graph-variation \
        --ref "${REPO}/test_data/chm13v2.0.chr20.renamed.fa" \
        --bam "${REPO}/test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam" \
        --sites "${REPO}/test_data/chr20.sites.striped.vcf.gz" \
        --gaf "${REPO}/test_data/HG002.chr20.annotated.coord.gaf.gz" \
        -r "CHM13#0#chr20:${lo}-${hi}" -t "${THREADS}" ${FLAGS} \
        -o "${d}/candidates.tsv" --phased-vcf-out "${d}/native.vcf" \
        --phased-bam-out "${d}/phased.bam" > "${d}/stdout.log" 2> "${d}/stderr.log" \
        || { echo "  FAILED ${gl}: $(tail -1 "${d}/stderr.log")"; continue; }
    echo "  w${gl} done"
done
echo "panel complete: ${OUT}"
