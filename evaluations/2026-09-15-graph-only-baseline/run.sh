#!/usr/bin/env bash
# Whole-chr20 graph-only phasing baseline: catalog sites and graph alignments
# only, no BAM channel. Produces the phase blocks that a gap subprocess would
# consume, plus read-level accuracy against the diplinator truth.
set -euo pipefail

readonly REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
readonly OUT="${OUT:-/tmp/graph-only-chr20}"
readonly DATA_ROOT="${DATA_ROOT:-${HOME}/Downloads/pgphase-eval-data}"
readonly PGPHASE="${PGPHASE:-${REPO}/pgphase}"
readonly THREADS="${THREADS:-8}"
# -q 1 is the adopted pass-1 floor: it leaves the MAPQ>=30 population's accuracy
# unchanged (0.856% -> 0.854%) while adding 5,041 clean het SNPs and 10% N50.
readonly MIN_MAPQ="${MIN_MAPQ:-1}"
readonly REGION='CHM13#0#chr20'
readonly TRUTH="${DATA_ROOT}/truth/chr20/diplinator_merged.bam"

mkdir -p "${OUT}"
"${PGPHASE}" collect-graph-variation \
    --ref "${REPO}/test_data/chm13v2.0.chr20.renamed.fa" \
    --sites "${REPO}/test_data/chr20.sites.striped.vcf.gz" \
    --gaf "${REPO}/test_data/HG002.chr20.annotated.coord.gaf.gz" \
    -r "${REGION}" -t "${THREADS}" -q "${MIN_MAPQ}" \
    -o "${OUT}/variants.tsv" --phased-vcf-out "${OUT}/native.vcf" \
    --phased-bam-out "${OUT}/phased.bam" \
    --phase-sites-out "${OUT}/phase_sites.tsv" \
    --filtered-sites-out "${OUT}/filtered_sites.tsv" \
    --phase-reads-out "${OUT}/phase_reads.tsv" \
    > "${OUT}/stdout.log" 2> "${OUT}/stderr.log"

mkdir -p "${OUT}/eval"
python3 "${REPO}/scripts/evaluate_phase_accuracy.py" "${OUT}/phased.bam" "${TRUTH}" \
    0 0 5 '' "${OUT}/eval" samtools '' '' '' '' > "${OUT}/eval.log" 2>&1
