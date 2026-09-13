#!/usr/bin/env bash
set -euo pipefail

readonly REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
readonly OUT_ROOT="${OUT_ROOT:-/tmp/pgphase_native_bam_graph_lock}"
readonly TRUTH_ROOT="${TRUTH_ROOT:-${HOME}/Downloads/pgphase-eval-data/truth}"
readonly THREADS="${THREADS:-20}"
readonly CHROM="chr18"
readonly REF="test_data/chm13v2.0.chr18.renamed.fa"
readonly SITES="test_data/chr18.sites.vcf.gz"
readonly GAF="test_data/HG002.chr18.annotated.coord.gaf.gz"
readonly BAM="test_data/HG002_chr18_hifi_mapped_to_CHM13_chr18_annotated.bam"
readonly TRUTH="${TRUTH_ROOT}/chr18/diplinator_merged.bam"
readonly OUT="${OUT_ROOT}/chr18"

cd "${REPO_ROOT}"
mkdir -p "${OUT}/private_gq10/eval" "${OUT}/bridge_gq10/eval"

# Negative control 1: every exact-private native GQ10 site in a graph gap.
python3 scripts/extract_private_gap_sites.py \
    --graph-phased-vcf "${OUT}/graph/phased.vcf.gz" \
    --graph-sites "${SITES}" \
    --linear-vcf "${OUT}/bam/phased.vcf.gz" \
    --contig "${CHROM}" \
    --output "${OUT}/private_gap_gq10.vcf" \
    --gaps-bed "${OUT}/graph_gaps_all_gq10.bed" \
    --min-gq 10

./pgphase collect-hybrid-variation \
    --ref "${REF}" --bam "${BAM}" --graph-sites "${SITES}" --gaf "${GAF}" \
    --private-sites "${OUT}/private_gap_gq10.vcf" \
    --min-read-margin 2 --min-phase-set-reads 0 \
    --stitch-min-margin 0 --stitch-rule 0 \
    -r "${CHROM}" -t "${THREADS}" \
    -o "${OUT}/private_gq10/candidates.tsv" \
    --phased-vcf-out "${OUT}/private_gq10/phased.vcf" \
    -b "${OUT}/private_gq10/phased.bam"

python3 scripts/evaluate_phase_accuracy.py \
    "${OUT}/private_gq10/phased.bam" "${TRUTH}" 0 0 5 "" \
    "${OUT}/private_gq10/eval" samtools "" "" "" ""

# Negative control 2: coordinate-spanning read-overlap paths only.
python3 scripts/extract_private_gap_sites.py \
    --graph-phased-vcf "${OUT}/graph/phased.vcf.gz" \
    --graph-sites "${SITES}" \
    --linear-vcf "${OUT}/bam/phased.vcf.gz" \
    --contig "${CHROM}" \
    --output "${OUT}/private_bridges_gq10.vcf" \
    --gaps-bed "${OUT}/graph_gaps_bridge_gq10.bed" \
    --min-gq 10 --bam "${BAM}" --min-bridge-reads 2 --min-mapq 20

./pgphase collect-hybrid-variation \
    --ref "${REF}" --bam "${BAM}" --graph-sites "${SITES}" --gaf "${GAF}" \
    --private-sites "${OUT}/private_bridges_gq10.vcf" \
    --min-read-margin 2 --min-phase-set-reads 0 \
    --stitch-min-margin 0 --stitch-rule 0 \
    -r "${CHROM}" -t "${THREADS}" \
    -o "${OUT}/bridge_gq10/candidates.tsv" \
    --phased-vcf-out "${OUT}/bridge_gq10/phased.vcf" \
    -b "${OUT}/bridge_gq10/phased.bam"

python3 scripts/evaluate_phase_accuracy.py \
    "${OUT}/bridge_gq10/phased.bam" "${TRUTH}" 0 0 5 "" \
    "${OUT}/bridge_gq10/eval" samtools "" "" "" ""
