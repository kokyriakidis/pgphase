#!/usr/bin/env bash
set -euo pipefail

readonly REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
readonly OUT_ROOT="${OUT_ROOT:-/tmp/pgphase_native_bam_graph_lock}"
readonly TRUTH_ROOT="${TRUTH_ROOT:-${HOME}/Downloads/pgphase-eval-data/truth}"
readonly SITES="${SITES:-${HOME}/Downloads/chr20.sites.vcf.gz}"
readonly THREADS="${THREADS:-20}"
readonly REF="test_data/chm13v2.0.chr20.renamed.fa"
readonly GAF="test_data/HG002.chr20.annotated.coord.gaf.gz"
readonly BAM="test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam"
readonly PRIVATE="evaluations/2026-09-12-native-bam-graph-lock/chr20/inputs/private_gap_gq10.vcf"
readonly TRUTH="${TRUTH_ROOT}/chr20/diplinator_merged.bam"
readonly OUT="${OUT_ROOT}/chr20"

cd "${REPO_ROOT}"
mkdir -p "${OUT}/graph/eval" "${OUT}/joint_gq10_ps0/eval" \
    "${OUT}/joint_gq10_ps50/eval" "${OUT}/graph_lock_ps0/eval" \
    "${OUT}/graph_lock_final50/eval"

./pgphase collect-graph-variation \
    --ref "${REF}" --sites "${SITES}" --gaf "${GAF}" \
    --min-read-margin 2 --anchor-af-margin 0.12 \
    -r chr20 -t "${THREADS}" \
    -o "${OUT}/graph/candidates.tsv" \
    --phased-vcf-out "${OUT}/graph/phased.vcf" \
    --phased-bam-out "${OUT}/graph/phased.bam"

./pgphase collect-hybrid-variation \
    --ref "${REF}" --bam "${BAM}" --graph-sites "${SITES}" --gaf "${GAF}" \
    --private-sites "${PRIVATE}" --min-read-margin 2 --min-phase-set-reads 0 \
    --stitch-min-margin 0 --stitch-rule 0 -r chr20 -t "${THREADS}" \
    -o "${OUT}/joint_gq10_ps0/candidates.tsv" \
    --phased-vcf-out "${OUT}/joint_gq10_ps0/phased.vcf" \
    -b "${OUT}/joint_gq10_ps0/phased.bam"

./pgphase collect-hybrid-variation \
    --ref "${REF}" --bam "${BAM}" --graph-sites "${SITES}" --gaf "${GAF}" \
    --private-sites "${PRIVATE}" --min-read-margin 2 --min-phase-set-reads 50 \
    --stitch-min-margin 0 --stitch-rule 0 -r chr20 -t "${THREADS}" \
    -o "${OUT}/joint_gq10_ps50/candidates.tsv" \
    --phased-vcf-out "${OUT}/joint_gq10_ps50/phased.vcf" \
    -b "${OUT}/joint_gq10_ps50/phased.bam"

python3 scripts/merge_graph_hybrid_tags.py \
    --graph-bam "${OUT}/graph/phased.bam" \
    --hybrid-bam "${OUT}/joint_gq10_ps0/phased.bam" \
    --output "${OUT}/graph_lock_ps0/phased.bam" \
    --min-shared-reads 10 --min-vote-margin 5 --min-purity 0.90 \
    --require-both-haplotypes --threads 8

python3 scripts/merge_graph_hybrid_tags.py \
    --graph-bam "${OUT}/graph/phased.bam" \
    --hybrid-bam "${OUT}/joint_gq10_ps0/phased.bam" \
    --output "${OUT}/graph_lock_final50/phased.bam" \
    --min-shared-reads 10 --min-vote-margin 5 --min-purity 0.90 \
    --require-both-haplotypes --min-output-phase-set-reads 50 --threads 8

for config in graph joint_gq10_ps0 joint_gq10_ps50 graph_lock_ps0 \
        graph_lock_final50; do
    python3 scripts/evaluate_phase_accuracy.py \
        "${OUT}/${config}/phased.bam" "${TRUTH}" 0 0 5 "" \
        "${OUT}/${config}/eval" samtools "" "" "" ""
done
