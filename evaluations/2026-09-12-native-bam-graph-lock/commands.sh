#!/usr/bin/env bash
set -euo pipefail

readonly REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
readonly OUT_ROOT="${OUT_ROOT:-/tmp/pgphase_native_bam_graph_lock}"
readonly TRUTH_ROOT="${TRUTH_ROOT:-${HOME}/Downloads/pgphase-eval-data/truth}"
readonly THREADS="${THREADS:-20}"

cd "${REPO_ROOT}"

run_chromosome() {
    local -r chrom="$1"
    local -r ref="test_data/chm13v2.0.${chrom}.renamed.fa"
    local -r sites="test_data/${chrom}.sites.vcf.gz"
    local -r gaf="test_data/HG002.${chrom}.annotated.coord.gaf.gz"
    local -r bam="test_data/HG002_${chrom}_hifi_mapped_to_CHM13_${chrom}_annotated.bam"
    local -r out="${OUT_ROOT}/${chrom}"
    local -r truth="${TRUTH_ROOT}/${chrom}/diplinator_merged.bam"

    mkdir -p "${out}/graph/eval" "${out}/bam" "${out}/hybrid/eval" \
        "${out}/graph_lock/eval"

    ./pgphase collect-graph-variation \
        --ref "${ref}" \
        --sites "${sites}" \
        --gaf "${gaf}" \
        --min-read-margin 2 \
        --anchor-af-margin 0.12 \
        -r "${chrom}" \
        -t "${THREADS}" \
        -o "${out}/graph/candidates.tsv" \
        --phased-vcf-out "${out}/graph/phased.vcf" \
        --phased-bam-out "${out}/graph/phased.bam"

    ./pgphase collect-bam-variation \
        --ref "${ref}" \
        --bam "${bam}" \
        -r "${chrom}" \
        -t "${THREADS}" \
        -o "${out}/bam/candidates.tsv" \
        --phased-vcf-out "${out}/bam/phased.vcf"

    bgzip -f -k "${out}/graph/phased.vcf"
    bgzip -f -k "${out}/bam/phased.vcf"
    tabix -f -p vcf "${out}/graph/phased.vcf.gz"
    tabix -f -p vcf "${out}/bam/phased.vcf.gz"

    python3 scripts/extract_private_gap_sites.py \
        --graph-phased-vcf "${out}/graph/phased.vcf.gz" \
        --graph-sites "${sites}" \
        --linear-vcf "${out}/bam/phased.vcf.gz" \
        --contig "${chrom}" \
        --output "${out}/private_absent_clean_snp_bridges.vcf" \
        --gaps-bed "${out}/graph_gaps.bed" \
        --min-gq 10 \
        --clean-snps-only \
        --exclude-graph-positions \
        --min-vaf 0.30 \
        --max-vaf 0.70 \
        --bam "${bam}" \
        --min-bridge-reads 2 \
        --min-mapq 20

    ./pgphase collect-hybrid-variation \
        --ref "${ref}" \
        --bam "${bam}" \
        --graph-sites "${sites}" \
        --gaf "${gaf}" \
        --private-sites "${out}/private_absent_clean_snp_bridges.vcf" \
        --min-read-margin 2 \
        --min-phase-set-reads 0 \
        -r "${chrom}" \
        -t "${THREADS}" \
        -o "${out}/hybrid/candidates.tsv" \
        --phased-vcf-out "${out}/hybrid/phased.vcf" \
        -b "${out}/hybrid/phased.bam"

    python3 scripts/merge_graph_hybrid_tags.py \
        --graph-bam "${out}/graph/phased.bam" \
        --hybrid-bam "${out}/hybrid/phased.bam" \
        --output "${out}/graph_lock/phased.bam" \
        --min-shared-reads 10 \
        --min-vote-margin 5 \
        --min-purity 0.90 \
        --require-both-haplotypes \
        --threads 8

    python3 scripts/evaluate_phase_accuracy.py \
        "${out}/graph/phased.bam" "${truth}" 0 0 5 "" \
        "${out}/graph/eval" samtools "" "" "" ""
    python3 scripts/evaluate_phase_accuracy.py \
        "${out}/hybrid/phased.bam" "${truth}" 0 0 5 "" \
        "${out}/hybrid/eval" samtools "" "" "" ""
    python3 scripts/evaluate_phase_accuracy.py \
        "${out}/graph_lock/phased.bam" "${truth}" 0 0 5 "" \
        "${out}/graph_lock/eval" samtools "" "" "" ""
}

make all -j"${THREADS}"
run_chromosome chr18
run_chromosome chr12
