#!/usr/bin/env bash
set -euo pipefail

readonly REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
readonly OUT="${OUT:-/tmp/pgphase-gap-rephasing}"
readonly DATA_ROOT="${DATA_ROOT:-${HOME}/Downloads/pgphase-eval-data}"
readonly PGPHASE="${PGPHASE:-${REPO}/pgphase}"
readonly REGION='CHM13#0#chr20:14950000-15180000'
readonly REF="${REPO}/test_data/chm13v2.0.chr20.renamed.fa"
readonly BAM="${REPO}/test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam"
readonly SITES="${REPO}/test_data/chr20.sites.striped.vcf.gz"
readonly GAF="${REPO}/test_data/HG002.chr20.annotated.coord.gaf.gz"
readonly PRIVATE="${REPO}/evaluations/2026-09-13-gap-rephasing/private_sites.vcf"

die() { echo "error: $*" >&2; exit 1; }
if [[ "${1:-}" == '--help' ]]; then
    echo 'Usage: test_chr20_gap_rephasing.sh [--help]'
    echo 'Run clean and MSA second-pass controls at chr20:14.95-15.18 Mb.'
    echo 'Environment: OUT, DATA_ROOT, PGPHASE. Requires local chr20 test data.'
    exit 0
fi
[[ $# -eq 0 ]] || die "unexpected argument: $1"
for input in "${PGPHASE}" "${REF}" "${BAM}" "${SITES}" "${GAF}" "${PRIVATE}"; do
    [[ -f "${input}" ]] || die "missing input: ${input}"
done
mkdir -p "${OUT}"
sha256sum "${PGPHASE}" "${PRIVATE}" > "${OUT}/fingerprints.txt"
for arm in clean msa msa_margin1; do
    extra=()
    if [[ "${arm}" != clean ]]; then
        extra=(--private-sites "${PRIVATE}" --private-msa
               --private-msa-admit-all-in-region --private-msa-snp-first)
    fi
    if [[ "${arm}" == msa_margin1 ]]; then
        extra+=(--private-msa-margin 1)
    fi
    command=("${PGPHASE}" collect-hybrid-variation --ref "${REF}" --bam "${BAM}"
             --graph-sites "${SITES}" --gaf "${GAF}" -r "${REGION}"
             --link-by-alleles --block-link-window 8 --min-read-margin 2
             "${extra[@]}" -o "${OUT}/${arm}.tsv"
             --phased-vcf-out "${OUT}/${arm}.vcf" -b "${OUT}/${arm}.bam"
             --phase-matrix-dump "${OUT}/${arm}.matrix")
    printf '%q ' "${command[@]}" > "${OUT}/${arm}.command.sh"
    printf '\n' >> "${OUT}/${arm}.command.sh"
    "${command[@]}" > "${OUT}/${arm}.log" 2>&1
    samtools index "${OUT}/${arm}.bam"
done

"${PGPHASE}" collect-bam-variation --ref "${REF}" --bam "${BAM}" -r "${REGION}" \
    -o "${OUT}/native.tsv" --phased-vcf-out "${OUT}/native.vcf" \
    -b "${OUT}/native.bam" > "${OUT}/native.log" 2>&1
readonly FROZEN="${DATA_ROOT}/results/chr12-18-20-comparison/chr20"
samtools view -b "${FROZEN}/graph/tagged_surjected.bam" "${REGION}" > "${OUT}/graph.bam"
samtools view -b "${FROZEN}/hiphase/phased.bam" 'chr20:14950000-15180000' > "${OUT}/hiphase.bam"

for arm in msa msa_margin1; do
    distance=300000
    # The frozen left PS starts at 14,729,749, so its ID is 385,638 bp from
    # the right PS although the audited gap is only 110,783 bp. This local
    # diagnostic explicitly admits that endpoint pair; it is not a panel default.
    if [[ "${arm}" == msa_margin1 ]]; then distance=400000; fi
    python3 "${REPO}/scripts/merge_graph_hybrid_tags.py" \
        --graph-bam "${OUT}/graph.bam" --hybrid-bam "${OUT}/${arm}.bam" \
        --output "${OUT}/stitched_${arm}.bam" --min-shared-reads 10 \
        --min-vote-margin 5 --min-purity 0.90 --require-both-haplotypes \
        --merge-graph-phase-sets --max-graph-bridge-distance "${distance}" \
        --graph-bridge-report "${OUT}/stitched_${arm}.edges.tsv" \
        > "${OUT}/stitched_${arm}.log"
done

# Restrict truth by read name only, after every phasing run. No truth enters
# candidate admission or stitching. Cache one truth scan for all local arms.
samtools view "${OUT}/clean.bam" | cut -f1 | sort -u > "${OUT}/read_names.txt"
samtools view -b -N "${OUT}/read_names.txt" \
    "${DATA_ROOT}/truth/chr20/diplinator_merged.bam" > "${OUT}/truth.bam"
for arm in graph hiphase native clean msa msa_margin1 stitched_msa stitched_msa_margin1; do
    mkdir -p "${OUT}/${arm}.eval"
    python3 "${REPO}/scripts/evaluate_phase_accuracy.py" "${OUT}/${arm}.bam" \
        "${OUT}/truth.bam" 0 0 5 '' "${OUT}/${arm}.eval" samtools '' '' '' '' \
        > "${OUT}/${arm}.eval.log" 2>&1
done
python3 "${REPO}/scripts/summarize_chr20_gap_rephasing.py" \
    --out "${OUT}" --repo "${REPO}" --data-root "${DATA_ROOT}"
