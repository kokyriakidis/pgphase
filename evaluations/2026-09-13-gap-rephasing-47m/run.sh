#!/usr/bin/env bash
set -euo pipefail

readonly REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
readonly OUT="${OUT:-/tmp/pgphase-gap-rephasing-47m}"
readonly DATA_ROOT="${DATA_ROOT:-${HOME}/Downloads/pgphase-eval-data}"
readonly PGPHASE="${PGPHASE:-${REPO}/pgphase}"
readonly REGION='CHM13#0#chr20:47600000-47820000'
readonly REF="${REPO}/test_data/chm13v2.0.chr20.renamed.fa"
readonly BAM="${REPO}/test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam"
readonly SITES="${REPO}/test_data/chr20.sites.striped.vcf.gz"
readonly GAF="${REPO}/test_data/HG002.chr20.annotated.coord.gaf.gz"
readonly PRIVATE="${OUT}/private_sites.vcf"

die() { echo "error: $*" >&2; exit 1; }
if [[ "${1:-}" == '--help' ]]; then
    echo 'Usage: run.sh [--help]'
    echo 'Run clean and MSA second-pass controls at chr20:47.60-47.82 Mb.'
    echo 'Environment: OUT, DATA_ROOT, PGPHASE. Requires local chr20 test data.'
    exit 0
fi
[[ $# -eq 0 ]] || die "unexpected argument: $1"
for input in "${PGPHASE}" "${REF}" "${BAM}" "${SITES}" "${GAF}"; do
    [[ -f "${input}" ]] || die "missing input: ${input}"
done
mkdir -p "${OUT}"
"${PGPHASE}" collect-bam-variation --ref "${REF}" --bam "${BAM}" -r "${REGION}" \
    -o "${OUT}/native.tsv" --phased-vcf-out "${OUT}/native.vcf" \
    -b "${OUT}/native.bam" > "${OUT}/native.log" 2>&1
python3 - "${OUT}/native.vcf" "${PRIVATE}" <<'PYVCF'
import sys
import pysam
with pysam.VariantFile(sys.argv[1]) as source, pysam.VariantFile(sys.argv[2], 'w', header=source.header) as dest:
    for record in source:
        if not 47666540 <= record.pos <= 47767233:
            continue
        sample = record.samples[0]
        gt = sample.get('GT')
        ad = sample.get('AD')
        if len(record.alleles) != 2 or gt is None or None in gt or gt[0] == gt[1]:
            continue
        if (sample.get('GQ') or 0) < 10 or not ad or sum(ad) == 0:
            continue
        if not 0.30 <= ad[1] / sum(ad) <= 0.70:
            continue
        sample['GT'] = (0, 1)
        sample.phased = False
        if 'PS' in sample:
            sample['PS'] = None
        dest.write(record)
PYVCF
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

readonly FROZEN="${DATA_ROOT}/results/chr12-18-20-comparison/chr20"
samtools view -b "${FROZEN}/graph/tagged_surjected.bam" "${REGION}" > "${OUT}/graph.bam"
samtools view -b "${FROZEN}/hiphase/phased.bam" 'chr20:47600000-47820000' > "${OUT}/hiphase.bam"

for arm in msa msa_margin1; do
    distance=300000
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
python3 "${REPO}/evaluations/2026-09-13-gap-rephasing-47m/summarize.py" \
    --out "${OUT}" --data-root "${DATA_ROOT}"
