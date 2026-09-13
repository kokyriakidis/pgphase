#!/usr/bin/env bash
set -euo pipefail

readonly REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
readonly OUT="${OUT:-/tmp/pgphase-auto-gap-regression}"
readonly DATA_ROOT="${DATA_ROOT:-${HOME}/Downloads/pgphase-eval-data}"
readonly PGPHASE="${PGPHASE:-${REPO}/pgphase}"
readonly REF="${REPO}/test_data/chm13v2.0.chr20.renamed.fa"
readonly BAM="${REPO}/test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam"
readonly SITES="${REPO}/test_data/chr20.sites.striped.vcf.gz"
readonly GAF="${REPO}/test_data/HG002.chr20.annotated.coord.gaf.gz"

die() { echo "error: $*" >&2; exit 1; }
if [[ "${1:-}" == --help ]]; then
    echo 'Usage: run.sh [--help]'
    echo 'Test automatic gap recovery on two chr20 regions and a split-chunk case.'
    echo 'Environment: OUT, DATA_ROOT, PGPHASE. Requires chr20 test inputs and read truth.'
    exit 0
fi
[[ $# -eq 0 ]] || die "unexpected argument: $1"
mkdir -p "${OUT}"
sha256sum "${PGPHASE}" > "${OUT}/fingerprint.txt"
for case_name in 15m 47m 15m_split; do
    region='CHM13#0#chr20:14950000-15180000'
    chunk_size=500000
    if [[ "${case_name}" == 47m ]]; then region='CHM13#0#chr20:47600000-47820000'; fi
    if [[ "${case_name}" == 15m_split ]]; then chunk_size=80000; fi
    case_out="${OUT}/${case_name}"
    mkdir -p "${case_out}"
    for arm in clean auto24 auto1; do
        extra=()
        if [[ "${arm}" != clean ]]; then
            extra=(--recover-gaps --gap-recovery-report "${case_out}/${arm}.tiers.tsv")
        fi
        if [[ "${arm}" == auto1 ]]; then extra+=(--private-msa-margin 1); fi
        command=("${PGPHASE}" collect-hybrid-variation --ref "${REF}" --bam "${BAM}"
                 --graph-sites "${SITES}" --gaf "${GAF}" -r "${region}"
                 --chunk-size "${chunk_size}" --threads 2
                 --link-by-alleles --block-link-window 8 --min-read-margin 2
                 "${extra[@]}" -o "${case_out}/${arm}.tsv"
                 --phased-vcf-out "${case_out}/${arm}.vcf" -b "${case_out}/${arm}.bam")
        printf '%q ' "${command[@]}" > "${case_out}/${arm}.command.sh"
        printf '\n' >> "${case_out}/${arm}.command.sh"
        "${command[@]}" > "${case_out}/${arm}.log" 2>&1
    done
    samtools view "${case_out}/clean.bam" | cut -f1 | sort -u > "${case_out}/read_names.txt"
    samtools view -b -N "${case_out}/read_names.txt" \
        "${DATA_ROOT}/truth/chr20/diplinator_merged.bam" > "${case_out}/truth.bam"
    for arm in clean auto24 auto1; do
        mkdir -p "${case_out}/${arm}.eval"
        python3 "${REPO}/scripts/evaluate_phase_accuracy.py" "${case_out}/${arm}.bam" \
            "${case_out}/truth.bam" 0 0 5 '' "${case_out}/${arm}.eval" samtools '' '' '' '' \
            > "${case_out}/${arm}.eval.log" 2>&1
    done
done
python3 "${REPO}/evaluations/2026-09-13-auto-gap-recovery/verify.py" --out "${OUT}"
