#!/usr/bin/env bash
# Validate --min-assign-mapq on the three chr20 gaps whose sub-floor reads were
# shown to carry correct phase (see evaluations/2026-09-15-mapq-starved-gaps).
# Each gap runs with the default single floor and with discovery opened to
# MAPQ 1 while haplotype assignment stays at 30 and at 5.
set -euo pipefail

readonly REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
readonly OUT="${OUT:-/tmp/pgphase-split-floor}"
readonly DATA_ROOT="${DATA_ROOT:-${HOME}/Downloads/pgphase-eval-data}"
readonly PGPHASE="${PGPHASE:-${REPO}/pgphase}"
readonly FLANK="${FLANK:-50000}"
readonly THREADS="${THREADS:-8}"
readonly REF="${REPO}/test_data/chm13v2.0.chr20.renamed.fa"
readonly BAM="${REPO}/test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam"
readonly SITES="${REPO}/test_data/chr20.sites.striped.vcf.gz"
readonly GAF="${REPO}/test_data/HG002.chr20.annotated.coord.gaf.gz"
readonly TRUTH="${DATA_ROOT}/truth/chr20/diplinator_merged.bam"
readonly GAPS=(25834662-25883079 25944471-25986123 26029591-26088679)

die() { echo "error: $*" >&2; exit 1; }
if [[ "${1:-}" == '--help' ]]; then
    echo 'Usage: run.sh [--help]'
    echo 'Run the three trustworthy chr20 gaps with single and split MAPQ floors.'
    echo 'Environment: OUT, DATA_ROOT, PGPHASE, FLANK, THREADS.'
    exit 0
fi
[[ $# -eq 0 ]] || die "unexpected argument: $1"
for input in "${PGPHASE}" "${REF}" "${BAM}" "${SITES}" "${GAF}" "${TRUTH}"; do
    [[ -f "${input}" ]] || die "missing input: ${input}"
done
mkdir -p "${OUT}"
sha256sum "${PGPHASE}" > "${OUT}/fingerprints.txt"

for gap in "${GAPS[@]}"; do
    left="${gap%-*}"; right="${gap#*-}"
    region="CHM13#0#chr20:$((left - FLANK))-$((right + FLANK))"
    for arm in default q1_assign30 q1_assign5; do
        dir="${OUT}/${gap}/${arm}"
        mkdir -p "${dir}"
        floors=(-q 30)
        case "${arm}" in
            q1_assign30) floors=(-q 1 --min-assign-mapq 30) ;;
            q1_assign5)  floors=(-q 1 --min-assign-mapq 5) ;;
        esac
        command=("${PGPHASE}" collect-hybrid-variation --ref "${REF}" --bam "${BAM}"
                 --graph-sites "${SITES}" --gaf "${GAF}" -r "${region}"
                 -t "${THREADS}" "${floors[@]}"
                 --link-by-alleles --block-link-window 8 --min-read-margin 2 --recover-gaps
                 -o "${dir}/candidates.tsv" --phased-vcf-out "${dir}/native.vcf"
                 -b "${dir}/phased.bam" --gap-recovery-report "${dir}/tiers.tsv")
        printf '%q ' "${command[@]}" > "${dir}/command.sh"; printf '\n' >> "${dir}/command.sh"
        "${command[@]}" > "${dir}/stdout.log" 2> "${dir}/stderr.log"
    done
    # One truth subset per gap, by read name, after every arm has run.
    samtools view "${OUT}/${gap}/q1_assign5/phased.bam" | cut -f1 | sort -u > "${OUT}/${gap}/read_names.txt"
    samtools view -b -N "${OUT}/${gap}/read_names.txt" "${TRUTH}" > "${OUT}/${gap}/truth.bam"
    for arm in default q1_assign30 q1_assign5; do
        dir="${OUT}/${gap}/${arm}"
        mkdir -p "${dir}/eval"
        python3 "${REPO}/scripts/evaluate_phase_accuracy.py" "${dir}/phased.bam" \
            "${OUT}/${gap}/truth.bam" 0 0 5 '' "${dir}/eval" samtools '' '' '' '' \
            > "${dir}/eval.log" 2>&1
    done
    echo "done ${gap}"
done
