#!/usr/bin/env bash
# Measure --gap-bridge-private-snps on whole chr20 against the same binary with
# the flag off. Both arms reuse the frozen gap evidence cache, so the only
# difference between them is the flag itself.
set -euo pipefail

readonly REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
readonly OUT="${OUT:-/tmp/pgphase-private-snp-bridge}"
readonly DATA_ROOT="${DATA_ROOT:-${HOME}/Downloads/pgphase-eval-data}"
readonly PGPHASE="${PGPHASE:-${REPO}/pgphase}"
readonly CACHE="${CACHE:-/tmp/chr20-evidence-v6.gapev}"
readonly REGION='CHM13#0#chr20'
readonly THREADS="${THREADS:-8}"
readonly REF="${REPO}/test_data/chm13v2.0.chr20.renamed.fa"
readonly BAM="${REPO}/test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam"
readonly SITES="${REPO}/test_data/chr20.sites.striped.vcf.gz"
readonly GAF="${REPO}/test_data/HG002.chr20.annotated.coord.gaf.gz"
readonly TRUTH="${DATA_ROOT}/truth/chr20/diplinator_merged.bam"

die() { echo "error: $*" >&2; exit 1; }
if [[ "${1:-}" == '--help' ]]; then
    echo 'Usage: run.sh [--help]'
    echo 'Run whole-chr20 gap recovery with --gap-bridge-private-snps off and on.'
    echo 'Environment: OUT, DATA_ROOT, PGPHASE, CACHE, THREADS.'
    exit 0
fi
[[ $# -eq 0 ]] || die "unexpected argument: $1"
for input in "${PGPHASE}" "${REF}" "${BAM}" "${SITES}" "${GAF}" "${TRUTH}" "${CACHE}"; do
    [[ -f "${input}" ]] || die "missing input: ${input}"
done

mkdir -p "${OUT}"
sha256sum "${PGPHASE}" > "${OUT}/fingerprints.txt"
(cd "${REPO}" && git rev-parse HEAD && git diff --stat -- src) > "${OUT}/source_state.txt"

for arm in off on; do
    dir="${OUT}/${arm}"
    mkdir -p "${dir}"
    command=("${PGPHASE}" collect-hybrid-variation
             --ref "${REF}" --bam "${BAM}" --graph-sites "${SITES}" --gaf "${GAF}"
             -r "${REGION}" --threads "${THREADS}"
             --link-by-alleles --block-link-window 8 --min-read-margin 2
             --recover-gaps --gap-evidence-cache "${CACHE}"
             -o "${dir}/candidates.tsv" --phased-vcf-out "${dir}/native.vcf"
             -b "${dir}/phased.bam" --gap-recovery-report "${dir}/tiers.tsv")
    if [[ "${arm}" == on ]]; then
        command+=(--gap-bridge-private-snps)
    fi
    printf '%q ' "${command[@]}" > "${dir}/command.sh"
    printf '\n' >> "${dir}/command.sh"
    echo "[$(date +%H:%M:%S)] phasing arm ${arm}"
    "${command[@]}" > "${dir}/stdout.log" 2> "${dir}/stderr.log"

    echo "[$(date +%H:%M:%S)] evaluating arm ${arm}"
    mkdir -p "${dir}/eval"
    python3 "${REPO}/scripts/evaluate_phase_accuracy.py" "${dir}/phased.bam" \
        "${TRUTH}" 0 0 5 '' "${dir}/eval" samtools '' '' '' '' \
        > "${dir}/eval.log" 2>&1
done
echo "[$(date +%H:%M:%S)] done"
