#!/usr/bin/env bash
set -euo pipefail

usage() {
    echo "Usage: $0 BEFORE_BINARY AFTER_BINARY OUTPUT_DIR"
    echo "Replay the bounded 9 Mb seam, paired 8-10 Mb owning chunks and full chr20."
    echo "Run from the pgphase repository with test_data and the truth map present."
}
die() { echo "error: $*" >&2; exit 1; }
if [[ "${1:-}" == "--help" ]]; then usage; exit 0; fi
[[ "$#" == 3 ]] || { usage >&2; exit 1; }
readonly before_binary="$(realpath "$1")"
readonly after_binary="$(realpath "$2")"
readonly output_dir="$(realpath -m "$3")"
readonly evaluator="evaluations/2026-10-05-complementary-insertion-boundary/evaluate.py"
readonly parity="evaluations/2026-10-02-bam-pair-mapq/measure_output_parity.py"
readonly python_bin="${PYTHON:-python3}"
[[ -x "${before_binary}" && -x "${after_binary}" ]] || die "Both binaries must be executable"
export LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
for stage in before after; do
    binary="${before_binary}"
    [[ "${stage}" != after ]] || binary="${after_binary}"
    for scope in bounded paired full; do
        region=()
        matrix=()
        case "${scope}" in
            bounded) region=(-r CHM13#0#chr20:8927829-9064032) ;;
            paired) region=(-r CHM13#0#chr20:8000001-10000000) ;;
        esac
        run_dir="${output_dir}/${stage}/${scope}"
        mkdir -p "${run_dir}"
        [[ "${scope}" != bounded ]] || matrix=(--phase-matrix-dump "${run_dir}/matrix")
        "${binary}" collect-graph-variation \
            --ref test_data/chm13v2.0.chr20.renamed.fa \
            --bam test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam \
            --sites test_data/chr20.sites.striped.vcf.gz \
            --gaf test_data/HG002.chr20.annotated.coord.gaf.gz \
            -t 8 "${region[@]}" "${matrix[@]}" -o "${run_dir}/candidates.tsv" \
            --phased-vcf-out "${run_dir}/phased.vcf" \
            --phased-bam-out "${run_dir}/phased.bam" \
            > "${run_dir}/stdout.log" 2> "${run_dir}/stderr.log"
    done
done
for scope in bounded paired full; do
    "${python_bin}" "${parity}" --before "${output_dir}/before/${scope}" \
        --after "${output_dir}/after/${scope}" --output "${output_dir}/${scope}-parity.json"
    "${python_bin}" "${evaluator}" --before "${output_dir}/before/${scope}" \
        --after "${output_dir}/after/${scope}" --output "${output_dir}/${scope}-audit.json"
done
"${python_bin}" evaluations/2026-10-05-complementary-insertion-boundary/evidence.py \
    --replay "${output_dir}/after/bounded" --output "${output_dir}/physical-evidence.json"
