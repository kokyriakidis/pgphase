#!/usr/bin/env bash
set -euo pipefail

usage() {
    echo "Usage: $0 BEFORE_BINARY AFTER_BINARY OUTPUT_DIR"
    echo "Replay owning 5 Mb, owning 23 Mb, paired 22-24 Mb and full chr20."
    echo "Run from the pgphase repository with test_data and the truth map present."
}
die() { echo "error: $*" >&2; exit 1; }
if [[ "${1:-}" == "--help" ]]; then usage; exit 0; fi
[[ "$#" == 3 ]] || { usage >&2; exit 1; }
readonly before_binary="$(realpath "$1")"
readonly after_binary="$(realpath "$2")"
readonly output_dir="$(realpath -m "$3")"
readonly evaluator="evaluations/2026-10-05-composed-physical-bridges/evaluate.py"
readonly parity="evaluations/2026-10-02-bam-pair-mapq/measure_output_parity.py"
readonly python_bin="${PYTHON:-python3}"
[[ -x "${before_binary}" && -x "${after_binary}" ]] || die "Both binaries must be executable"
export LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
for stage in before after; do
    binary="${before_binary}"
    [[ "${stage}" != after ]] || binary="${after_binary}"
    for scope in owning5 owning23 paired23 full; do
        region=()
        case "${scope}" in
            owning5) region=(-r CHM13#0#chr20:5000001-6000000) ;;
            owning23) region=(-r CHM13#0#chr20:23000001-24000000) ;;
            paired23) region=(-r CHM13#0#chr20:22000001-24000000) ;;
        esac
        run_dir="${output_dir}/${stage}/${scope}"
        mkdir -p "${run_dir}"
        "${binary}" collect-graph-variation \
            --ref test_data/chm13v2.0.chr20.renamed.fa \
            --bam test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam \
            --sites test_data/chr20.sites.striped.vcf.gz \
            --gaf test_data/HG002.chr20.annotated.coord.gaf.gz \
            -t 8 "${region[@]}" -o "${run_dir}/candidates.tsv" \
            --phased-vcf-out "${run_dir}/phased.vcf" \
            --phased-bam-out "${run_dir}/phased.bam" \
            > "${run_dir}/stdout.log" 2> "${run_dir}/stderr.log"
    done
done
for scope in owning5 owning23 paired23 full; do
    "${python_bin}" "${parity}" --before "${output_dir}/before/${scope}" \
        --after "${output_dir}/after/${scope}" --output "${output_dir}/${scope}-parity.json"
    if [[ "${scope}" != owning5 ]]; then
        "${python_bin}" "${evaluator}" --before "${output_dir}/before/${scope}" \
            --after "${output_dir}/after/${scope}" --output "${output_dir}/${scope}-audit.json"
    fi
done
