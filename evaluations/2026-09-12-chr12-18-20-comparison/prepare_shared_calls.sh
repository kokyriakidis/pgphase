#!/usr/bin/env bash
set -euo pipefail

readonly REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
readonly DATA_ROOT="${DATA_ROOT:-${HOME}/Downloads/pgphase-eval-data}"
readonly OUT_ROOT="${SHARED_VCF_ROOT:-${DATA_ROOT}/shared_calls}"
readonly DEEPVARIANT_IMAGE="${DEEPVARIANT_IMAGE:-google/deepvariant:1.10.0}"
readonly THREADS="${THREADS:-20}"
readonly CHROMS="${CHROMS:-chr12 chr18 chr20}"

die() {
    echo "error: $*" >&2
    exit 1
}

usage() {
    cat <<'EOF'
Usage: prepare_shared_calls.sh

Generate one DeepVariant PACBIO callset for each chromosome. Environment:
  CHROMS              space-separated chromosomes [chr12 chr18 chr20]
  DATA_ROOT            durable evaluation root [~/Downloads/pgphase-eval-data]
  SHARED_VCF_ROOT      output root [$DATA_ROOT/shared_calls]
  DEEPVARIANT_IMAGE    Docker image [google/deepvariant:1.10.0]
  THREADS              DeepVariant shards [20]
EOF
}

if [[ "${1:-}" == "--help" ]]; then
    usage
    exit 0
fi
[[ $# -eq 0 ]] || die "unexpected argument: $1"
command -v docker >/dev/null || die "docker is required"
command -v samtools >/dev/null || die "samtools is required"

mkdir -p "${OUT_ROOT}"
for chrom in ${CHROMS}; do
    source_ref="${REPO_ROOT}/test_data/chm13v2.0.${chrom}.renamed.fa"
    source_bam="${REPO_ROOT}/test_data/HG002_${chrom}_hifi_mapped_to_CHM13_${chrom}_annotated.bam"
    out="${OUT_ROOT}/${chrom}"
    [[ -f "${source_ref}" ]] || die "missing reference: ${source_ref}"
    [[ -f "${source_bam}" ]] || die "missing BAM: ${source_bam}"
    mkdir -p "${out}"

    docker_ref="/repo/test_data/chm13v2.0.${chrom}.renamed.fa"
    docker_bam="/repo/test_data/HG002_${chrom}_hifi_mapped_to_CHM13_${chrom}_annotated.bam"
    if [[ "${chrom}" == "chr20" ]]; then
        normalized_ref="${out}/chm13v2.0.${chrom}.normalized.fa"
        normalized_bam="${out}/HG002.${chrom}.normalized.bam"
        if [[ ! -s "${normalized_ref}" ]]; then
            echo "[${chrom}] normalizing CHM13 graph-path FASTA header"
            sed '1s/^>CHM13#0#chr20$/>chr20/' "${source_ref}" > "${normalized_ref}"
            samtools faidx "${normalized_ref}"
        fi
        if [[ ! -s "${normalized_bam}" ]]; then
            echo "[${chrom}] normalizing CHM13 graph-path BAM headers"
            samtools reheader -c 'sed s/SN:CHM13#0#chr/SN:chr/g' \
                "${source_bam}" > "${normalized_bam}"
            samtools index -@ "${THREADS}" "${normalized_bam}"
        fi
        [[ "$(cut -f1 "${normalized_ref}.fai")" == "chr20" ]] || \
            die "normalized FASTA does not contain chr20"
        samtools view -H "${normalized_bam}" | grep -q $'SN:chr20\t' || \
            die "normalized BAM does not contain chr20"
        docker_ref="/output/${chrom}/chm13v2.0.${chrom}.normalized.fa"
        docker_bam="/output/${chrom}/HG002.${chrom}.normalized.bam"
    fi
    if [[ -s "${out}/deepvariant.vcf.gz" && -s "${out}/deepvariant.vcf.gz.tbi" ]]; then
        echo "[${chrom}] shared VCF exists; skipping"
        continue
    fi
    echo "[${chrom}] DeepVariant PACBIO"
    /usr/bin/time -v -o "${out}/deepvariant.resources.txt" \
        docker run --rm \
        -v "${REPO_ROOT}:/repo:ro" \
        -v "${OUT_ROOT}:/output" \
        "${DEEPVARIANT_IMAGE}" \
        /opt/deepvariant/bin/run_deepvariant \
        --model_type=PACBIO \
        --ref="${docker_ref}" \
        --reads="${docker_bam}" \
        --regions="${chrom}" \
        --output_vcf="/output/${chrom}/deepvariant.vcf.gz" \
        --num_shards="${THREADS}" \
        --intermediate_results_dir="/output/${chrom}/intermediate"
done
