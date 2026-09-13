#!/usr/bin/env bash
set -euo pipefail

readonly REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
readonly DATA_ROOT="${DATA_ROOT:-${HOME}/Downloads/pgphase-eval-data}"
readonly SHARED_VCF_ROOT="${SHARED_VCF_ROOT:-${DATA_ROOT}/shared_calls}"
readonly TRUTH_VCF="${TRUTH_VCF:-${DATA_ROOT}/truth_chm13/HG002_CHM13v2.0_v5.0q_smvar.vcf.gz}"
readonly OUT_ROOT="${OUT_ROOT:-${DATA_ROOT}/results/chr12-18-20-comparison}"
readonly PHASER_BIN="${PHASER_BIN:-${HOME}/micromamba/envs/bench-phasers/bin}"
readonly THREADS="${THREADS:-20}"
readonly CHROMS="${CHROMS:-chr12 chr18 chr20}"

die() {
    echo "error: $*" >&2
    exit 1
}

usage() {
    cat <<'EOF'
Usage: commands.sh

Run the chr12/chr18/chr20 shared-call phasing comparison. Environment:
  CHROMS           space-separated chromosomes [chr12 chr18 chr20]
  DATA_ROOT         durable evaluation root [~/Downloads/pgphase-eval-data]
  SHARED_VCF_ROOT   DeepVariant callset root [$DATA_ROOT/shared_calls]
  OUT_ROOT          result root [$DATA_ROOT/results/chr12-18-20-comparison]
  PHASER_BIN        directory containing whatshap/hiphase/longphase/bcftools
  THREADS           worker threads [20]
EOF
}

if [[ "${1:-}" == "--help" ]]; then
    usage
    exit 0
fi
[[ $# -eq 0 ]] || die "unexpected argument: $1"

readonly BCFTOOLS="${PHASER_BIN}/bcftools"
readonly BGZIP="${PHASER_BIN}/bgzip"
readonly TABIX="${PHASER_BIN}/tabix"
readonly WHATSHAP="${PHASER_BIN}/whatshap"
readonly HIPHASE="${PHASER_BIN}/hiphase"
readonly LONGPHASE="${PHASER_BIN}/longphase"

for tool in "${BCFTOOLS}" "${BGZIP}" "${TABIX}" "${WHATSHAP}" "${HIPHASE}" "${LONGPHASE}"; do
    [[ -x "${tool}" ]] || die "missing executable: ${tool}"
done
[[ -x "${REPO_ROOT}/pgphase" ]] || die "missing pgphase binary; run make all"
[[ -f "${TRUTH_VCF}" ]] || die "missing truth VCF: ${TRUTH_VCF}"

run_timed() {
    local -r resource_file="$1"
    shift
    /usr/bin/time -v -o "${resource_file}" "$@"
}

compress_vcf() {
    local -r input="$1"
    local -r output="$2"
    "${BGZIP}" -f -c "${input}" > "${output}"
    "${TABIX}" -f -p vcf "${output}"
}

evaluate_reads() {
    local -r bam="$1"
    local -r truth="$2"
    local -r out="$3"
    mkdir -p "${out}"
    python3 "${REPO_ROOT}/scripts/evaluate_phase_accuracy.py" \
        "${bam}" "${truth}" 0 0 5 "" "${out}" samtools "" "" "" ""
}

evaluate_vcf() {
    local -r chrom="$1"
    local -r length="$2"
    local -r label="$3"
    local -r vcf="$4"
    local -r out="$5"
    "${WHATSHAP}" compare \
        --tsv-pairwise "${out}/${label}.compare.tsv" \
        --switch-error-bed "${out}/${label}.switch_errors.bed" \
        --names truth,"${label}" \
        "${out}/truth.vcf.gz" "${vcf}" \
        > "${out}/${label}.compare.txt"
    "${WHATSHAP}" stats \
        --tsv "${out}/${label}.stats.tsv" \
        --block-list "${out}/${label}.blocks.tsv" \
        "${vcf}" > "${out}/${label}.stats.txt"
    python3 "${REPO_ROOT}/scripts/compute_ngc50.py" \
        "${out}/${label}.blocks.tsv" "${out}/${label}.switch_errors.bed" \
        "${out}/${label}.ngc50.json" --genome-size "${length}"
    printf '%s\t%s\t%s\n' "${chrom}" "${label}" "${vcf}" >> "${out}/evaluated_vcfs.tsv"
}

run_chromosome() {
    local -r chrom="$1"
    local -r ref="${REPO_ROOT}/test_data/chm13v2.0.${chrom}.renamed.fa"
    local -r sites_name="$([[ "${chrom}" == chr20 ]] && echo chr20.sites.striped.vcf.gz || echo "${chrom}.sites.vcf.gz")"
    local -r sites="${REPO_ROOT}/test_data/${sites_name}"
    local -r gaf="${REPO_ROOT}/test_data/HG002.${chrom}.annotated.coord.gaf.gz"
    local -r bam="${REPO_ROOT}/test_data/HG002_${chrom}_hifi_mapped_to_CHM13_${chrom}_annotated.bam"
    local linear_ref="${ref}"
    local linear_bam="${bam}"
    local bam_region="${chrom}"
    local graph_contig="${chrom}"
    local -r truth_bam="${DATA_ROOT}/truth/${chrom}/diplinator_merged.bam"
    local -r shared_vcf="${SHARED_VCF_ROOT}/${chrom}/deepvariant.vcf.gz"
    local -r out="${OUT_ROOT}/${chrom}"
    local length
    length="$(cut -f2 "${ref}.fai")"
    if [[ "${chrom}" == "chr20" ]]; then
        linear_ref="${SHARED_VCF_ROOT}/${chrom}/chm13v2.0.${chrom}.normalized.fa"
        linear_bam="${SHARED_VCF_ROOT}/${chrom}/HG002.${chrom}.normalized.bam"
        bam_region="CHM13#0#chr20"
        graph_contig="CHM13#0#chr20"
    fi

    for input in "${ref}" "${sites}" "${gaf}" "${bam}" "${linear_ref}" "${linear_bam}" \
        "${truth_bam}" "${shared_vcf}"; do
        [[ -f "${input}" ]] || die "[${chrom}] missing input: ${input}"
    done
    mkdir -p "${out}"/{graph,bam,hybrid,graph_lock,whatshap,whatshap_opt,hiphase,longphase,eval}
    : > "${out}/evaluated_vcfs.tsv"

    "${BCFTOOLS}" view -r "${chrom}" -Oz -o "${out}/truth.vcf.gz" "${TRUTH_VCF}"
    "${TABIX}" -f -p vcf "${out}/truth.vcf.gz"

    if [[ ! -s "${out}/graph/phased.bam" ]]; then
        run_timed "${out}/graph/resources.txt" "${REPO_ROOT}/pgphase" collect-graph-variation \
            --ref "${ref}" --sites "${sites}" --gaf "${gaf}" -r "${chrom}" \
            --min-read-margin 2 --anchor-af-margin 0.12 -t "${THREADS}" \
            -o "${out}/graph/candidates.tsv" \
            --filtered-sites-out "${out}/graph/filtered_sites.tsv" \
            --phase-sites-out "${out}/graph/phase_sites.tsv" \
            --phased-vcf-out "${out}/graph/native.vcf" \
            --phased-bam-out "${out}/graph/phased.bam"
    fi
    samtools index -@ "${THREADS}" "${out}/graph/phased.bam" 2>/dev/null || true
    compress_vcf "${out}/graph/native.vcf" "${out}/graph/native.vcf.gz"

    if [[ ! -s "${out}/bam/native.vcf" ]]; then
        run_timed "${out}/bam/resources.txt" "${REPO_ROOT}/pgphase" collect-bam-variation \
            --hifi --ref "${ref}" --bam "${bam}" -r "${bam_region}" -t "${THREADS}" \
            -o "${out}/bam/candidates.tsv" --phased-vcf-out "${out}/bam/native.vcf" \
            -b "${out}/bam/phased.bam"
    fi
    compress_vcf "${out}/bam/native.vcf" "${out}/bam/native.vcf.gz"

    python3 "${REPO_ROOT}/scripts/extract_private_gap_sites.py" \
        --graph-phased-vcf "${out}/graph/native.vcf.gz" --graph-sites "${sites}" \
        --linear-vcf "${out}/bam/native.vcf.gz" --contig "${graph_contig}" \
        --output "${out}/private_gap_sites.vcf" --gaps-bed "${out}/graph_gaps.bed" \
        --min-gq 10 --clean-snps-only --exclude-graph-positions \
        --min-vaf 0.30 --max-vaf 0.70 --bam "${bam}" \
        --min-bridge-reads 2 --min-mapq 20

    local private_sites_hash
    private_sites_hash="$(sha256sum "${out}/private_gap_sites.vcf" | cut -d' ' -f1)"
    local previous_private_sites_hash=""
    if [[ -s "${out}/hybrid/private_sites.sha256" ]]; then
        previous_private_sites_hash="$(<"${out}/hybrid/private_sites.sha256")"
    fi
    if [[ ! -s "${out}/hybrid/phased.bam" || \
          "${private_sites_hash}" != "${previous_private_sites_hash}" ]]; then
        run_timed "${out}/hybrid/resources.txt" "${REPO_ROOT}/pgphase" collect-hybrid-variation \
            --ref "${ref}" --bam "${bam}" --graph-sites "${sites}" --gaf "${gaf}" \
            --private-sites "${out}/private_gap_sites.vcf" -r "${bam_region}" \
            --min-read-margin 2 --min-phase-set-reads 50 -t "${THREADS}" \
            -o "${out}/hybrid/candidates.tsv" --phased-vcf-out "${out}/hybrid/native.vcf" \
            -b "${out}/hybrid/phased.bam"
        printf '%s\n' "${private_sites_hash}" > "${out}/hybrid/private_sites.sha256"
    fi
    samtools index -@ "${THREADS}" "${out}/hybrid/phased.bam" 2>/dev/null || true

    python3 "${REPO_ROOT}/scripts/merge_graph_hybrid_tags.py" \
        --graph-bam "${out}/graph/phased.bam" --hybrid-bam "${out}/hybrid/phased.bam" \
        --output "${out}/graph_lock/phased.bam" --min-shared-reads 10 \
        --min-vote-margin 5 --min-purity 0.90 --require-both-haplotypes --threads 8
    samtools index -@ "${THREADS}" "${out}/graph_lock/phased.bam" 2>/dev/null || true

    if [[ ! -s "${out}/graph/tagged_surjected.bam" ]]; then
        python3 "${REPO_ROOT}/scripts/inject_hp_tags.py" \
            "${out}/graph/phased.bam" "${bam}" "${out}/graph/tagged_surjected.bam" samtools
    fi
    [[ -s "${out}/graph/tagged_surjected.bam.bai" ]] || \
        samtools index -@ "${THREADS}" "${out}/graph/tagged_surjected.bam"
    python3 "${REPO_ROOT}/scripts/phase_vcf_from_hp.py" \
        "${out}/graph/tagged_surjected.bam" "${shared_vcf}" "${out}/graph/shared.vcf" \
        --region "${chrom}" --min-reads 2 --min-ratio 0.70 \
        --support-cache "${out}/graph/shared.support.tsv"
    compress_vcf "${out}/graph/shared.vcf" "${out}/graph/shared.vcf.gz"
    samtools index -@ "${THREADS}" "${out}/bam/phased.bam" 2>/dev/null || true
    python3 "${REPO_ROOT}/scripts/phase_vcf_from_hp.py" \
        "${out}/bam/phased.bam" "${shared_vcf}" "${out}/bam/shared.vcf" \
        --region "${chrom}" --min-reads 2 --min-ratio 0.70 \
        --support-cache "${out}/bam/shared.support.tsv"
    compress_vcf "${out}/bam/shared.vcf" "${out}/bam/shared.vcf.gz"
    python3 "${REPO_ROOT}/scripts/phase_vcf_from_hp.py" \
        "${out}/hybrid/phased.bam" "${shared_vcf}" "${out}/hybrid/shared.vcf" \
        --region "${chrom}" --min-reads 2 --min-ratio 0.70 \
        --support-cache "${out}/hybrid/shared.support.tsv"
    compress_vcf "${out}/hybrid/shared.vcf" "${out}/hybrid/shared.vcf.gz"
    python3 "${REPO_ROOT}/scripts/phase_vcf_from_hp.py" \
        "${out}/graph_lock/phased.bam" "${shared_vcf}" "${out}/graph_lock/shared.vcf" \
        --region "${chrom}" --min-reads 2 --min-ratio 0.70 \
        --support-cache "${out}/graph_lock/shared.support.tsv"
    compress_vcf "${out}/graph_lock/shared.vcf" "${out}/graph_lock/shared.vcf.gz"

    if [[ ! -s "${out}/whatshap/phased.vcf.gz" ]]; then
        run_timed "${out}/whatshap/resources.txt" "${WHATSHAP}" phase --reference "${linear_ref}" \
            --ignore-read-groups --indels -o "${out}/whatshap/phased.vcf.gz" "${shared_vcf}" "${linear_bam}"
    fi
    "${TABIX}" -f -p vcf "${out}/whatshap/phased.vcf.gz"
    if [[ ! -s "${out}/whatshap_opt/phased.vcf.gz" ]]; then
        run_timed "${out}/whatshap_opt/resources.txt" "${WHATSHAP}" phase --reference "${linear_ref}" \
            --ignore-read-groups --indels --distrust-genotypes \
            -o "${out}/whatshap_opt/phased.vcf.gz" "${shared_vcf}" "${linear_bam}"
    fi
    "${TABIX}" -f -p vcf "${out}/whatshap_opt/phased.vcf.gz"
    if [[ ! -s "${out}/hiphase/phased.vcf.gz" || ! -s "${out}/hiphase/phased.bam" ]]; then
        run_timed "${out}/hiphase/resources.txt" "${HIPHASE}" --reference "${linear_ref}" \
            --bam "${linear_bam}" --vcf "${shared_vcf}" --threads "${THREADS}" \
            --output-bam "${out}/hiphase/phased.bam" --output-vcf "${out}/hiphase/phased.vcf.gz" \
            --stats-file "${out}/hiphase/stats.native.tsv" \
            --blocks-file "${out}/hiphase/blocks.native.tsv" \
            --summary-file "${out}/hiphase/summary.native.tsv"
    fi
    "${TABIX}" -f -p vcf "${out}/hiphase/phased.vcf.gz"
    if [[ ! -s "${out}/longphase/phased.vcf" ]]; then
        run_timed "${out}/longphase/resources.txt" "${LONGPHASE}" phase -s "${shared_vcf}" \
            -b "${linear_bam}" -r "${linear_ref}" -t "${THREADS}" --pb -o "${out}/longphase/phased"
    fi
    compress_vcf "${out}/longphase/phased.vcf" "${out}/longphase/phased.vcf.gz"

    if [[ ! -s "${out}/whatshap/phased.bam" ]]; then
        "${WHATSHAP}" haplotag --reference "${linear_ref}" --output-threads "${THREADS}" \
            -o "${out}/whatshap/phased.bam" "${out}/whatshap/phased.vcf.gz" "${linear_bam}"
    fi
    if [[ ! -s "${out}/whatshap_opt/phased.bam" ]]; then
        "${WHATSHAP}" haplotag --reference "${linear_ref}" --output-threads "${THREADS}" \
            -o "${out}/whatshap_opt/phased.bam" "${out}/whatshap_opt/phased.vcf.gz" "${linear_bam}"
    fi
    if [[ ! -s "${out}/longphase/phased.bam" ]]; then
        "${LONGPHASE}" haplotag -s "${out}/longphase/phased.vcf.gz" -b "${linear_bam}" \
            -r "${linear_ref}" -t "${THREADS}" -o "${out}/longphase/phased"
    fi

    evaluate_reads "${out}/graph/phased.bam" "${truth_bam}" "${out}/eval/graph_reads"
    evaluate_reads "${out}/bam/phased.bam" "${truth_bam}" "${out}/eval/bam_reads"
    evaluate_reads "${out}/hybrid/phased.bam" "${truth_bam}" "${out}/eval/hybrid_reads"
    evaluate_reads "${out}/graph_lock/phased.bam" "${truth_bam}" "${out}/eval/graph_lock_reads"
    evaluate_reads "${out}/whatshap/phased.bam" "${truth_bam}" "${out}/eval/whatshap_reads"
    evaluate_reads "${out}/whatshap_opt/phased.bam" "${truth_bam}" "${out}/eval/whatshap_opt_reads"
    evaluate_reads "${out}/hiphase/phased.bam" "${truth_bam}" "${out}/eval/hiphase_reads"
    evaluate_reads "${out}/longphase/phased.bam" "${truth_bam}" "${out}/eval/longphase_reads"

    evaluate_vcf "${chrom}" "${length}" graph "${out}/graph/shared.vcf.gz" "${out}"
    evaluate_vcf "${chrom}" "${length}" bam "${out}/bam/shared.vcf.gz" "${out}"
    evaluate_vcf "${chrom}" "${length}" hybrid "${out}/hybrid/shared.vcf.gz" "${out}"
    evaluate_vcf "${chrom}" "${length}" graph_lock "${out}/graph_lock/shared.vcf.gz" "${out}"
    evaluate_vcf "${chrom}" "${length}" whatshap "${out}/whatshap/phased.vcf.gz" "${out}"
    evaluate_vcf "${chrom}" "${length}" whatshap_opt "${out}/whatshap_opt/phased.vcf.gz" "${out}"
    evaluate_vcf "${chrom}" "${length}" hiphase "${out}/hiphase/phased.vcf.gz" "${out}"
    evaluate_vcf "${chrom}" "${length}" longphase "${out}/longphase/phased.vcf.gz" "${out}"

    python3 "${REPO_ROOT}/scripts/analyze_graph_gap_site_loss.py" \
        --truth-vcf "${out}/truth.vcf.gz" --phased-vcf "${out}/graph/shared.vcf.gz" \
        --catalog-vcf "${sites}" --candidates-tsv "${out}/graph/candidates.tsv" \
        --filtered-tsv "${out}/graph/filtered_sites.tsv" \
        --phase-sites-tsv "${out}/graph/phase_sites.tsv" --contig "${chrom}" \
        --competitor-vcf "hiphase=${out}/hiphase/phased.vcf.gz" \
        --competitor-vcf "longphase=${out}/longphase/phased.vcf.gz" \
        --competitor-vcf "whatshap=${out}/whatshap/phased.vcf.gz" \
        --recovery-vcf "hybrid=${out}/hybrid/shared.vcf.gz" \
        --recovery-vcf "graph_lock=${out}/graph_lock/shared.vcf.gz" \
        --out "${out}/gap_site_loss.tsv" > "${out}/gap_site_loss.summary.txt"
}

cd "${REPO_ROOT}"
mkdir -p "${OUT_ROOT}"
for chrom in ${CHROMS}; do
    echo "===== ${chrom} ====="
    run_chromosome "${chrom}"
done

completed_chroms=()
for chrom in chr12 chr18 chr20; do
    if [[ -s "${OUT_ROOT}/${chrom}/gap_site_loss.summary.txt" && \
          -s "${OUT_ROOT}/${chrom}/eval/longphase_reads/summary.json" ]]; then
        completed_chroms+=("${chrom}")
    fi
done
python3 "$(dirname "${BASH_SOURCE[0]}")/summarize.py" \
    "${OUT_ROOT}" "${completed_chroms[@]}"
