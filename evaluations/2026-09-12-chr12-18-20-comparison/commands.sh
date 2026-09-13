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
readonly RUN_COMPETITORS="${RUN_COMPETITORS:-0}"
readonly FORCE_PGPHASE="${FORCE_PGPHASE:-0}"
readonly LONGCALLD="${LONGCALLD:-${HOME}/Downloads/longcallD/bin/longcallD}"

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
  RUN_COMPETITORS   explicitly allow missing competitor artifacts to run [0]
  FORCE_PGPHASE     bypass pgphase stage signatures [0]
  LONGCALLD         LongcallD executable [~/Downloads/longcallD/bin/longcallD]
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

for tool in "${BCFTOOLS}" "${BGZIP}" "${TABIX}" "${WHATSHAP}"; do
    [[ -x "${tool}" ]] || die "missing executable: ${tool}"
done
if [[ "${RUN_COMPETITORS}" == "1" ]]; then
    for tool in "${HIPHASE}" "${LONGPHASE}"; do
        [[ -x "${tool}" ]] || die "missing competitor executable: ${tool}"
    done
fi
[[ -x "${REPO_ROOT}/pgphase" ]] || die "missing pgphase binary; run make all"
[[ -f "${TRUTH_VCF}" ]] || die "missing truth VCF: ${TRUTH_VCF}"

run_timed() {
    local -r resource_file="$1"
    shift
    /usr/bin/time -v -o "${resource_file}" "$@"
}

cached_step() {
    local force_args=()
    if [[ "${FORCE_PGPHASE}" == "1" ]]; then
        force_args+=(--force)
    fi
    python3 "${REPO_ROOT}/scripts/run_cached_step.py" "${force_args[@]}" "$@"
}

allow_competitor_run() {
    local -r artifact="$1"
    [[ "${RUN_COMPETITORS}" == "1" ]] || die \
        "frozen competitor artifact missing: ${artifact}; restore it or explicitly set RUN_COMPETITORS=1"
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
    cached_step --state "${out}/step.json" \
        --input "${REPO_ROOT}/scripts/evaluate_phase_accuracy.py" \
        --input "${bam}" --input "${truth}" \
        --output "${out}/summary.json" --output "${out}/per_phase_set.tsv" -- \
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
    mkdir -p "${out}"/{graph,bam,hybrid,graph_lock,graph_bridge,whatshap,whatshap_opt,hiphase,longphase,longcalld,eval}
    : > "${out}/evaluated_vcfs.tsv"

    if [[ ! -s "${out}/truth.vcf.gz" ]]; then
        "${BCFTOOLS}" view -r "${chrom}" -Oz -o "${out}/truth.vcf.gz" "${TRUTH_VCF}"
        "${TABIX}" -f -p vcf "${out}/truth.vcf.gz"
    fi

    cached_step --state "${out}/graph/step.json" \
        --input "${REPO_ROOT}/pgphase" --input "${ref}" --input "${sites}" --input "${gaf}" \
        --output "${out}/graph/phased.bam" --output "${out}/graph/native.vcf" \
        --output "${out}/graph/candidates.tsv" --output "${out}/graph/filtered_sites.tsv" \
        --output "${out}/graph/phase_sites.tsv" -- \
        /usr/bin/time -v -o "${out}/graph/resources.txt" \
        "${REPO_ROOT}/pgphase" collect-graph-variation \
            --ref "${ref}" --sites "${sites}" --gaf "${gaf}" -r "${chrom}" \
            --min-read-margin 2 --anchor-af-margin 0.12 -t "${THREADS}" \
            -o "${out}/graph/candidates.tsv" \
            --filtered-sites-out "${out}/graph/filtered_sites.tsv" \
            --phase-sites-out "${out}/graph/phase_sites.tsv" \
            --phased-vcf-out "${out}/graph/native.vcf" \
            --phased-bam-out "${out}/graph/phased.bam"
    samtools index -@ "${THREADS}" "${out}/graph/phased.bam" 2>/dev/null || true
    compress_vcf "${out}/graph/native.vcf" "${out}/graph/native.vcf.gz"

    cached_step --state "${out}/bam/step.json" \
        --input "${REPO_ROOT}/pgphase" --input "${ref}" --input "${bam}" \
        --output "${out}/bam/phased.bam" --output "${out}/bam/native.vcf" \
        --output "${out}/bam/candidates.tsv" -- \
        /usr/bin/time -v -o "${out}/bam/resources.txt" \
        "${REPO_ROOT}/pgphase" collect-bam-variation \
            --hifi --ref "${ref}" --bam "${bam}" -r "${bam_region}" -t "${THREADS}" \
            -o "${out}/bam/candidates.tsv" --phased-vcf-out "${out}/bam/native.vcf" \
            -b "${out}/bam/phased.bam"
    samtools index -@ "${THREADS}" "${out}/bam/phased.bam" 2>/dev/null || true
    compress_vcf "${out}/bam/native.vcf" "${out}/bam/native.vcf.gz"

    cached_step --state "${out}/private_gap_sites.step.json" \
        --input "${REPO_ROOT}/scripts/extract_private_gap_sites.py" \
        --input "${out}/graph/native.vcf.gz" --input "${sites}" \
        --input "${out}/bam/native.vcf.gz" --input "${bam}" \
        --output "${out}/private_gap_sites.vcf" --output "${out}/graph_gaps.bed" -- \
        python3 "${REPO_ROOT}/scripts/extract_private_gap_sites.py" \
        --graph-phased-vcf "${out}/graph/native.vcf.gz" --graph-sites "${sites}" \
        --linear-vcf "${out}/bam/native.vcf.gz" --contig "${graph_contig}" \
        --output "${out}/private_gap_sites.vcf" --gaps-bed "${out}/graph_gaps.bed" \
        --min-gq 10 --clean-snps-only --exclude-graph-positions \
        --min-vaf 0.30 --max-vaf 0.70 --bam "${bam}" \
        --min-bridge-reads 2 --min-mapq 20

    cached_step --state "${out}/hybrid/step.json" \
        --input "${REPO_ROOT}/pgphase" --input "${ref}" --input "${bam}" \
        --input "${sites}" --input "${gaf}" --input "${out}/private_gap_sites.vcf" \
        --output "${out}/hybrid/phased.bam" --output "${out}/hybrid/native.vcf" \
        --output "${out}/hybrid/candidates.tsv" -- \
        /usr/bin/time -v -o "${out}/hybrid/resources.txt" \
        "${REPO_ROOT}/pgphase" collect-hybrid-variation \
            --ref "${ref}" --bam "${bam}" --graph-sites "${sites}" --gaf "${gaf}" \
            --private-sites "${out}/private_gap_sites.vcf" -r "${bam_region}" \
            --min-read-margin 2 --min-phase-set-reads 50 -t "${THREADS}" \
            -o "${out}/hybrid/candidates.tsv" --phased-vcf-out "${out}/hybrid/native.vcf" \
            -b "${out}/hybrid/phased.bam"
    samtools index -@ "${THREADS}" "${out}/hybrid/phased.bam" 2>/dev/null || true

    cached_step --state "${out}/graph_lock/step.json" \
        --input "${REPO_ROOT}/scripts/merge_graph_hybrid_tags.py" \
        --input "${out}/graph/phased.bam" --input "${out}/hybrid/phased.bam" \
        --output "${out}/graph_lock/phased.bam" -- \
        python3 "${REPO_ROOT}/scripts/merge_graph_hybrid_tags.py" \
        --graph-bam "${out}/graph/phased.bam" --hybrid-bam "${out}/hybrid/phased.bam" \
        --output "${out}/graph_lock/phased.bam" --min-shared-reads 10 \
        --min-vote-margin 5 --min-purity 0.90 --require-both-haplotypes --threads 8
    samtools index -@ "${THREADS}" "${out}/graph_lock/phased.bam" 2>/dev/null || true

    cached_step --state "${out}/graph_bridge/step.json" \
        --input "${REPO_ROOT}/scripts/merge_graph_hybrid_tags.py" \
        --input "${out}/graph/phased.bam" --input "${out}/hybrid/phased.bam" \
        --output "${out}/graph_bridge/phased.bam" \
        --output "${out}/graph_bridge/phased.bam.bai" \
        --output "${out}/graph_bridge/bridge_edges.tsv" -- \
        /usr/bin/time -v -o "${out}/graph_bridge/resources.txt" \
        python3 "${REPO_ROOT}/scripts/merge_graph_hybrid_tags.py" \
        --graph-bam "${out}/graph/phased.bam" --hybrid-bam "${out}/hybrid/phased.bam" \
        --output "${out}/graph_bridge/phased.bam" --min-shared-reads 10 \
        --min-vote-margin 5 --min-purity 0.90 --require-both-haplotypes \
        --merge-graph-phase-sets --max-graph-bridge-distance 300000 \
        --graph-bridge-report "${out}/graph_bridge/bridge_edges.tsv" --threads 8

    cached_step --state "${out}/graph/inject_tags.step.json" \
        --input "${REPO_ROOT}/scripts/inject_hp_tags.py" \
        --input "${out}/graph/phased.bam" --input "${bam}" \
        --output "${out}/graph/tagged_surjected.bam" -- \
        python3 "${REPO_ROOT}/scripts/inject_hp_tags.py" \
            "${out}/graph/phased.bam" "${bam}" "${out}/graph/tagged_surjected.bam" samtools
    samtools index -@ "${THREADS}" "${out}/graph/tagged_surjected.bam"
    for label in graph bam hybrid graph_lock graph_bridge; do
        local phased_bam="${out}/${label}/phased.bam"
        local transfer_args=(
            --region "${chrom}" --min-reads 2 --min-ratio 0.70
            --support-cache "${out}/${label}/shared.support.tsv"
            --rebuild-support-cache
        )
        local transfer_outputs=(
            --output "${out}/${label}/shared.vcf"
            --output "${out}/${label}/shared.support.tsv"
        )
        if [[ "${label}" == "graph" ]]; then
            phased_bam="${out}/graph/tagged_surjected.bam"
        fi
        if [[ "${label}" == "graph_bridge" ]]; then
            transfer_args+=(
                --merge-phase-sets --merge-single-hap --merge-require-full-side
                --merge-min-sites 3 --merge-min-margin 2
                --merge-min-read-support 10 --merge-edge-order reads
                --merge-edge-report "${out}/graph_bridge/shared_merge_edges.tsv"
            )
            transfer_outputs+=(
                --output "${out}/graph_bridge/shared_merge_edges.tsv"
            )
        fi
        cached_step --state "${out}/${label}/shared_transfer.step.json" \
            --input "${REPO_ROOT}/scripts/phase_vcf_from_hp.py" \
            --input "${phased_bam}" --input "${shared_vcf}" \
            "${transfer_outputs[@]}" -- \
            python3 "${REPO_ROOT}/scripts/phase_vcf_from_hp.py" \
            "${phased_bam}" "${shared_vcf}" "${out}/${label}/shared.vcf" \
            "${transfer_args[@]}"
        compress_vcf "${out}/${label}/shared.vcf" "${out}/${label}/shared.vcf.gz"
    done

    if [[ ! -s "${out}/whatshap/phased.vcf.gz" ]]; then
        allow_competitor_run "${out}/whatshap/phased.vcf.gz"
        run_timed "${out}/whatshap/resources.txt" "${WHATSHAP}" phase --reference "${linear_ref}" \
            --ignore-read-groups --indels -o "${out}/whatshap/phased.vcf.gz" "${shared_vcf}" "${linear_bam}"
        "${TABIX}" -f -p vcf "${out}/whatshap/phased.vcf.gz"
    fi
    if [[ ! -s "${out}/whatshap_opt/phased.vcf.gz" ]]; then
        allow_competitor_run "${out}/whatshap_opt/phased.vcf.gz"
        run_timed "${out}/whatshap_opt/resources.txt" "${WHATSHAP}" phase --reference "${linear_ref}" \
            --ignore-read-groups --indels --distrust-genotypes \
            -o "${out}/whatshap_opt/phased.vcf.gz" "${shared_vcf}" "${linear_bam}"
        "${TABIX}" -f -p vcf "${out}/whatshap_opt/phased.vcf.gz"
    fi
    if [[ ! -s "${out}/hiphase/phased.vcf.gz" || ! -s "${out}/hiphase/phased.bam" ]]; then
        allow_competitor_run "${out}/hiphase/phased.vcf.gz and phased.bam"
        run_timed "${out}/hiphase/resources.txt" "${HIPHASE}" --reference "${linear_ref}" \
            --bam "${linear_bam}" --vcf "${shared_vcf}" --threads "${THREADS}" \
            --output-bam "${out}/hiphase/phased.bam" --output-vcf "${out}/hiphase/phased.vcf.gz" \
            --stats-file "${out}/hiphase/stats.native.tsv" \
            --blocks-file "${out}/hiphase/blocks.native.tsv" \
            --summary-file "${out}/hiphase/summary.native.tsv"
        "${TABIX}" -f -p vcf "${out}/hiphase/phased.vcf.gz"
    fi
    if [[ ! -s "${out}/longphase/phased.vcf.gz" ]]; then
        allow_competitor_run "${out}/longphase/phased.vcf.gz"
        if [[ ! -s "${out}/longphase/phased.vcf" ]]; then
            run_timed "${out}/longphase/resources.txt" "${LONGPHASE}" phase -s "${shared_vcf}" \
                -b "${linear_bam}" -r "${linear_ref}" -t "${THREADS}" --pb \
                -o "${out}/longphase/phased"
        fi
        compress_vcf "${out}/longphase/phased.vcf" "${out}/longphase/phased.vcf.gz"
    fi

    if [[ "${chrom}" == "chr20" ]]; then
        if [[ ! -s "${out}/longcalld/native.vcf" || \
              ! -s "${out}/longcalld/phased.bam" ]]; then
            allow_competitor_run "${out}/longcalld/native.vcf and phased.bam"
            [[ -x "${LONGCALLD}" ]] || die "missing LongcallD executable: ${LONGCALLD}"
            run_timed "${out}/longcalld/resources.txt" "${LONGCALLD}" call \
                --hifi -o "${out}/longcalld/native.vcf" \
                -b "${out}/longcalld/phased.bam" -t "${THREADS}" \
                "${ref}" "${bam}" "${graph_contig}"
        fi
        cached_step --state "${out}/longcalld/import.step.json" \
            --input "${out}/longcalld/native.vcf" \
            --input "$(dirname "${BASH_SOURCE[0]}")/longcalld_chr20_contigs.tsv" \
            --output "${out}/longcalld/phased.vcf.gz" -- \
            "${BCFTOOLS}" annotate \
                --rename-chrs "$(dirname "${BASH_SOURCE[0]}")/longcalld_chr20_contigs.tsv" \
                -Oz -o "${out}/longcalld/phased.vcf.gz" \
                "${out}/longcalld/native.vcf"
        "${TABIX}" -f -p vcf "${out}/longcalld/phased.vcf.gz"
    fi

    if [[ ! -s "${out}/whatshap/phased.bam" ]]; then
        allow_competitor_run "${out}/whatshap/phased.bam"
        "${WHATSHAP}" haplotag --reference "${linear_ref}" --output-threads "${THREADS}" \
            -o "${out}/whatshap/phased.bam" "${out}/whatshap/phased.vcf.gz" "${linear_bam}"
    fi
    if [[ ! -s "${out}/whatshap_opt/phased.bam" ]]; then
        allow_competitor_run "${out}/whatshap_opt/phased.bam"
        "${WHATSHAP}" haplotag --reference "${linear_ref}" --output-threads "${THREADS}" \
            -o "${out}/whatshap_opt/phased.bam" "${out}/whatshap_opt/phased.vcf.gz" "${linear_bam}"
    fi
    if [[ ! -s "${out}/longphase/phased.bam" ]]; then
        allow_competitor_run "${out}/longphase/phased.bam"
        "${LONGPHASE}" haplotag -s "${out}/longphase/phased.vcf.gz" -b "${linear_bam}" \
            -r "${linear_ref}" -t "${THREADS}" -o "${out}/longphase/phased"
    fi

    evaluate_reads "${out}/graph/phased.bam" "${truth_bam}" "${out}/eval/graph_reads"
    evaluate_reads "${out}/bam/phased.bam" "${truth_bam}" "${out}/eval/bam_reads"
    evaluate_reads "${out}/hybrid/phased.bam" "${truth_bam}" "${out}/eval/hybrid_reads"
    evaluate_reads "${out}/graph_lock/phased.bam" "${truth_bam}" "${out}/eval/graph_lock_reads"
    evaluate_reads "${out}/graph_bridge/phased.bam" "${truth_bam}" "${out}/eval/graph_bridge_reads"
    if [[ "${RUN_COMPETITORS}" == "1" ]]; then
        evaluate_reads "${out}/whatshap/phased.bam" "${truth_bam}" "${out}/eval/whatshap_reads"
        evaluate_reads "${out}/whatshap_opt/phased.bam" "${truth_bam}" "${out}/eval/whatshap_opt_reads"
        evaluate_reads "${out}/hiphase/phased.bam" "${truth_bam}" "${out}/eval/hiphase_reads"
        evaluate_reads "${out}/longphase/phased.bam" "${truth_bam}" "${out}/eval/longphase_reads"
    else
        for label in whatshap whatshap_opt hiphase longphase; do
            [[ -s "${out}/eval/${label}_reads/summary.json" ]] || \
                die "missing frozen read evaluation: ${out}/eval/${label}_reads/summary.json"
        done
    fi
    if [[ "${chrom}" == "chr20" ]]; then
        evaluate_reads "${out}/longcalld/phased.bam" "${truth_bam}" \
            "${out}/eval/longcalld_reads"
    fi

    evaluate_vcf "${chrom}" "${length}" graph "${out}/graph/shared.vcf.gz" "${out}"
    evaluate_vcf "${chrom}" "${length}" bam "${out}/bam/shared.vcf.gz" "${out}"
    evaluate_vcf "${chrom}" "${length}" hybrid "${out}/hybrid/shared.vcf.gz" "${out}"
    evaluate_vcf "${chrom}" "${length}" graph_lock "${out}/graph_lock/shared.vcf.gz" "${out}"
    evaluate_vcf "${chrom}" "${length}" graph_bridge "${out}/graph_bridge/shared.vcf.gz" "${out}"
    if [[ "${RUN_COMPETITORS}" == "1" ]]; then
        evaluate_vcf "${chrom}" "${length}" whatshap "${out}/whatshap/phased.vcf.gz" "${out}"
        evaluate_vcf "${chrom}" "${length}" whatshap_opt "${out}/whatshap_opt/phased.vcf.gz" "${out}"
        evaluate_vcf "${chrom}" "${length}" hiphase "${out}/hiphase/phased.vcf.gz" "${out}"
        evaluate_vcf "${chrom}" "${length}" longphase "${out}/longphase/phased.vcf.gz" "${out}"
    else
        for label in whatshap whatshap_opt hiphase longphase; do
            for suffix in compare.txt stats.tsv ngc50.json; do
                [[ -s "${out}/${label}.${suffix}" ]] || \
                    die "missing frozen VCF evaluation: ${out}/${label}.${suffix}"
            done
        done
    fi
    if [[ "${chrom}" == "chr20" ]]; then
        evaluate_vcf "${chrom}" "${length}" longcalld \
            "${out}/longcalld/phased.vcf.gz" "${out}"
    fi

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
        --recovery-vcf "graph_bridge=${out}/graph_bridge/shared.vcf.gz" \
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
