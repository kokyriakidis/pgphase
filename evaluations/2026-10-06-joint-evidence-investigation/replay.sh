#!/usr/bin/env bash
set -euo pipefail

if [[ "${1:-}" == "--help" ]]; then
    echo "Usage: BENCH_PYTHON=python-with-pysam HIPHASE=hiphase replay.sh"
    echo "Replay three current owning chunks, six same-read HiPhase controls, and diagnostic audits."
    exit 0
fi
die() { echo "$*" >&2; exit 1; }
[[ $# -eq 0 ]] || die "Unknown argument: $1 (use --help)"
readonly repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
readonly bench_python="${BENCH_PYTHON:-python3}"
readonly hiphase="${HIPHASE:-hiphase}"
readonly evaluation="evaluations/2026-10-06-joint-evidence-investigation"
readonly replays="test_data/tmp_joint_evidence_investigation"
cd "${repo_root}"
for chunk in 4 5 6; do
    directory="${replays}/${chunk}"
    mkdir -p "${directory}"
    LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu ./pgphase collect-graph-variation \
        --ref test_data/chm13v2.0.chr20.renamed.fa \
        --bam test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam \
        --sites test_data/chr20.sites.striped.vcf.gz \
        --gaf test_data/HG002.chr20.annotated.coord.gaf.gz \
        -t 1 -r "CHM13#0#chr20:$((chunk*1000000+1))-$(((chunk+1)*1000000))" \
        -o "${directory}/candidates.tsv" \
        --phased-vcf-out "${directory}/phased.vcf" \
        --phased-bam-out "${directory}/phased.bam" \
        --phase-matrix-dump "${directory}/matrix" \
        --recovery-audit-out "${directory}/recovery.tsv" \
        > "${directory}/stdout.log" 2> "${directory}/stderr.log"
done
"${bench_python}" "${evaluation}/prepare_competitor.py"
for mode in default local dv add-marker diploid both; do
    input="${replays}/hiphase5/${mode}-input.vcf.gz"
    extra=()
    if [[ "${mode}" == default || "${mode}" == local ]]; then
        input="${replays}/hiphase5/input.vcf.gz"
    fi
    if [[ "${mode}" == local ]]; then extra=(--disable-global-realignment); fi
    "${hiphase}" --bam "${replays}/hiphase5/input.bam" --vcf "${input}" \
        --output-bam "${replays}/hiphase5/${mode}.bam" \
        --output-vcf "${replays}/hiphase5/${mode}.vcf.gz" \
        --reference test_data/chm13v2.0.chr20.renamed.fa \
        --ignore-read-groups -t 1 -vvv "${extra[@]}" \
        > "${replays}/hiphase5/${mode}.log" 2>&1
done
"${bench_python}" "${evaluation}/audit.py" --contract "${evaluation}/native-gap-contract.tsv"
"${bench_python}" "${evaluation}/audit_matrices.py"
"${bench_python}" "${evaluation}/replay_joint.py"
"${bench_python}" "${evaluation}/replay_joint.py" --snp-priority
"${bench_python}" "${evaluation}/audit_competitor.py"
