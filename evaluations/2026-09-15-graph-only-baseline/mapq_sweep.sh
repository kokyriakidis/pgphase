#!/usr/bin/env bash
# Graph-only whole-chr20 at several MAPQ floors. The surjected BAM's MAPQ is the
# graph aligner's own MAPQ (propagated by `pgbam annotate`), so this floor gates
# the same placement ambiguity on both channels.
set -euo pipefail
readonly REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
readonly OUT="${OUT:-/tmp/graph-only-mapq}"
readonly DATA_ROOT="${DATA_ROOT:-${HOME}/Downloads/pgphase-eval-data}"
readonly TRUTH="${DATA_ROOT}/truth/chr20/diplinator_merged.bam"
readonly THREADS="${THREADS:-8}"
mkdir -p "${OUT}"
for q in "$@"; do
    d="${OUT}/q${q}"; mkdir -p "${d}"
    echo "[$(date +%H:%M:%S)] graph-only at -q ${q}"
    "${REPO}/pgphase" collect-graph-variation \
        --ref "${REPO}/test_data/chm13v2.0.chr20.renamed.fa" \
        --sites "${REPO}/test_data/chr20.sites.striped.vcf.gz" \
        --gaf "${REPO}/test_data/HG002.chr20.annotated.coord.gaf.gz" \
        -r 'CHM13#0#chr20' -t "${THREADS}" -q "${q}" \
        -o "${d}/variants.tsv" --phased-vcf-out "${d}/native.vcf" \
        --phased-bam-out "${d}/phased.bam" \
        --phase-sites-out "${d}/phase_sites.tsv" \
        --filtered-sites-out "${d}/filtered_sites.tsv" \
        --phase-reads-out "${d}/phase_reads.tsv" \
        > "${d}/stdout.log" 2> "${d}/stderr.log"
    mkdir -p "${d}/eval"
    python3 "${REPO}/scripts/evaluate_phase_accuracy.py" "${d}/phased.bam" "${TRUTH}" \
        0 0 5 '' "${d}/eval" samtools '' '' '' '' > "${d}/eval.log" 2>&1
    python3 "${REPO}/evaluations/2026-09-15-graph-only-baseline/gap_inventory.py" \
        --phase-sites "${d}/phase_sites.tsv" --phase-reads "${d}/phase_reads.tsv" \
        --bam "${REPO}/test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam" \
        --output "${d}/gaps.tsv" > "${d}/gaps.log" 2>&1
done
echo "[$(date +%H:%M:%S)] done"
