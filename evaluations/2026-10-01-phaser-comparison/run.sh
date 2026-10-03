#!/usr/bin/env bash
set -euo pipefail

readonly REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
readonly PHASER_BIN="${PHASER_BIN:-${HOME}/micromamba/envs/bench-phasers/bin}"
readonly SHARED_ROOT="${SHARED_ROOT:-${HOME}/Downloads/pgphase-eval-data/shared_calls/chr20}"
readonly OUT_ROOT="${OUT_ROOT:-/tmp/pgphase-hiphase-comparison-2026-10-01-replay}"
readonly THREADS="${THREADS:-8}"

die() { echo "error: $*" >&2; exit 1; }
usage() {
    cat <<'EOF'
Usage: run.sh [--help]
Run fresh, sequential chr20 pgphase and HiPhase benchmarks.
Environment: PHASER_BIN, SHARED_ROOT, OUT_ROOT, THREADS (default 8).
OUT_ROOT must not exist; outputs are never reused as fresh timings.
Includes native DeepVariant and unphased pgphase-callset HiPhase arms.
EOF
}
if [[ "${1:-}" == --help ]]; then usage; exit 0; fi
[[ $# -eq 0 ]] || die "unexpected argument: $1"
[[ ! -e "${OUT_ROOT}" ]] || die "output already exists: ${OUT_ROOT}"
mkdir -p "${OUT_ROOT}"/{pgphase,hiphase_dv,hiphase_pg_calls}
cd "${REPO_ROOT}"

/usr/bin/time -v -o "${OUT_ROOT}/pgphase/resources.txt" \
    ./pgphase collect-graph-variation \
    --ref test_data/chm13v2.0.chr20.renamed.fa \
    --bam test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam \
    --sites test_data/chr20.sites.striped.vcf.gz \
    --gaf test_data/HG002.chr20.annotated.coord.gaf.gz \
    -t "${THREADS}" -r 'CHM13#0#chr20' \
    -o "${OUT_ROOT}/pgphase/candidates.tsv" \
    --phased-vcf-out "${OUT_ROOT}/pgphase/phased.vcf" \
    --phased-bam-out "${OUT_ROOT}/pgphase/phased.bam" \
    > "${OUT_ROOT}/pgphase/stdout.log" 2> "${OUT_ROOT}/pgphase/stderr.log"

/usr/bin/time -v -o "${OUT_ROOT}/hiphase_dv/resources.txt" \
    "${PHASER_BIN}/hiphase" --threads "${THREADS}" --ignore-read-groups \
    --reference "${SHARED_ROOT}/chm13v2.0.chr20.normalized.fa" \
    --bam "${SHARED_ROOT}/HG002.chr20.normalized.bam" \
    --vcf "${SHARED_ROOT}/deepvariant.vcf.gz" \
    --output-bam "${OUT_ROOT}/hiphase_dv/phased.bam" \
    --output-vcf "${OUT_ROOT}/hiphase_dv/phased.vcf.gz" \
    > "${OUT_ROOT}/hiphase_dv/stdout.log" 2> "${OUT_ROOT}/hiphase_dv/stderr.log"

"${PHASER_BIN}/python" evaluations/2026-10-01-phaser-comparison/prepare_calls.py \
    "${OUT_ROOT}/pgphase/phased.vcf" "${OUT_ROOT}/hiphase_pg_calls/input.vcf.gz"
/usr/bin/time -v -o "${OUT_ROOT}/hiphase_pg_calls/resources.txt" \
    "${PHASER_BIN}/hiphase" --threads "${THREADS}" --ignore-read-groups \
    --reference test_data/chm13v2.0.chr20.renamed.fa \
    --bam test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam \
    --vcf "${OUT_ROOT}/hiphase_pg_calls/input.vcf.gz" \
    --output-bam "${OUT_ROOT}/hiphase_pg_calls/phased.bam" \
    --output-vcf "${OUT_ROOT}/hiphase_pg_calls/phased.vcf.gz" \
    > "${OUT_ROOT}/hiphase_pg_calls/stdout.log" 2> "${OUT_ROOT}/hiphase_pg_calls/stderr.log"

/usr/bin/time -v -o "${OUT_ROOT}/pgphase/materialize_resources.txt" \
    "${PHASER_BIN}/python" evaluations/2026-10-01-phaser-comparison/materialize_tags.py \
    "${OUT_ROOT}/pgphase/phased.bam" \
    test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam \
    "${OUT_ROOT}/pgphase/haplotagged.bam" \
    > "${OUT_ROOT}/pgphase/materialize_stdout.log" \
    2> "${OUT_ROOT}/pgphase/materialize_stderr.log"

"${PHASER_BIN}/python" evaluations/2026-10-01-phaser-comparison/score.py "${OUT_ROOT}"
