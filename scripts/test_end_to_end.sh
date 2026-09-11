#!/usr/bin/env bash
#
# End-to-end smoke test of the full evaluation chain, using only pinned tools:
#
#   reads  -> minimap2 (third_party/minimap2, pinned to a release tag)
#          -> hiphap   (third_party/hiphap, pinned commit + vendored Cargo.lock)
#          -> truth BAM
#   GAF    -> pgphase collect-graph-variation -> phased BAM
#   both   -> evaluate_phase_accuracy.py -> accuracy summary
#
# This is the test that would have caught, in order: the hiphap rename flipping
# merged output to the default (-p now required), the rust-htslib/hts-sys build
# break, the spurious paftools.js requirement, and hiphap's -A auto-estimate
# failing on small inputs.  It runs in about a minute on the checked-in
# 500 kb chr20 fixture.
#
# The only input not in the repo is the HG002 v1.1 diploid assembly, which the
# truth BAM is built against.  Point --asm at it; the test skips cleanly if it
# is absent so CI without the assembly still passes the rest.
#
# Usage:
#   ./scripts/test_end_to_end.sh [--asm /path/to/hg002v1.1.fasta] [--outdir DIR]
#                                [--eval-data DIR]   # default $PGPHASE_EVAL_DATA
#
# Requires: samtools, and `make eval-tools` already run (or hiphap/minimap2 on PATH).

set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ASM="${HG002_ASM:-/home/kokyriakidis/Downloads/hg002v1.1.fasta}"
# Reusable eval artifacts live OUTSIDE the repo so re-cloning, updating or
# git-cleaning pgphase cannot delete them.  A truth BAM costs 35-90 min to
# rebuild and depends only on (reads, assembly, minimap2, hiphap) -- never on a
# pgphase setting -- so it is cached here and shared by every run.
EVAL_DATA="${PGPHASE_EVAL_DATA:-/home/kokyriakidis/Downloads/pgphase-eval-data}"
OUTDIR=""
THREADS=8

while [[ $# -gt 0 ]]; do
    case "$1" in
        --asm)       ASM="$2";       shift 2 ;;
        --eval-data) EVAL_DATA="$2"; shift 2 ;;
        --outdir)  OUTDIR="$2";  shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        -h|--help) sed -n '2,25p' "${BASH_SOURCE[0]}"; exit 0 ;;
        *) echo "unknown option: $1" >&2; exit 1 ;;
    esac
done

FIXTURE="$REPO/test_data/graph_chr20"
PGPHASE="$REPO/pgphase"
status=0
fail() { echo "  FAIL: $*"; status=1; }
ok()   { echo "  ok:   $*"; }

echo "== preflight =="
[[ -x "$PGPHASE" ]] || { echo "pgphase not built; run make all"; exit 1; }
for f in ref.fa.gz chr20_25M.sites.vcf.gz HG002_chr20_25M.coord.gaf.gz HG002_chr20_25M.bam; do
    [[ -f "$FIXTURE/$f" ]] || { echo "missing fixture: $FIXTURE/$f"; exit 1; }
done
command -v samtools >/dev/null || { echo "samtools not on PATH"; exit 1; }

MM2="$REPO/third_party/minimap2/minimap2"
HIPHAP="$REPO/third_party/hiphap/target/release/hiphap"
[[ -x "$MM2" ]]    || MM2="$(command -v minimap2 || true)"
[[ -x "$HIPHAP" ]] || HIPHAP="$(command -v hiphap || true)"
[[ -n "$MM2"    && -x "$MM2"    ]] && ok "minimap2 $("$MM2" --version)" || fail "minimap2 missing (make eval-tools)"
[[ -n "$HIPHAP" && -x "$HIPHAP" ]] && ok "$("$HIPHAP" --version)"       || fail "hiphap missing (make eval-tools)"
[[ "$status" -eq 0 ]] || exit 1

if [[ ! -f "$ASM" ]]; then
    echo
    echo "SKIP: diploid assembly not found at $ASM"
    echo "      Pass --asm /path/to/hg002v1.1.fasta to run the truth-BAM half."
    echo "      Tooling preflight passed."
    exit 0
fi

WORK="${OUTDIR:-$(mktemp -d -t pgphase_e2e.XXXXXX)}"
mkdir -p "$WORK"
[[ -n "$OUTDIR" ]] || trap 'rm -rf "$WORK"' EXIT
echo "  work dir: $WORK"

echo
echo "== 1/4 extracting haplotype references and reads =="
[[ -f "$ASM.fai" ]] || samtools faidx "$ASM"
samtools faidx "$ASM" chr20_MATERNAL > "$WORK/mat.fa"
samtools faidx "$ASM" chr20_PATERNAL > "$WORK/pat.fa"
samtools fastq -@ "$THREADS" "$FIXTURE/HG002_chr20_25M.bam" > "$WORK/reads.fq" 2>/dev/null
n_reads=$(( $(wc -l < "$WORK/reads.fq") / 4 ))
[[ "$n_reads" -gt 0 ]] && ok "$n_reads reads extracted" || fail "no reads extracted"

echo
echo "== 2/4 building truth BAM (minimap2 + hiphap) =="
# Cached in the external store: the truth BAM is a function of the reads and the
# assembly only, so it is reused across runs rather than rebuilt each time.
TRUTH_DIR="$EVAL_DATA/truth/chr20_25M"
TRUTH="$TRUTH_DIR/diplinator_merged.bam"
if [[ -s "$TRUTH" ]]; then
    ok "reusing cached truth BAM: $TRUTH"
else
    mkdir -p "$TRUTH_DIR"
    "$REPO/scripts/build_truth_bam.sh" \
        --reads "$WORK/reads.fq" --mat-ref "$WORK/mat.fa" --pat-ref "$WORK/pat.fa" \
        --outdir "$TRUTH_DIR" --threads "$THREADS" > "$WORK/truth.log" 2>&1 || {
            fail "build_truth_bam.sh failed; see $WORK/truth.log"; tail -5 "$WORK/truth.log"; exit 1; }
    # Keep only the merged BAM; the per-haplotype SAMs are large and derivable.
    rm -f "$TRUTH_DIR"/vs_*.sam "$TRUTH_DIR"/hiphap_*.sam "$TRUTH_DIR"/*_span_chrom.fastq
fi
[[ -s "$TRUTH" ]] && ok "truth BAM: $(samtools view -c "$TRUTH") reads" || fail "no truth BAM produced"
# Both haplotype tags must be present, or the haplotype-of-origin stamping broke
# -- which is exactly what a silent revert to hiphap's merged default would do.
for tag in MAT PAT; do
    c=$(samtools view "$TRUTH" 2>/dev/null | grep -c "HO:Z:$tag" || true)
    [[ "$c" -gt 0 ]] && ok "HO:Z:$tag present ($c reads)" || fail "no HO:Z:$tag reads -- hiphap -p missing?"
done

echo
echo "== 3/4 phasing the same fixture =="
"$PGPHASE" collect-graph-variation \
    --ref "$FIXTURE/ref.fa.gz" --sites "$FIXTURE/chr20_25M.sites.vcf.gz" \
    --gaf "$FIXTURE/HG002_chr20_25M.coord.gaf.gz" \
    -o "$WORK/cands.tsv" --phased-bam-out "$WORK/phased.bam" \
    -t "$THREADS" > "$WORK/phase.log" 2>&1 || { fail "collect-graph-variation failed"; exit 1; }
n_hp=$(samtools view "$WORK/phased.bam" | grep -c "HP:i:" || true)
[[ "$n_hp" -gt 0 ]] && ok "$n_hp reads carry HP tags" || fail "no HP tags in phased BAM"

echo
echo "== 4/4 evaluating =="
mkdir -p "$WORK/eval"
python3 "$REPO/scripts/evaluate_phase_accuracy.py" \
    "$WORK/phased.bam" "$TRUTH" 0 0 5 "" "$WORK/eval" samtools "" "" "" "" \
    > "$WORK/eval.log" 2>&1 || { fail "evaluate_phase_accuracy.py failed"; tail -5 "$WORK/eval.log"; exit 1; }

python3 - "$WORK/eval/summary.json" <<'PY'
import json, sys
d = json.load(open(sys.argv[1]))
ev, disc = d["total_reads_evaluated"], d["discordant_reads"]
print(f"  reads evaluated : {ev:,}")
print(f"  concordant      : {d['concordant_reads']:,}")
print(f"  discordant      : {disc}")
print(f"  hamming         : {d['hamming_error_rate']:.6f}")
print(f"  phase sets      : {d['phase_sets_evaluated']}")
# This fixture is a clean, well-covered 500 kb slice; anything but near-perfect
# phasing means a real regression somewhere in the chain, not a hard region.
sys.exit(0 if ev > 500 and disc / max(ev, 1) < 0.02 else 1)
PY
[[ $? -eq 0 ]] && ok "accuracy within expected range for this fixture" \
               || fail "accuracy outside expected range"

echo
if [[ "$status" -eq 0 ]]; then echo "PASS: end-to-end chain works with pinned tools"; else echo "FAILURES above"; fi
exit "$status"
