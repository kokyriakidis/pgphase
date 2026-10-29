#!/usr/bin/env bash
# Regression gates: collect-bam-variation and collect-graph-variation against fixture golden outputs, once per
# noisy-region MSA this pgphase has: plegma (the default; goldens in test_data/expected/plegma) and abPOA (--msa
# abpoa; test_data/expected, longcallD's outputs: the port's parity check). UPDATE_GOLDEN=plegma rewrites plegma's
# goldens from this build's outputs.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PGPHASE="${PGPHASE:-$ROOT/pgphase}"
# Defaults: HG002 HiFi chr11 slice matching the committed golden outputs.
HIFI_FA="${HIFI_FA:-$ROOT/test_data/chr11_2M.fa}"
HIFI_BAM="${HIFI_BAM:-$ROOT/test_data/HG002_chr11_hifi_test.bam}"
ONT_FA="${ONT_FA:-$ROOT/test_data/chr11_2M.fa}"
ONT_BAM="${ONT_BAM:-$ROOT/test_data/HG002_chr11_ont_test.bam}"

# Backward-compatible fallback for older mixed-case fixture paths.
if [[ ! -f "$ONT_FA" && -f "$ROOT/test_data/longCallD/chr11_2M.fa" ]]; then
  ONT_FA="$ROOT/test_data/longCallD/chr11_2M.fa"
fi
if [[ ! -f "$ONT_BAM" && -f "$ROOT/test_data/longCallD/HG002_chr11_ont_test.bam" ]]; then
  ONT_BAM="$ROOT/test_data/longCallD/HG002_chr11_ont_test.bam"
fi
TMPDIR="${TMPDIR:-/tmp}"

if [[ ! -x "$PGPHASE" ]]; then
  echo "error: $PGPHASE not found or not executable (build with: make)" >&2
  exit 1
fi

if [[ ! -f "$HIFI_FA" || ! -f "$HIFI_BAM" ]]; then
  echo "error: HiFi fixture missing (FA=$HIFI_FA BAM=$HIFI_BAM)" >&2
  exit 1
fi

hash_file() {
  local f="$1"
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum < "$f" | awk '{print $1}'
  elif command -v shasum >/dev/null 2>&1; then
    shasum -a 256 "$f" | awk '{print $1}'
  else
    echo "error: no sha256 tool found (need sha256sum or shasum)" >&2
    exit 1
  fi
}

# Phased VCF headers include ##fileDate=YYYYMMDD from run day; strip for stable golden compares.
hash_vcf_stable() {
  local f="$1"
  if command -v sha256sum >/dev/null 2>&1; then
    grep -v '^##fileDate=' "$f" | sha256sum | awk '{print $1}'
  elif command -v shasum >/dev/null 2>&1; then
    grep -v '^##fileDate=' "$f" | shasum -a 256 | awk '{print $1}'
  else
    echo "error: no sha256 tool found (need sha256sum or shasum)" >&2
    exit 1
  fi
}

run_collect() {
  local threads="$1"
  local out_tsv="$2"
  local out_vcf="$3"
  local fa="$4"
  local bam="$5"
  local mode="$6"
  local msa="$7"
  local log
  log="$(mktemp)"
  set +e
  "$PGPHASE" collect-bam-variation -t "$threads" "$mode" --include-filtered --msa "$msa" \
    --ref "$fa" --bam "$bam" -o "$out_tsv" --phased-vcf-out "$out_vcf" >"$log" 2>&1
  local st=$?
  set -e
  if [[ "$st" -ne 0 ]]; then
    cat "$log" >&2
    rm -f "$log"
    exit "$st"
  fi
  grep -oE 'Collected [0-9]+ candidate variant sites' "$log" | grep -oE '[0-9]+' | head -1
  rm -f "$log"
}

run_graph() {
  local log
  log="$(mktemp)"
  set +e
  "$PGPHASE" collect-graph-variation "$@" >"$log" 2>&1
  local st=$?
  set -e
  if [[ "$st" -ne 0 ]]; then
    cat "$log" >&2
    rm -f "$log"
    exit "$st"
  fi
  grep -oE 'Collected [0-9]+ candidate variant sites' "$log" | grep -oE '[0-9]+' | head -1
  rm -f "$log"
}

# A bgzipped file by its decompressed content.
hash_gz() {
  gzip -dc "$1" | if command -v sha256sum >/dev/null 2>&1; then sha256sum; else shasum -a 256; fi | awk '{print $1}'
}

# what, actual, expected, hasher[, msa]: with UPDATE_GOLDEN set to msa, expected becomes actual.
check_golden() {
  local what="$1" actual="$2" expected="$3" hasher="$4" msa="${5:-}"
  local a e
  if [[ -n "$msa" && "${UPDATE_GOLDEN:-}" == "$msa" ]]; then
    mkdir -p "$(dirname "$expected")"
    cp "$actual" "$expected"
    echo "updated: $what golden ($expected)"
    return
  fi
  if [[ ! -f "$expected" ]]; then
    echo "error: $what golden missing: $expected" >&2
    exit 1
  fi
  a="$("$hasher" "$actual")"
  e="$("$hasher" "$expected")"
  if [[ "$a" != "$e" ]]; then
    echo "error: $what golden mismatch (sha256 $a expected $e)" >&2
    echo "       expected: $expected" >&2
    echo "       actual:   $actual" >&2
    exit 1
  fi
  echo "ok: $what golden match ($a)"
}

same_output() {
  local what="$1" first="$2" second="$3" hasher="$4"
  if [[ "$("$hasher" "$first")" != "$("$hasher" "$second")" ]]; then
    echo "error: $what not deterministic ($first vs $second)" >&2
    exit 1
  fi
}

# collect-bam-variation on the chr11 HiFi slice (1 and 4 threads) and the ONT one, with one MSA.
bam_gates() {
  local msa="$1" dir="$2"
  local t1_tsv="$TMPDIR/pgphase_validate_${msa}_hifi_t1.tsv" t1_vcf="$TMPDIR/pgphase_validate_${msa}_hifi_t1.phased.vcf"
  local t4_tsv="$TMPDIR/pgphase_validate_${msa}_hifi_t4.tsv" t4_vcf="$TMPDIR/pgphase_validate_${msa}_hifi_t4.phased.vcf"
  local count_t1
  count_t1="$(run_collect 1 "$t1_tsv" "$t1_vcf" "$HIFI_FA" "$HIFI_BAM" --hifi "$msa")"
  run_collect 4 "$t4_tsv" "$t4_vcf" "$HIFI_FA" "$HIFI_BAM" --hifi "$msa" > /dev/null
  check_golden "$msa HiFi TSV (sites=$count_t1)" "$t1_tsv" "$dir/hifi_collect_expected.tsv" hash_file "$msa"
  check_golden "$msa HiFi phased VCF" "$t1_vcf" "$dir/hifi_collect_expected.phased.vcf" hash_vcf_stable "$msa"
  same_output "$msa HiFi TSV, 1 and 4 threads," "$t1_tsv" "$t4_tsv" hash_file
  same_output "$msa HiFi phased VCF, 1 and 4 threads," "$t1_vcf" "$t4_vcf" hash_vcf_stable
  if [[ -f "$ONT_FA" && -f "$ONT_BAM" ]]; then
    local ont_tsv="$TMPDIR/pgphase_validate_${msa}_ont.tsv" ont_vcf="$TMPDIR/pgphase_validate_${msa}_ont.phased.vcf"
    local count_ont
    count_ont="$(run_collect 1 "$ont_tsv" "$ont_vcf" "$ONT_FA" "$ONT_BAM" --ont "$msa")"
    check_golden "$msa ONT TSV (sites=$count_ont)" "$ont_tsv" "$dir/ont_collect_expected.tsv" hash_file "$msa"
    check_golden "$msa ONT phased VCF" "$ont_vcf" "$dir/ont_collect_expected.phased.vcf" hash_vcf_stable "$msa"
  else
    echo "skip: ONT fixture not present (FA=$ONT_FA BAM=$ONT_BAM)"
  fi
}

# The truth-free contracts of the union run's phased VCF (scripts/check_contracts.py): VCF order and REF bases must
# hold; overlapping records on one haplotype and MSA regions the records do not spell must be among the accepted ones
# ($dir/union_chr20_contracts.tsv, rewritten with the goldens by UPDATE_GOLDEN).
union_contracts() {
  local msa="$1" dir="$2" diag="$3"
  local args=(--vcf "$TMPDIR/pgphase_validate_${msa}_union_t1.phased.vcf" --ref "$CHR20_DIR/ref.fa.gz"
    --msa-dump "$diag/msa_dump.tsv" --trace "$diag/trace.tsv")
  if [[ "${UPDATE_GOLDEN:-}" == "$msa" ]]; then
    python3 "$ROOT/scripts/check_contracts.py" "${args[@]}" --write-known "$dir/union_chr20_contracts.tsv"
    echo "updated: $msa union contracts ($dir/union_chr20_contracts.tsv)"
    return
  fi
  if ! python3 "$ROOT/scripts/check_contracts.py" "${args[@]}" --known "$dir/union_chr20_contracts.tsv"; then
    echo "error: $msa union contracts: a violation not in $dir/union_chr20_contracts.tsv" >&2
    exit 1
  fi
  echo "ok: $msa union contracts"
}

# Union gap phasing (collect-graph-variation --bam) on the chr20:25.0-25.5 Mb slice, 1 and 4 threads, with one MSA.
# The 1-thread run also writes the aligned output with --refine-aln, the 4-thread run without: the TSV and VCF must not
# change, and the two alignments must hold the same records and labels (scripts/check_aligned_output.py).
union_gates() {
  local msa="$1" dir="$2" t count_u
  local diag="$TMPDIR/pgphase_validate_${msa}_union_diag"
  rm -rf "$diag"
  mkdir -p "$diag"
  for t in 1 4; do
    # The 1-thread run also writes its stage trace and MSA dump, for the contracts below.
    local refine=()
    [[ "$t" == 1 ]] && refine=(--refine-aln)
    count_u="$(if [[ "$t" == 1 ]]; then export PGPHASE_TRACE="$diag" PGPHASE_MSA_DUMP="$diag/msa_dump.tsv"; fi
      run_graph -t "$t" --msa "$msa" --ref "$CHR20_DIR/ref.fa.gz" --sites "$CHR20_DIR/chr20_25M.sites.vcf.gz" \
      --gaf "$CHR20_DIR/HG002_chr20_25M.coord.gaf.gz" --bam "$CHR20_DIR/HG002_chr20_25M.bam" \
      -o "$TMPDIR/pgphase_validate_${msa}_union_t$t.tsv" \
      --phased-vcf-out "$TMPDIR/pgphase_validate_${msa}_union_t$t.phased.vcf" \
      --out-sam "$TMPDIR/pgphase_validate_${msa}_union_t$t.sam" "${refine[@]}")"
    if [[ "$t" == 1 ]]; then
      check_golden "$msa union TSV (sites=$count_u)" "$TMPDIR/pgphase_validate_${msa}_union_t1.tsv" \
        "$dir/union_chr20_expected.tsv" hash_file "$msa"
      check_golden "$msa union phased VCF" "$TMPDIR/pgphase_validate_${msa}_union_t1.phased.vcf" \
        "$dir/union_chr20_expected.phased.vcf" hash_vcf_stable "$msa"
      union_contracts "$msa" "$dir" "$diag"
      # The trace's own contract: columns, documented reasons, parseable details (scripts/check_trace.py).
      if ! python3 "$ROOT/scripts/check_trace.py" --trace "$diag" --doc "$ROOT/docs/TRACE.md"; then
        echo "error: $msa union trace: see scripts/check_trace.py and docs/TRACE.md" >&2
        exit 1
      fi
      echo "ok: $msa union trace"
    fi
  done
  same_output "$msa union TSV, 1 and 4 threads," "$TMPDIR/pgphase_validate_${msa}_union_t1.tsv" \
    "$TMPDIR/pgphase_validate_${msa}_union_t4.tsv" hash_file
  same_output "$msa union phased VCF, 1 and 4 threads," "$TMPDIR/pgphase_validate_${msa}_union_t1.phased.vcf" \
    "$TMPDIR/pgphase_validate_${msa}_union_t4.phased.vcf" hash_vcf_stable
  if ! python3 "$ROOT/scripts/check_aligned_output.py" --refined "$TMPDIR/pgphase_validate_${msa}_union_t1.sam" \
      --plain "$TMPDIR/pgphase_validate_${msa}_union_t4.sam" --ref "$CHR20_DIR/ref.fa.gz"; then
    echo "error: $msa union aligned output (--refine-aln): see scripts/check_aligned_output.py" >&2
    exit 1
  fi
  echo "ok: $msa union aligned output"
}

EXPECTED_DIR="${EXPECTED_DIR:-$ROOT/test_data/expected}"
MSAS=(abpoa)
if "$PGPHASE" --version | grep -q "MSA: plegma"; then
  MSAS=(plegma abpoa)
fi
for msa in "${MSAS[@]}"; do
  dir="$EXPECTED_DIR"
  [[ "$msa" == plegma ]] && dir="$EXPECTED_DIR/plegma"
  bam_gates "$msa" "$dir"
done

# --- collect-graph-variation ---------------------------------------------------
# Graph-only on the MICB-KIR3DL1 graph slice, and union gap phasing (--bam) on
# the chr20:25.0-25.5 Mb slice. Both are compared with committed goldens, the
# union run at 1 and 4 threads.
GRAPH_DIR="${GRAPH_DIR:-$ROOT/test_data/graph}"
CHR20_DIR="${CHR20_DIR:-$ROOT/test_data/graph_chr20}"
GRAPH_EXPECTED_TSV="${GRAPH_EXPECTED_TSV:-$ROOT/test_data/expected/graph_micb_expected.tsv}"

g_tsv="$TMPDIR/pgphase_validate_graph.tsv"
count_g="$(run_graph -t 2 --ref "$GRAPH_DIR/ref.fa.gz" --sites "$GRAPH_DIR/micb-kir3dl1.sites.vcf.gz" \
  --gaf "$GRAPH_DIR/micb-kir3dl1_HG003.coord.gaf.gz" -o "$g_tsv")"
check_golden "graph-only TSV (sites=$count_g)" "$g_tsv" "$GRAPH_EXPECTED_TSV" hash_file

# The same run on the .pgsites catalog folder pgcat build writes.
GRAPH_CATALOG="$GRAPH_DIR/micb-kir3dl1.catalog"
count_g="$(run_graph -t 2 --ref "$GRAPH_DIR/ref.fa.gz" --sites "$GRAPH_CATALOG" \
  --gaf "$GRAPH_DIR/micb-kir3dl1_HG003.coord.gaf.gz" -o "$g_tsv")"
check_golden "graph-only TSV on .pgsites (sites=$count_g)" "$g_tsv" "$GRAPH_EXPECTED_TSV" hash_file

# pgcat rebuilds that catalog from its committed snarls (so without vg).
PGCAT="${PGCAT:-$ROOT/pgcat}"
if [[ -x "$PGCAT" ]]; then
  cat_out="$TMPDIR/pgphase_validate_catalog"
  rm -rf "$cat_out"
  mkdir -p "$cat_out"
  cp -r "$GRAPH_CATALOG/snarls" "$cat_out/snarls"
  "$PGCAT" build --graph "$GRAPH_DIR/micb-kir3dl1.gbz" --reference-sample CHM13 --out "$cat_out" \
    --region chr6 --region chr19 --threads 2 2>/dev/null
  for c in chr6 chr19; do
    check_golden "pgcat $c catalog" "$cat_out/$c.pgsites" "$GRAPH_CATALOG/$c.pgsites" hash_file
  done
  # ... and the reference-indexed GAF from the raw alignments: the same
  # records and index (compared decompressed).
  "$PGCAT" index-gaf --graph "$GRAPH_DIR/micb-kir3dl1.gbz" --gaf "$GRAPH_DIR/micb-kir3dl1_HG003.gaf" \
    --reference-sample CHM13 --out "$cat_out/reads.coord.gaf.gz" -t 2 2>/dev/null
  check_golden "pgcat index-gaf records" "$cat_out/reads.coord.gaf.gz" "$GRAPH_DIR/micb-kir3dl1_HG003.coord.gaf.gz" hash_gz
  check_golden "pgcat index-gaf index" "$cat_out/reads.coord.gaf.gz.tbi" "$GRAPH_DIR/micb-kir3dl1_HG003.coord.gaf.gz.tbi" hash_gz
else
  echo "skip: pgcat not built (make pgcat); catalog and GAF rebuilds not checked"
fi

for msa in "${MSAS[@]}"; do
  dir="$EXPECTED_DIR"
  [[ "$msa" == plegma ]] && dir="$EXPECTED_DIR/plegma"
  union_gates "$msa" "$dir"
done
