#!/usr/bin/env bash
# Regression gates: collect-bam-variation and collect-graph-variation against fixture golden outputs.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PGPHASE="${PGPHASE:-$ROOT/pgphase}"
# Defaults: HG002 HiFi chr11 slice matching the committed golden outputs.
HIFI_FA="${HIFI_FA:-$ROOT/test_data/chr11_2M.fa}"
HIFI_BAM="${HIFI_BAM:-$ROOT/test_data/HG002_chr11_hifi_test.bam}"
ONT_FA="${ONT_FA:-$ROOT/test_data/chr11_2M.fa}"
ONT_BAM="${ONT_BAM:-$ROOT/test_data/HG002_chr11_ont_test.bam}"
HIFI_EXPECTED_TSV="${HIFI_EXPECTED_TSV:-$ROOT/test_data/expected/hifi_collect_expected.tsv}"
ONT_EXPECTED_TSV="${ONT_EXPECTED_TSV:-$ROOT/test_data/expected/ont_collect_expected.tsv}"
HIFI_EXPECTED_VCF="${HIFI_EXPECTED_VCF:-$ROOT/test_data/expected/hifi_collect_expected.phased.vcf}"
ONT_EXPECTED_VCF="${ONT_EXPECTED_VCF:-$ROOT/test_data/expected/ont_collect_expected.phased.vcf}"

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
if [[ ! -f "$HIFI_EXPECTED_TSV" ]]; then
  echo "error: HiFi expected TSV missing: $HIFI_EXPECTED_TSV" >&2
  exit 1
fi
if [[ ! -f "$HIFI_EXPECTED_VCF" ]]; then
  echo "error: HiFi expected phased VCF missing: $HIFI_EXPECTED_VCF" >&2
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
  local log
  log="$(mktemp)"
  set +e
  "$PGPHASE" collect-bam-variation -t "$threads" "$mode" --include-filtered \
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

t1_tsv="$TMPDIR/pgphase_validate_hifi_t1.tsv"
t1_vcf="$TMPDIR/pgphase_validate_hifi_t1.phased.vcf"
t4_tsv="$TMPDIR/pgphase_validate_hifi_t4.tsv"
t4_vcf="$TMPDIR/pgphase_validate_hifi_t4.phased.vcf"
count_t1="$(run_collect 1 "$t1_tsv" "$t1_vcf" "$HIFI_FA" "$HIFI_BAM" --hifi)"
count_t4="$(run_collect 4 "$t4_tsv" "$t4_vcf" "$HIFI_FA" "$HIFI_BAM" --hifi)"

hash_exp_hifi_tsv="$(hash_file "$HIFI_EXPECTED_TSV")"
hash_t1_tsv="$(hash_file "$t1_tsv")"
hash_t4_tsv="$(hash_file "$t4_tsv")"
hash_exp_hifi_vcf="$(hash_vcf_stable "$HIFI_EXPECTED_VCF")"
hash_t1_vcf="$(hash_vcf_stable "$t1_vcf")"
hash_t4_vcf="$(hash_vcf_stable "$t4_vcf")"

if [[ "$hash_t1_tsv" != "$hash_exp_hifi_tsv" ]]; then
  echo "error: HiFi TSV golden mismatch (t1 sha256 $hash_t1_tsv expected $hash_exp_hifi_tsv)" >&2
  echo "       expected: $HIFI_EXPECTED_TSV" >&2
  echo "       actual:   $t1_tsv" >&2
  exit 1
fi
if [[ "$hash_t1_tsv" != "$hash_t4_tsv" ]]; then
  echo "error: deterministic TSV mismatch (t1 sha256 $hash_t1_tsv vs t4 $hash_t4_tsv)" >&2
  exit 1
fi

if [[ "$hash_t1_vcf" != "$hash_exp_hifi_vcf" ]]; then
  echo "error: HiFi phased VCF golden mismatch (t1 sha256 $hash_t1_vcf expected $hash_exp_hifi_vcf)" >&2
  echo "       expected: $HIFI_EXPECTED_VCF" >&2
  echo "       actual:   $t1_vcf" >&2
  exit 1
fi
if [[ "$hash_t1_vcf" != "$hash_t4_vcf" ]]; then
  echo "error: deterministic phased VCF mismatch (t1 sha256 $hash_t1_vcf vs t4 $hash_t4_vcf)" >&2
  exit 1
fi

echo "ok: HiFi TSV golden match ($hash_t1_tsv), t1 vs t4 deterministic; sites=$count_t1"
echo "ok: HiFi phased VCF golden match ($hash_t1_vcf), t1 vs t4 deterministic"

if [[ -f "$ONT_FA" && -f "$ONT_BAM" ]]; then
  ont_tsv="$TMPDIR/pgphase_validate_ont.tsv"
  ont_vcf="$TMPDIR/pgphase_validate_ont.phased.vcf"
  count_ont="$(run_collect 1 "$ont_tsv" "$ont_vcf" "$ONT_FA" "$ONT_BAM" --ont)"
  if [[ -f "$ONT_EXPECTED_TSV" ]]; then
    hash_exp_ont_tsv="$(hash_file "$ONT_EXPECTED_TSV")"
    hash_ont_tsv="$(hash_file "$ont_tsv")"
    if [[ "$hash_ont_tsv" != "$hash_exp_ont_tsv" ]]; then
      echo "error: ONT TSV golden mismatch (sha256 $hash_ont_tsv expected $hash_exp_ont_tsv)" >&2
      echo "       expected: $ONT_EXPECTED_TSV" >&2
      echo "       actual:   $ont_tsv" >&2
      exit 1
    fi
    echo "ok: ONT TSV golden match ($hash_ont_tsv); sites=$count_ont"
  else
    echo "ok: ONT fixture ran ($count_ont sites) -> $ont_tsv"
    echo "warn: ONT expected TSV not found: $ONT_EXPECTED_TSV (golden check skipped)"
  fi
  if [[ -f "$ONT_EXPECTED_VCF" ]]; then
    hash_exp_ont_vcf="$(hash_vcf_stable "$ONT_EXPECTED_VCF")"
    hash_ont_vcf="$(hash_vcf_stable "$ont_vcf")"
    if [[ "$hash_ont_vcf" != "$hash_exp_ont_vcf" ]]; then
      echo "error: ONT phased VCF golden mismatch (sha256 $hash_ont_vcf expected $hash_exp_ont_vcf)" >&2
      echo "       expected: $ONT_EXPECTED_VCF" >&2
      echo "       actual:   $ont_vcf" >&2
      exit 1
    fi
    echo "ok: ONT phased VCF golden match ($hash_ont_vcf)"
  else
    echo "warn: ONT expected phased VCF not found: $ONT_EXPECTED_VCF (golden check skipped)"
  fi
else
  echo "skip: ONT fixture not present (FA=$ONT_FA BAM=$ONT_BAM)"
fi

# --- collect-graph-variation ---------------------------------------------------
# Graph-only on the MICB-KIR3DL1 graph slice, and union gap phasing (--bam) on
# the chr20:25.0-25.5 Mb slice. Both are compared with committed goldens, the
# union run at 1 and 4 threads.
GRAPH_DIR="${GRAPH_DIR:-$ROOT/test_data/graph}"
CHR20_DIR="${CHR20_DIR:-$ROOT/test_data/graph_chr20}"
GRAPH_EXPECTED_TSV="${GRAPH_EXPECTED_TSV:-$ROOT/test_data/expected/graph_micb_expected.tsv}"
UNION_EXPECTED_TSV="${UNION_EXPECTED_TSV:-$ROOT/test_data/expected/union_chr20_expected.tsv}"
UNION_EXPECTED_VCF="${UNION_EXPECTED_VCF:-$ROOT/test_data/expected/union_chr20_expected.phased.vcf}"

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

check_golden() {
  local what="$1" actual="$2" expected="$3" hasher="$4"
  local a e
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

for t in 1 4; do
  count_u="$(run_graph -t "$t" --ref "$CHR20_DIR/ref.fa.gz" --sites "$CHR20_DIR/chr20_25M.sites.vcf.gz" \
    --gaf "$CHR20_DIR/HG002_chr20_25M.coord.gaf.gz" --bam "$CHR20_DIR/HG002_chr20_25M.bam" \
    -o "$TMPDIR/pgphase_validate_union_t$t.tsv" \
    --phased-vcf-out "$TMPDIR/pgphase_validate_union_t$t.phased.vcf")"
  check_golden "union TSV t$t (sites=$count_u)" "$TMPDIR/pgphase_validate_union_t$t.tsv" "$UNION_EXPECTED_TSV" hash_file
  check_golden "union phased VCF t$t" "$TMPDIR/pgphase_validate_union_t$t.phased.vcf" "$UNION_EXPECTED_VCF" hash_vcf_stable
done
