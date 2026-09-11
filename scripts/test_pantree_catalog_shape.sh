#!/usr/bin/env bash
#
# Integration test: does collect-graph-variation ingest a pantree-shaped catalog?
#
# pantree_to_pgphase_catalog.py emits a catalog whose *shape* differs from a
# vg deconstruct one in three ways that matter to src/graph_sites.cpp:
#   1. every record is biallelic (one variant edge = one binary site),
#   2. there is no LV/PS/PA nesting, so no parent gating,
#   3. the ID is the oriented variant edge, and it is the site key.
#
# This test reshapes the checked-in deconstruct catalog into exactly that form
# without needing pantree installed, runs collect-graph-variation on both the
# original and the reshaped catalog over the same GAF, and checks the reshaped
# run produces phased output. It validates the integration contract, not the
# biology -- the reshaped catalog contains no variants deconstruct missed.
#
# Usage: ./scripts/test_pantree_catalog_shape.sh [pgphase_binary]

set -euo pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
PGPHASE="${1:-$REPO/pgphase}"
DATA="$REPO/test_data/graph_chr20"
WORK="$(mktemp -d -t pantree_shape.XXXXXX)"
trap 'rm -rf "$WORK"' EXIT

for tool in bgzip tabix; do
    command -v "$tool" >/dev/null || { echo "error: $tool not on PATH"; exit 1; }
done
[[ -x "$PGPHASE" ]] || { echo "error: pgphase binary not found at $PGPHASE"; exit 1; }

echo "== reshaping deconstruct catalog into pantree form =="
python3 - "$DATA/chr20_25M.sites.vcf.gz" "$WORK/reshaped.sites.vcf" <<'PY'
import gzip, sys

src, dst = sys.argv[1], sys.argv[2]


def parse_walk(text):
    steps, i, n = [], 0, len(text)
    while i < n:
        if text[i] not in "<>":
            return None
        j = i + 1
        while j < n and text[j] not in "<>":
            j += 1
        steps.append((text[i + 1:j], text[i] == "<"))
        i = j
    return steps or None


def canonical_edge(a, b):
    render = lambda s: ("<" if s[1] else ">") + s[0]
    fwd = render(a) + render(b)
    rev = render((b[0], not b[1])) + render((a[0], not a[1]))
    return min(fwd, rev)


def edges(steps):
    return {canonical_edge(steps[i], steps[i + 1]) for i in range(len(steps) - 1)}


def classify(ref, alt):
    if len(ref) == len(alt):
        return "SNP" if len(ref) == 1 else "MNP"
    return "INS" if len(alt) > len(ref) else "DEL"


records, seen_ids = [], set()
with gzip.open(src, "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        info = dict(x.partition("=")[::2] for x in f[7].split(";") if x)
        at = info.get("AT")
        if not at:
            continue
        walks = at.split(",")
        alts = f[4].split(",")
        if len(walks) != len(alts) + 1:
            continue
        ref_steps = parse_walk(walks[0])
        if ref_steps is None:
            continue
        ref_edges = edges(ref_steps)
        for i, alt in enumerate(alts, start=1):
            alt_steps = parse_walk(walks[i])
            if alt_steps is None or alt == "*":
                continue
            novel = sorted(edges(alt_steps) - ref_edges)
            if not novel:
                continue
            # Mimic pantree: the site ID *is* the variant edge.
            site_id = novel[0]
            if site_id in seen_ids:
                continue
            seen_ids.add(site_id)
            out_info = [f"VT={classify(f[3], alt)}", f"AT={walks[0]},{walks[i]}"]
            for key in ("AC", "AN"):
                if key in info:
                    out_info.append(f"{key}={info[key].split(',')[i - 1] if key == 'AC' else info[key]}")
            records.append((f[0], int(f[1]), site_id, f[3], alt, ";".join(out_info)))

records.sort(key=lambda r: (r[0], r[1], r[2]))
with open(dst, "w") as out:
    out.write("##fileformat=VCFv4.2\n")
    out.write('##INFO=<ID=AT,Number=R,Type=String,Description="Allele Traversal as path in graph">\n')
    out.write('##INFO=<ID=VT,Number=1,Type=String,Description="Variant type">\n')
    out.write('##INFO=<ID=AC,Number=A,Type=Integer,Description="Alt allele count">\n')
    out.write('##INFO=<ID=AN,Number=1,Type=Integer,Description="Total alleles">\n')
    out.write(f"##contig=<ID={records[0][0]}>\n")
    out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for chrom, pos, site_id, ref, alt, info in records:
        out.write(f"{chrom}\t{pos}\t{site_id}\t{ref}\t{alt}\t60\tPASS\t{info}\n")
print(f"reshaped {len(records)} biallelic edge sites")
PY

bgzip -f "$WORK/reshaped.sites.vcf"
tabix -f -p vcf "$WORK/reshaped.sites.vcf.gz"

run_pgphase() {
    local sites="$1" tag="$2"
    "$PGPHASE" collect-graph-variation \
        --ref "$DATA/ref.fa.gz" \
        --sites "$sites" \
        --gaf "$DATA/HG002_chr20_25M.coord.gaf.gz" \
        -o "$WORK/$tag.tsv" \
        --phased-vcf-out "$WORK/$tag.vcf" \
        -t 4 > "$WORK/$tag.log" 2>&1 || {
            echo "FAIL: collect-graph-variation ($tag) exited non-zero"
            tail -20 "$WORK/$tag.log"
            exit 1
        }
}

echo "== baseline: deconstruct catalog =="
run_pgphase "$DATA/chr20_25M.sites.vcf.gz" baseline
echo "== candidate: pantree-shaped catalog =="
run_pgphase "$WORK/reshaped.sites.vcf.gz" reshaped

base_vars=$(grep -vc '^#' "$WORK/baseline.vcf" || true)
resh_vars=$(grep -vc '^#' "$WORK/reshaped.vcf" || true)
base_ps=$(grep -v '^#' "$WORK/baseline.vcf" | grep -c 'PS' || true)
resh_ps=$(grep -v '^#' "$WORK/reshaped.vcf" | grep -c 'PS' || true)

echo
printf '%-28s %12s %12s\n' "metric" "deconstruct" "pantree-form"
printf '%-28s %12s %12s\n' "phased VCF records" "$base_vars" "$resh_vars"
printf '%-28s %12s %12s\n' "records carrying PS" "$base_ps" "$resh_ps"
echo

status=0
if [[ "$resh_vars" -eq 0 ]]; then
    echo "FAIL: pantree-shaped catalog produced no variant records"
    status=1
fi
if [[ "$resh_ps" -eq 0 ]]; then
    echo "FAIL: pantree-shaped catalog produced no phased records"
    status=1
fi
# A flat biallelic catalog decomposes multiallelic snarls, so it should not lose
# a large share of the calls. A big shortfall means sites are being rejected.
if [[ "$base_vars" -gt 0 ]] && [[ $((resh_vars * 2)) -lt "$base_vars" ]]; then
    echo "FAIL: pantree-shaped catalog lost more than half the calls"
    echo "  check $WORK/reshaped.log for site rejection reasons"
    status=1
fi

if [[ "$status" -eq 0 ]]; then
    echo "PASS: collect-graph-variation ingests a flat, biallelic, edge-keyed catalog"
fi
exit "$status"
