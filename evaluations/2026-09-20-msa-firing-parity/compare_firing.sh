#!/usr/bin/env bash
# Diff where the noisy MSA fires, pgphase against longcallD, over one span.
#
# Both tools log one line per noisy region with a three-way branch label chosen
# by the same rule: guided when a phase set carries both haplotypes, unguided
# when n_full_reads >= min_dp, otherwise skipped. Upstream prints Hap / NoHap /
# Skipped at -V 1 (align.c:1792, :1798, :1800); ours prints MsaFire at
# --verbose 1 (align.cpp, collect_noisy_reg_aln_strs).
#
# Parsing note: upstream ALSO emits "Skipped region: <chrom>:<beg>-<end> ..."
# from the noisy-region ratio filter in collect_var.c, which is a different
# message and is not an MSA firing. Both patterns below therefore require a
# <chrom>:<beg>-<end> field immediately after the label.
#
# Usage: LCD=/path/to/longcallD/bin/longcallD ./compare_firing.sh REF BAM REGION
set -euo pipefail
REF=${1:?ref fasta}; BAM=${2:?bam}; REGION=${3:?region}
LCD=${LCD:?set LCD to the longcallD binary}
PG=${PG:-./pgphase}
W=$(mktemp -d)
trap 'rm -rf "$W"' EXIT
TAB=$(printf '\t')

"$LCD" call -t 1 -V 1 -o "$W/up.vcf" "$REF" "$BAM" "$REGION" 2> "$W/up.log"
"$PG" collect-bam-variation --ref "$REF" --bam "$BAM" -r "$REGION" -t 1 --verbose 1 \
      -o "$W/our.tsv" --phased-vcf-out "$W/our.vcf" > /dev/null 2> "$W/our.log"

awk '$1 ~ /^(Hap|NoHap|Skipped)$/ && $2 ~ /^[^ ]+:[0-9]+-[0-9]+$/ {
        split($2, a, ":"); print a[2] "\t" $1 }' "$W/up.log" | sort -u > "$W/up.fire"
awk '$1 == "MsaFire" && $3 ~ /^[0-9]+-[0-9]+$/ { print $3 "\t" $2 }' "$W/our.log" \
    | sort -u > "$W/our.fire"

cut -f1 "$W/up.fire"  | sort -u > "$W/up.reg"
cut -f1 "$W/our.fire" | sort -u > "$W/our.reg"

echo "MSA firings on $REGION"
echo "  upstream regions: $(wc -l < "$W/up.reg")   ours: $(wc -l < "$W/our.reg")"
join -t"$TAB" -j1 "$W/up.fire" "$W/our.fire" > "$W/both"
echo "  identical region bounds: $(wc -l < "$W/both")"
echo "  upstream-only regions:   $(comm -23 "$W/up.reg" "$W/our.reg" | wc -l)"
comm -23 "$W/up.reg" "$W/our.reg" | sed 's/^/      only upstream: /'
echo "  ours-only regions:       $(comm -13 "$W/up.reg" "$W/our.reg" | wc -l)"
comm -13 "$W/up.reg" "$W/our.reg" | sed 's/^/      only ours:     /'
echo "  branch counts upstream: $(cut -f2 "$W/up.fire" | sort | uniq -c | tr -d '\n')"
echo "  branch counts ours    : $(cut -f2 "$W/our.fire" | sort | uniq -c | tr -d '\n')"
echo "  branch disagreements: $(awk -F'\t' '$2 != $3' "$W/both" | wc -l)"
awk -F'\t' '$2 != $3 { printf "      %-24s upstream=%-8s ours=%s\n", $1, $2, $3 }' "$W/both"
