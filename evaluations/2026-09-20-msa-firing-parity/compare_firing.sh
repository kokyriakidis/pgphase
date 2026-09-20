#!/usr/bin/env bash
# Diff where the noisy MSA fires, pgphase against longcallD, over one span.
#
# Both tools log one line per noisy region with a three-way branch label chosen
# by the same rule: guided when a phase set carries both haplotypes, unguided
# when n_full_reads >= min_dp, otherwise skipped. Upstream prints Hap / NoHap /
# Skipped at -V 1 (align.c:1792, :1798, :1800); ours prints MsaFire at
# --verbose 1 (align.cpp, collect_noisy_reg_aln_strs).
#
# Usage: LCD=/path/to/longcallD/bin/longcallD ./compare_firing.sh REF BAM REGION
set -euo pipefail
REF=${1:?ref fasta}; BAM=${2:?bam}; REGION=${3:?region}
LCD=${LCD:?set LCD to the longcallD binary}
PG=${PG:-./pgphase}
W=$(mktemp -d)
"$LCD" call -t 1 -V 1 -o "$W/up.vcf" "$REF" "$BAM" "$REGION" 2> "$W/up.log"
"$PG" collect-bam-variation --ref "$REF" --bam "$BAM" -r "$REGION" -t 1 --verbose 1 \
      -o "$W/our.tsv" --phased-vcf-out "$W/our.vcf" > /dev/null 2> "$W/our.log"
awk '/^(Hap|NoHap|Skipped) /{split($2,a,":"); print a[2]"\t"$1}' "$W/up.log" | sort -u > "$W/up.fire"
awk '/^MsaFire /{print $3"\t"$2}' "$W/our.log" | sort -u > "$W/our.fire"
echo "upstream regions: $(wc -l < "$W/up.fire")   ours: $(wc -l < "$W/our.fire")"
join -t$'\t' -j1 "$W/up.fire" "$W/our.fire" > "$W/both"
echo "identical region bounds: $(wc -l < "$W/both")"
cut -f1 "$W/up.fire" | sort -u > "$W/up.reg"; cut -f1 "$W/our.fire" | sort -u > "$W/our.reg"
echo "upstream-only regions:   $(comm -23 "$W/up.reg" "$W/our.reg" | wc -l)"
comm -23 "$W/up.reg" "$W/our.reg" | sed 's/^/    only upstream: /'
echo "ours-only regions:       $(comm -13 "$W/up.reg" "$W/our.reg" | wc -l)"
comm -13 "$W/up.reg" "$W/our.reg" | sed 's/^/    only ours:     /'
echo "regions upstream logs more than once (processed in two chunks):"
cut -f1 "$W/up.fire" | sort | uniq -d | sed 's/^/    /'"
echo "branch disagreements:"
awk -F'\t' '$2!=$3 {printf "  %-24s upstream=%-8s ours=%s\n",$1,$2,$3}' "$W/both"
awk -F'\t' '$2!=$3' "$W/both" | wc -l | xargs printf "  total: %s\n"
rm -rf "$W"
