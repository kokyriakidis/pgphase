#!/usr/bin/env bash
# Usage: iterate.sh NAME BINARY [pgphase args...]   (env vars pass through)
# Runs the tier-1 gap panel and full chr20, prints one summary block.
set -euo pipefail
readonly name="$1" bin="$2"; shift 2
readonly root=/home/kokyriakidis/Downloads/pgphase
readonly inv="${root}/test_data/tmp_hybrid_investigation"
readonly py=/home/kokyriakidis/micromamba/envs/bench-phasers/bin/python
cd "${root}"
"${py}" -I "${inv}/gapbench.py" run "${inv}" "${bin}" "${inv}/gb/${name}" 1 -- --union-gap-phasing "$@" | tail -1
"${py}" -I "${inv}/gapbench.py" score "${inv}" "${inv}/gb/${name}" 1 | head -5
out="${inv}/full/${name}"
mkdir -p "${out}"
start=$(date +%s)
"${bin}" collect-graph-variation --ref test_data/chm13v2.0.chr20.renamed.fa \
  --bam test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam \
  --sites test_data/chr20.sites.striped.vcf.gz --gaf test_data/HG002.chr20.annotated.coord.gaf.gz \
  -r "CHM13#0#chr20" --threads 10 --union-gap-phasing "$@" -o "${out}/candidates.tsv" \
  --phased-vcf-out "${out}/phased.vcf" --phased-bam-out "${out}/phased.bam" > /dev/null 2> "${out}/stderr.log"
echo "full chr20 run: $(( $(date +%s) - start )) s"
"${py}" -I "${inv}/audit_hybrid.py" test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam \
  test_data/derived/chr20_truth_hap.tsv "${out}" u="${out}" 2>/dev/null > "${out}/audit.json"
python3 -c "
import json; d=json.load(open('${out}/audit.json'))['arms']['u']; r=d['reads']
print('FULL correct', r['correct'], 'discordant', r['discordant'], 'unphased', r['unphased'], 'PS', d['phase_sets'], 'N50', d['span_n50'], 'largest', d['largest_block'])"
truth_vcf=/home/kokyriakidis/Downloads/pgphase-eval-data/results/chr12-18-20-comparison/chr20/truth.vcf.gz
python3 -I "${inv}/switch_points.py" "${truth_vcf}" "${out}/phased.vcf" VARIANT --list | awk "/^SWITCH/ {p=\$4; if (p>=26000000 && p<32000000) c++; else o++} !/^SWITCH/ {line=\$0} END {print line, \"\tswitches_cen\", c+0, \"switches_arms\", o+0}"
pinned=/home/kokyriakidis/Downloads/pgphase-eval-data/results/pinned_chr20_scores/pinned_chr20.tsv
awk -F'\t' '$1=="hiphase" || $1=="longphase" || $1=="legacy_hybrid" {printf "PINNED %-14s correct %s discordant %s unphased %s PS %s N50 %s switches %s (arms %s, cen %s)\n", $1, $2, $3, $4, $5, $6, $8, $10, $11}' "${pinned}"
# Arms / centromere, and reads crossing a truth het ("phaseable") vs not.
truth_vcf=/home/kokyriakidis/Downloads/pgphase-eval-data/results/chr12-18-20-comparison/chr20/truth.vcf.gz
"${py}" -I "${inv}/region_reads.py" test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam \
  test_data/derived/chr20_truth_hap.tsv "${truth_vcf}" "u=${out}" > "${out}/regions.tsv"
awk -F'\t' '$1=="u" {printf "REGION %-15s correct %7s discordant %6s unphased %6s  phased %5.1f%%  correct %5.1f%%\n", $2, $3, $4, $5, 100*($3+$4)/$6, 100*$3/$6}' "${out}/regions.tsv"
awk -F'\t' '($2=="arms" || $2=="arms_phaseable" || $2=="arms_nohet") && ($1=="hiphase" || $1=="longphase" || $1=="legacy_hybrid") {printf "PINNED %-14s %-15s correct %7s discordant %6s unphased %6s  phased %5.1f%%\n", $1, $2, $3, $4, $5, 100*($3+$4)/$6}' \
  /home/kokyriakidis/Downloads/pgphase-eval-data/results/pinned_chr20_scores/pinned_chr20_regions.tsv
