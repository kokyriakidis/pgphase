#!/usr/bin/env bash
# Derive the read -> parental-haplotype map the window tests score against.
#
# The truth BAM is aligned to the two parental assemblies, so a read's haplotype
# is the contig it aligns to (chr20_MATERNAL / chr20_PATERNAL) with the HO:Z:
# tag as a fallback for a merge that dropped the suffix. The rule below is the
# same one scripts/evaluate_phase_accuracy.py applies, kept identical on purpose
# so a window test and a whole-chromosome evaluation cannot disagree about truth:
#
#   -F 0x904   drop secondary, supplementary and unmapped records
#   MAPQ       >= MIN_MAPQ (0 by default, matching how this project's truth has
#              been built: a MAPQ floor of 1 drops 804 records whose placement
#              on the parental assemblies is ambiguous, and excluding them here
#              while the whole-chromosome evaluation keeps them would make the
#              two disagree about 0.3% of reads)
#   hq:i:      >= MIN_HAPQ when the tag is present (60 assumed when absent)
#   first qualifying record per read name wins
#
# Output is a two-column TSV, read name then MATERNAL or PATERNAL. It is derived
# rather than committed because the truth BAM is 1.5 GB and the map is ~9 MB;
# regenerate it after changing the truth BAM.
set -euo pipefail
readonly REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
readonly TRUTH_BAM="${TRUTH_BAM:-${REPO}/test_data/diplinator_eval_chr20/diplinator_merged.tagged.bam}"
readonly OUT="${OUT:-${REPO}/test_data/derived/chr20_truth_hap.tsv}"
readonly MIN_MAPQ="${MIN_MAPQ:-0}"
readonly MIN_HAPQ="${MIN_HAPQ:-0}"

[[ -f "${TRUTH_BAM}" ]] || { echo "truth BAM not found: ${TRUTH_BAM}" >&2; exit 1; }
mkdir -p "$(dirname "${OUT}")"

samtools view -F 0x904 -q "${MIN_MAPQ}" "${TRUTH_BAM}" \
| awk -v min_hapq="${MIN_HAPQ}" 'BEGIN { FS = OFS = "\t" }
    {
        hap = ""
        if ($3 ~ /_MATERNAL$/)      hap = "MATERNAL"
        else if ($3 ~ /_PATERNAL$/) hap = "PATERNAL"
        hapq = 60
        for (i = 12; i <= NF; ++i) {
            if ($i ~ /^HO:Z:/ && hap == "") { sub(/^HO:Z:/, "", $i); hap = $i }
            else if ($i ~ /^hq:i:/)         { sub(/^hq:i:/, "", $i); hapq = $i + 0 }
        }
        if (hap != "MATERNAL" && hap != "PATERNAL") next
        if (hapq < min_hapq) next
        if ($1 in seen) next          # first qualifying record wins
        seen[$1] = 1
        print $1, hap
    }' > "${OUT}.tmp"

mv "${OUT}.tmp" "${OUT}"
echo "wrote ${OUT}: $(wc -l < "${OUT}") reads"
