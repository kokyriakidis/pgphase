#!/usr/bin/env bash
# THE arm. One configuration, iterated on -- not a sweep.
#
# This is the alignment machinery run over the union of the alignment's own
# candidates and the graph's catalog sites, in one solve per chunk, with nothing
# added on top:
#
#   * no --recover-gaps      -- the second pass is an addition, and it zeroes
#                               every non-graph candidate's category before
#                               phasing and skips the noisy-region MSA outright.
#   * no --retry-unphased-with-bam -- same reason, it exists only because the
#                               above withholds the sites.
#   * no --min-read-margin   -- min_read_hap_margin defaults to 0 and the hybrid
#                               subcommand does not override it; passing 2 was
#                               our own addition and on
#                               chr20:36,217,274-36,268,291 it stripped 271 reads
#                               of which truth says 254 (93.7%) were correct.
#   * --keep-noisy-kmeans    -- skip_noisy_kmeans = true is a hybrid-specific
#                               override; the alignment pipeline does not apply
#                               it, and the evidence inside gaps is exactly that
#                               class.
#   * -q 1                   -- the adopted MAPQ floor.
#
# gap_link_by_alleles is already true by default, so --link-by-alleles is
# omitted rather than passed: it changes nothing and reads as if it did.
set -euo pipefail

readonly REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
readonly OUT="${OUT:?set OUT}"
readonly REGION="${REGION:?set REGION}"
readonly THREADS="${THREADS:-1}"

mkdir -p "${OUT}"
"${REPO}/pgphase" collect-hybrid-variation \
    --ref "${REPO}/test_data/chm13v2.0.chr20.renamed.fa" \
    --bam "${REPO}/test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam" \
    --graph-sites "${REPO}/test_data/chr20.sites.striped.vcf.gz" \
    --gaf "${REPO}/test_data/HG002.chr20.annotated.coord.gaf.gz" \
    -r "${REGION}" -t "${THREADS}" -q 1 --keep-noisy-kmeans \
    -o "${OUT}/candidates.tsv" \
    --phased-vcf-out "${OUT}/native.vcf" \
    -b "${OUT}/phased.bam" \
    > "${OUT}/stdout.log" 2> "${OUT}/stderr.log"
