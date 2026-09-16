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
#   * --joint-het-orientation -- iteration 1: without it a site whose own allele
#                               depths call it heterozygous can still be emitted
#                               homozygous, because each haplotype's consensus is
#                               recomputed independently by majority. chr20:48,204,383
#                               (DP 71, 30 ref / 41 alt, AF 0.577) is emitted 1|1
#                               without it and 0|1 with it, matching hiphase.
#   * -q 1                   -- the adopted MAPQ floor.
#
#   * --link-by-alleles     -- REQUIRED, and it is not the same option as
#                               gap_link_by_alleles. Options::link_by_alleles
#                               defaults to FALSE (phasing_types.hpp:430) while
#                               gap_link_by_alleles defaults to true (:310); the
#                               block-link vote reads the former. With it off the
#                               vote uses check_agree_haps, which needs the read
#                               to already carry a haplotype AND its allele at
#                               the left site to match that haplotype's own
#                               consensus, so at a noisy site nearly every read
#                               returns -1: on chr20 the 48,202,056 -> 48,204,383
#                               link had 57 reads with usable alleles at both
#                               sites and scored agree=1 conflict=0. With it on
#                               the vote uses check_agree_alleles, which scores
#                               the read's allele pattern directly.
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
    -r "${REGION}" -t "${THREADS}" -q 1 --keep-noisy-kmeans --link-by-alleles \
    --joint-het-orientation \
    -o "${OUT}/candidates.tsv" \
    --phased-vcf-out "${OUT}/native.vcf" \
    -b "${OUT}/phased.bam" \
    > "${OUT}/stdout.log" 2> "${OUT}/stderr.log"
