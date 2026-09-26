# Seam-local BAM phase-set transfer experiment (chr20)

The implemented transfer carries selected BAM source phase-set membership
through exact shared graph sites after the existing recovery stitch. It keeps
BAM-private candidate insertion inside the detected seam. A complete source
PS is evidence for the seam; sites outside the seam keep their graph PS.

## Controls and rejected arms

- **Positive:** BAM PS 51,262,081 has 155 oriented sites spanning
  51.262–51.414 Mb. Consecutive clean SNP pairs around the missed 51 Mb
  seam have 30, 4, and 49 informative observations. The final full-chunk
  VCF puts 51,262,081, 51,270,774, 51,286,463, and 51,287,372 in
  PS 51,262,082. The unsupported 51,235,063 site remains in PS 51,165,532.
- **Negative:** an unconditional whole-source-PS overlay joined graph blocks
  around 19.395 and 19.415 Mb through BAM PS 19,272,142, although one source
  cut had only one-haplotype support and another had 6 consistent versus 20
  conflicting links. The 19–20 Mb regional truth comparison changed from
  4,071 phased / 4,020 correct / 51 discordant to
  4,071 / 3,699 / 372 under that rejected overlay. The final transfer checks
  both haplotypes and conflict at every source cut; those two graph blocks
  stay separate.
- Materializing every private source-PS flank site across a 1 Mb solve
  reduced full-chromosome tagged reads by 61 and did not improve N50. The
  final transfer retains the original strict-seam candidate import.

## Matched full chr20

Both arms used the same graph catalog, GAF, annotated BAM, reference and
command settings (`--link-by-alleles --block-link-window 8
--chunk-size 1000000 -t8`). Truth is
`test_data/derived/chr20_truth_hap.tsv`; read accuracy counts the majority
parental assignment independently within each BAM PS. VCF N50 uses spans
of phase blocks with at least two records.

| Arm | Truth-evaluable phased reads | Truth-correct | Discordant | Read accuracy | Phased VCF hets | VCF blocks | N50 |
|---|---:|---:|---:|---:|---:|---:|---:|
| Baseline | 236,520 | 228,184 | 8,336 | 96.47556% | 61,645 | 457 | 456,233 bp |
| Seam-local transfer | 236,526 | 228,189 | 8,337 | 96.47523% | 61,646 | 464 | 456,233 bp |

The local connectivity bug is fixed without a large wrong join, but the
chromosome-wide gain is small: six additional phased reads, five correct and
one discordant, and no N50 gain. This result does not establish that every
connected BAM source PS can be safely translated through every graph seam.

Commands used for the A/B: `python3 /tmp/pgphase_n50_read_score.py` on the
phased BAMs and `python3 scripts/phase_block_stats.py` on the native VCFs.
Outputs are in `/tmp/pgphase-stitch-nextseam-default-chr20/` and
`/tmp/pgphase-wholeps-final2-chr20/` on this workspace.

## Validation

`make -j8 pgphase`, `make unit-tests`, `make window-tests`,
`make predicate-tests`, and `make check` passed. The window panel includes
the 51 Mb positive and 19 Mb negative-control assertions.
