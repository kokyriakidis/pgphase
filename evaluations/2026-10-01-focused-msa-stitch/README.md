# Focused MSA admission, complete transfer, and stitch continuity

## Retained fixes

1. **Use the observed binomial majority for conflicting boundary pairs.**
   The sparse-indel admission test used `2^-n`, the tail for unanimous votes,
   even when both parities were present. At 54,684,758–54,702,072 the seven
   callable pairs are 6:1. Their one-sided fair-parity tail is 0.0625, rather
   than 0.0078125. The corrected test admits the existing focused MSA retry.
   This is admission evidence, not a stitch certificate.
2. **Certify the final backfilled focused matrix.** The source-path check ran
   before exact-CIGAR backfill. At 4.767 Mb the raw focused matrix passed,
   whereas the matrix actually transferred had a weak cut at 4,778,792.
   Backfill now precedes the certificate; that retry is rejected and the
   original source keeps its correct connection. No weak cut is discarded.
3. **Keep private context rows of validated focused source blocks.** Source
   selection already kept whole BAM phase-set metadata, but candidate transfer
   dropped private rows outside the original seam. That could truncate a path
   and remove sites when a successful join eliminated a later discovery seam.
   Selected, cut-free focused blocks now transfer all their admitted phased
   private rows. Ordinary solves retain their existing seam scope. Source-path
   evidence is cached once for transfer and stitch metadata.
4. **Corroborate a single weak left graph edge with direct BAM SNP calls.**
   Connecting the 54.685 Mb gap makes its left graph block part of the block
   used at 54.894 Mb. An internal edge at 54,607,517–54,607,533 has GAF votes
   0/7/2 (haplotype 1 agreement / haplotype 2 agreement / reversal). Direct
   Q30 BAM pairs have 31 agreements and six reversals. A physical clean-SNP
   stitch can now validate this first failed edge using both haplotypes,
   one-sided fair-parity `p <= 0.001`, and signed base/mapping-quality odds.
   The graph prefix has already passed; the entire suffix must still pass.
   Dominant GAF reversal is a veto. Boundary evidence and both block paths
   remain required. Existing unanimous insertion-edge checks retain their
   earlier behavior; conflicting calls need the stronger count certificate.

The production changes contain no coordinates, competitor output, or parental
truth. Existing MSA/WFA is reused; no new realignment stage is added.

## Owning-chunk and chromosome checks

The 54–55 Mb replay joins the target deletion `CCT>C` and SNP `G>A` with the
same ALT polarity. Five private source rows outside that seam remain present
and phased, including separate complementary insertion/deletion rows at
54,585,084. The previously connected 54,894,127–54,912,022 pair remains joined
in that orientation. The owning chunk changes from 4,315 phased / 4,306
correct / 9 discordant reads to 4,323 / 4,312 / 11.

Full chr20, same BAM + graph catalog + GAF + reference, defaults, eight threads:

| Metric | Accepted baseline | Current |
|---|---:|---:|
| Output reads | 256,570 | 256,586 |
| Phased/truth-scored reads | 237,071 | 237,101 |
| Truth-correct reads | 229,892 | 229,920 |
| Discordant reads | 7,179 | 7,181 |
| Per-read-PS concordance | 96.971793% | 96.971333% |
| Read phase sets | 690 | 686 |
| VCF phase blocks | 340 | 340 |
| VCF records | 62,368 | 62,798 |
| VCF span N50 | 684,798 bp | 739,888 bp |
| Connected committed panel gaps in chr20 | 73/92 | 74/92 |

All 73 earlier chromosome panel spans survive. The only added panel connection
is 54,684,758–54,702,072; 14 tracked competitor targets and four intentional
controls remain open. N50 improves by 55,090 bp (8.04%). The four replaced
VCF keys and nearby new representations are listed in
`changed-representations.tsv`: three substitutions acquire competing inserted
alleles at their anchor, and an MNP changes its span. The existing output
conflict filter removes incompatible haplotype descriptions. No new writer
rule or row-merging rule was introduced. This evaluation is against the
accepted pgphase state; native competitor runs were not repeated.

Final command:

```bash
./pgphase collect-graph-variation \
  --ref test_data/chm13v2.0.chr20.renamed.fa \
  --bam test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam \
  --sites test_data/chr20.sites.striped.vcf.gz \
  --gaf test_data/HG002.chr20.annotated.coord.gaf.gz \
  -t 8 -r 'CHM13#0#chr20' \
  -o /tmp/pgphase-gap-next5/final-full/candidates.tsv \
  --phased-vcf-out /tmp/pgphase-gap-next5/final-full/phased.vcf \
  --phased-bam-out /tmp/pgphase-gap-next5/final-full/phased.bam
```

The accepted baseline is `/tmp/pgphase-gap-next4/anchor-full`. Final production
run time was 150.05 seconds. `metrics.json` and `chr20-panel-spans.tsv` retain
measurements without depending on these temporary output paths.

## Regression checks

- The existing 54.685 Mb panel case now replays its owning chunk and requires
  `spans=1`, concordance >=98%, separation >=83%, and the phased boundary
  deletion at canonical 54,684,759. The TOTAL span floor rises from 73 to 74;
  all earlier floors and required sites remain unchanged.
- A new owning-chunk case checks nine exact allele descriptions, including
  five private context rows, both target boundaries, and both downstream
  54.894 Mb boundaries. It checks allele polarity, common PS, read parental
  orientation, and absence of incompatible haplotype alleles. The old
  implementation fails the target join and separation checks.
- A new 4–5 Mb owning-chunk case checks that the 4.767 Mb insertion/SNP pair
  stays connected with its original allele relation and parental orientation.
  Certification before backfill fails the PS, allele-relation, and separation
  checks. Its new separation floor is the measured accepted owning-chunk
  baseline, 64%; the chromosome-wide output has a different denominator.

The first full-panel run used an incorrectly estimated 69% floor for this new
owning-chunk case. Both the accepted baseline and corrected production measure
64%. Only that new floor was corrected; the erroneous-certificate trial still
fails both phase and separation checks. The new test is rerun separately after
calibration. No pre-existing expectation is weakened.

## Rejected experiments

- Whole private-row transfer for all source blocks cost 421 correct chr20
  assignments and added 458 discordant ones; it is not retained.
- Broad MSA admission and converting all newly admitted reads to local-only
  allele calls were tested. Neither replaces the targeted admission rule.
- Correct admission alone closed the target but split existing 4.767 and
  54.894 Mb joins. Better aggregate read totals did not justify retaining
  those splits; final-matrix certification and the physical edge check fix
  them independently.
- Padding the focused region by another 50 kb increased local discordance;
  the original complete adjacent graph phase-set extents remain in use.
- An extra early physical stitch pass, retaining stale source metadata across
  recovery passes, freezing source-block polarities, and trying every earlier
  MSA indel anchor did not repair these failures and are not retained.

The five context keys initially described as MSA consensus losses were present
in the BAM source. Their disappearance was in candidate transfer. The source
consensus is not rebuilt when a confidently placed extra read is appended;
this corrects the earlier proposed explanation of that experiment.

Final validation: GCC build has no new warnings; shared unit binaries pass;
294 predicate, 27 port-parity, and 168,696 upstream-parity assertions pass.
The full native pipeline panel completes 51 cases / 4,010 assertions, with
only the incorrectly estimated new separation floor failing. The calibrated
new case passes in fresh native runs against both the accepted baseline and
final production (17 assertions each); the new 54 Mb case also passes in a
fresh native run (49 assertions). All 4,010 current assertions in 51 cases are
then rechecked successfully against the unchanged completed native outputs,
with the production binary SHA-256 checked before each cached replay.
`native-panel.log`, `panel.log`, and `validation.json` distinguish these runs.
`git diff --check` passes.
