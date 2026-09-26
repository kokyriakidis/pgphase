# Chr20 graph and BAM site-representation audit

Date: 2026-09-23

## Inputs and method

The current standalone BAM caller and graph+BAM command were run on the same
annotated HG002 chr20 BAM and CHM13 reference. The graph command also used
`test_data/chr20.sites.striped.vcf.gz` and the annotated GAF, at its default
1 Mb chunk size. The BAM caller used its default 500 kb chunks. The graph run
emitted candidate, retained-site, filtered-site, and targeted-recovery audits.
The initial diagnostic graph VCF was byte-identical to the committed baseline
run (SHA-256 `44d510498d03dce5a1ba6c469d28dae615d0bc490894f281fd2bb7ff60939d21`).

Phased heterozygous VCF alleles were compared by the internal normalized tuple
`(sort_pos, type, ref_len, alt)`. This is an **exact representation** test, not
a haplotype-equivalence test: a multibase substitution can contain the same
SNP as a separate BAM record, and indels in repeats can shift or rotate. A
five-base proximity check with the same type, length, and ALT is only a triage
proxy, not proof of equivalence. Counts below are unique allele keys, so one
multiallelic record may contribute more than one key.

## Bugs found and corrected

1. **Split graph ALT misindexing.** The targeted recovery sequence index mapped
   every ALT from `site_meta` to each biallelic candidate. A split candidate
   carries the source VCF's full ALT list, but represents only the ALT selected
   by `site_allele_orig_idx`. The same selected-ALT rule now serves targeted
   recovery and whole-chunk BAM observation attachment. Whole multiallelic
   candidates do not accept an exact binary BAM match. In the complete chr20
   recovery audit, 24 candidates inside a recovery window changed from falsely
   “known graph” to correctly appended BAM candidates. The phased VCF gained
   15 rows and lost none. Another 124 shared VCF rows changed GT and PS as
   local block orientation changed, while their depths stayed the same. The
   corrected run phased 236,378 reads, with 227,799 truth-correct and 8,579
   discordant, versus 236,338 / 227,754 / 8,584 before this fix. It gained 55
   tags, lost 15, and changed the truth verdict of 29
   common tagged reads (19 improved, 10 worsened). Conditional truth accuracy
   moved from 96.3679% to 96.3706%.

2. **Truncated multibase substitution REF.** Graph candidates with an
   equal-length multibase replacement carry internal type `Snp` and
   `ref_len > 1`. The VCF writer used one REF base. For example, catalog
   `chr20:5,648,115 TG>CA` became `T>CA`, which describes a different
   haplotype. The writer now fetches the full `ref_len` span. Across chr20,
   exactly 952 VCF rows changed their REF and derived END; all other fields
   on those rows were identical.

3. **Lost deletion-replacement ALT.** A length-decreasing replacement can
   retain bases in `key.alt`. The writer emitted only the left anchor; for
   example, `chr20:29,573,098 TATAATAAAA>TT` became `TATAATAAAA>T`.
   The writer now appends `key.alt` and computes SVLEN as the net ALT minus REF
   length, matching the upstream longcallD definition. Forty chr20 rows
   changed ALT, 31 changed INFO/SVLEN, with 62 changed rows in their union.
   GT, PS, POS, and REF were unchanged by this writer fix.

The writer-only changes left the complete candidate TSV byte-identical. Every
one of the 62,153 final chr20 VCF records has a REF string matching the CHM13
FASTA at its reported position. All 333 records carrying SVLEN report the net
length of their emitted REF and ALT alleles.

## Remaining exact-key gap

The standalone BAM VCF contains 83,270 phased heterozygous allele keys. The
final graph+BAM VCF contains 61,742; 54,348 are shared exactly. Of the 28,922
BAM keys absent from the graph VCF, a candidate-first comparison gives:

| Exact-key status | Alleles | Interpretation |
|---|---:|---|
| No exact catalog ALT or final candidate | 19,301 | These need sample-private discovery or a representation-aware comparison. A different graph allele can still encode the same haplotype. |
| Exact catalog ALT, no exact final candidate | 4,008 | Graph observation, allele selection, or candidate filtering removed the exact allele. Filter-reason memberships include 1,892 `ref_only`, 1,259 `high_af`, 213 `no_reads_in_chunk`, 193 `low_af`, and 20 `low_depth`; 449 lack a filter row at a matching source ID. Memberships can overlap across source snarls. |
| Exact final candidate, no phased heterozygous VCF key | 5,613 | 5,534 are `REP_HET_INDEL` candidates excluded from clean output. The other 79 are clean candidates with a non-heterozygous genotype: 46 are 0/0 and 33 are 1|1. A positive internal PS label alone is not evidence of a heterozygous output call. |

Independently of that prioritization, **21,806** BAM keys have no exact catalog
ALT: 9,544 SNPs, 6,257 deletions, and 6,005 insertions. Of these, 2,505 did
reach an exact final candidate through BAM recovery; 19,301 have neither
source. A nearby same-type, same-length, same-ALT catalog key exists within five
bases for 317, but proximity does not prove equivalence. At least 153 other
apparently missing BAM SNP keys are substituted bases inside phased graph
multibase alleles. The candidate TSV does not expose `ref_len` for a one-base
complex insertion, so a small number of insertion classifications may remain
ambiguous.

The `high_af` group includes alternate-vs-alternate snarls whose per-pair
reference depth is low although two alternate alleles split the reads. The
previous controlled `--snarl-allele-phasing` experiment reduced such misses but
more than doubled discordant read assignments (248 to 518 on its recorded
chr20 comparison), so this audit does not enable that mode by default. The
internal MNP category is still `CLEAN_HET_SNP`; repeat-context screening of
these anchors remains a separate accuracy question documented in `CHECKPOINT.md`.

The previous complementary-deletion loss at 55,336,460 does not recur in the
current VCF: both the insertion and deletion rows are present.

## Multiallelic graph-site follow-up

The chr20 catalog contains 977,275 records, including 80,373 records with
multiple ALT alleles. Every record has one traversal per VCF allele (REF plus
ALTs), so this fixture shows no ALT/traversal-index mismatch at catalog load.
The compact walk matcher preserves the original allele index. After its
per-allele `min_alt_depth` filter, the default graph conversion creates one
biallelic REF/ALT candidate for each surviving ALT and remaps read observations
to those candidates. The source ALT index remains in `site_allele_orig_idx`;
the selected-ALT recovery fix above uses it.

In the current default chr20 run, 17,133 multiallelic catalog sites supply at
least one retained candidate, 2,416 supply at least two, and 63,240 supply no
retained candidate. These are **site counts**, not missing biological variants:
the last group includes sites without ALT support and intentionally filtered
alleles. The filtered-site dump contains 3,272 unique source IDs with at least
two distinct `high_af` ALT rows. This is direct evidence that the default
REF/ALT denominator can discard both sides of an ALT-versus-ALT site, although
it does not prove every one of those sites is a true diploid heterozygote.

Among the 4,008 BAM-phased exact allele keys present in the catalog but absent
from the final graph candidates, 1,689 occur in a multiallelic catalog record.
Of these, 794 have a `high_af` filter row at a matching source ID. The
remaining 2,319 occur only in biallelic catalog records. Exact allele matching
remains representation-sensitive, and a key may occur in multiple overlapping
snarls; these counts are triage, not a biological recall estimate.

`--snarl-allele-phasing` changes the comparison to ALT versus all other alleles;
`--snarl-keep-whole` can retain one n-allelic anchor. Neither is the default.
The earlier controlled chr20 experiment recorded 248 discordant reads at the
default, 518 with split ALT-versus-other phasing, and 393 with whole-snarl
phasing. Thus the default loses some real graph signal, but changing it for
every multiallelic site has measured accuracy cost. This follow-up found no
additional ALT-index or traversal-count bug and does not change that default.

## Verification

- `make -j16 pgphase test_graph_bam_adapter`
- `make unit-tests`, `make predicate-tests`, `make window-tests`, `make check`
- 232 window assertions and 153 predicate assertions passed
- Direct synthetic VCF regressions for split ALTs, multibase REF,
  length-decreasing replacement ALT, and net SVLEN
- Full-chromosome candidate-table identity before and after writer changes
- Full-chromosome REF verification against the indexed CHM13 FASTA and
  net-length verification of every emitted SVLEN

The exact-key categories are a prioritization tool. This audit does not certify
that the graph catalog and BAM callset contain the same biological variants,
or that every remaining unphased site should be emitted as a call.
