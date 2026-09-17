# Is every site injected, and represented correctly, on stock defaults?

`verify_injection.py` runs six independent checks per panel window. Truth is
applied only to decide what a site **is**, never to select which sites to look
at, so a pass means the pipeline reached the right answer on its own evidence.

| check | what it flags |
|---|---|
| dropped | a candidate the alignment channel has and the hybrid does not |
| missing | a catalog claim read truth calls heterozygous, absent from our candidates |
| duplicated | two emitted records at one position |
| allele set | a locus with two clean parental modes and no reference read, described by one allele |
| depth | a record whose DP is under half the reads covering its window |
| verdict | a record classified homozygous that read truth calls heterozygous |

## Result over the six panel windows, 1,575 candidates

| window | candidates | dropped | duplicated | missing | alleleSet | depth | verdict | attributes |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 48,183,976 | 154 | 0 | 0 | 0 | 0 | 1 | 0 | 0 |
| 55,843,827 | 245 | 0 | 0 | 0 | 2 | 2 | 0 | 0 |
| 24,105,188 | 261 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 5,309,406 | 308 | 0 | 0 | 0 | 1 | 3 | 0 | 0 |
| 12,717,796 | 362 | 0 | 0 | 0 | 0 | 1 | 1 | 0 |
| 39,838,293 | 245 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| **total** | **1,575** | **0** | **0** | **0** | **3** | **7** | **1** | **0** |

**Injection is complete.** Nothing the alignment channel finds is dropped, and
every catalog claim read truth calls heterozygous is present. The `dropped`
column read 1 before this session's merge-preference fix: `55,896,396`, a 16 bp
homopolymer deletion whose merged two-allele record was displaced by the
graph-claimed single-allele version.

The `allele set` check was over-sensitive at first and read 51. Requiring **two
clean parental modes** before calling a single-allele record a failure brings it
to 3: of the 51 loci carrying no reference read, 48 have no allele pair to
express at all -- in a repeat tract the reads scatter across many net lengths and
neither haplotype has a modal one, so a single-allele record is not the wrong
description of a two-allele locus, it is the only description available. Flagging
those buried the three that are real.

## The remaining defects, each with its mechanism

## Fixed: the duplicate, and it was not where I said it was

The claim above that this needed the allele-assignment path was wrong, and
checking the one thing that settles it showed why: the hybrid's copy of the
alignment channel's record is **byte-identical** -- `DP=63`, `5/58`,
`AF=0.920635`, `CLEAN_HOM`, the same numbers the alignment channel reports. The
representation is copied directly and correctly. The defect was purely
**additive**: a second record beside it.

`augment_chunk_with_graph_sites` tries every catalog ALT for an *exact* key match
against the existing candidates, and when none matches it adds the first ALT as a
new graph-only candidate. At this locus the catalog carries a 2 bp insertion and
the alignment channel called a 1 bp one, so no key matched and the 2 bp claim was
added -- then counted as though the reads carried it.

The catalog's claim is that the locus varies. *Which* length is present there is
a read measurement, and the alignment channel has already made it. So before
adding, injection looks for an indel of the same type already called at that
position, and where it finds one the catalog's site is **dropped** -- the
alignment's record is left entirely alone, not even marked `graph_site`.

Marking it was the first version of this fix and it is the more dangerous one.
`graph_site` is what routes a record through `classify_graph_only_candidates`,
which re-derives the category from the allele fraction -- and that is precisely
the mechanism that turns a usable `NoisyCandHet` into `CleanHom`. At
`55,896,396` what made the record usable in the alignment channel was its
category, not its fraction, and replacing it with an AF-derived `CleanHom`
verdict dropped the locus from phasing altogether. A claim that cannot name the
allele the reads carry has nothing to add to a record that already measured it.

Both versions measure identically on the panel; the drop is the one without that
exposure. Verified at the locus: the hybrid's row is `DP=63`, `5/58`,
`AF=0.920635`, `CLEAN_HOM` -- the alignment channel's row, unmodified.

Result: one record at the locus, `T>TA` at `1|1:63:5,58`, matching the alignment
channel and matching truth. Panel candidates 1,576 -> 1,575, stock defaults
byte-identical, `--retry-unphased-with-bam` 4 of 6 spanned at 99.49% with 0
concordant-to-discordant.

### The original diagnosis, kept for the record

**`55,903,460` -- duplicated, and the deeper bug is the allele count.** The
alignment channel emits one record, `T>TA`, a 1 bp insertion, and read truth
agrees: 58 of 61 covering reads carry +1 at purity 0.603, so the locus is
homozygous for it. The catalog claims `T>TTA`, a 2 bp insertion, which injection
adds as a second record -- and that record accumulates **59 of 63 reads as alt
support** when truth puts exactly **one** read at +2. So the injected allele is
not merely redundant, it is counted as though the reads carried it.
`backfill_graph_candidate_counts` is faithful; it copies `prof.alleles[]`, so the
over-permissive attribution is in the per-read allele assignment for an injected
indel.

The deletion analogue of this is already solved: `make_colocated_deletions_exclusive`
attributes a read to one co-located record. **The insertion mirror is not safe**,
which is why it is not implemented here: for deletions "keep the longest record
the read supports" is correct because a longer deletion subsumes a shorter one,
but at this locus that rule keeps the +2 record and zeroes the correct +1. The
read's own insertion length is not in the profile, so resolving it means fixing
the assignment, not attributing after the fact.

**`55,890,331` -- a real het described as two single-allele records of different
type.** Truth is -28 paternal on 31 reads and +4 maternal on 20, both at purity
1.000, over 66 covering reads. We hold an `INS` record (DP 57) and a `DEL` record
(DP 36) at the same position, both `NOISY_CAND_HET`. The merges built this session
pair records of the same type; this locus needs a cross-type pair, where one
allele inserts and the other deletes relative to the reference.

**`5,339,364` -- the emitted alleles are not the locus' alleles.** Truth is +34
paternal and +38 maternal, no reference read; we emit a 5 bp deletion. No
representation fix reaches it, because neither allele was ever proposed.

**`12,735,895` -- verdict.** Counts were corrected this session (DP 8 to DP 61,
0 reference to 22), but the het/hom verdict still reads homozygous: re-deriving it
from the corrected counts is measured to cost 380 concordant-to-discordant reads
on the panel and was reverted.

## Is the alignment channel's record the right one to copy? Measured both ways

The record does survive injection intact -- that is now verified rather than
assumed. But "intact" and "right" are different claims, so the same six checks
were run against the **alignment channel's own candidates**, with no graph
channel involved at all:

| check | alignment channel alone | hybrid, after injection |
|---|---:|---:|
| candidates | 1,517 | 1,575 |
| dropped | -- | **0** |
| missing | 0 | **0** |
| duplicated | **8** | **0** |
| allele set | **3** | **3** |
| depth | **21** | **7** |
| verdict | **1** | **1** |

Two things follow, and they pull in opposite directions.

**Injection is not the problem any more; it is a net improvement.** The hybrid
carries *fewer* duplicated positions than the channel it copies from (0 against
8) and *fewer* starved records (7 against 21). The graph claim is what does it:
a claim promotes a locus and its counts are backfilled from the graph reads, so
records the alignment channel left at a fraction of their coverage reach it. On
the two checks where the hybrid is not better it is exactly equal.

**And the three remaining allele-set failures are the alignment channel's own**
-- the same loci, `55,890,331` (INS) and (DEL), and `5,339,364`. So is the one
wrong verdict, `12,735,895`. Copying the alignment channel more faithfully cannot
fix any of them, because at those loci its own record is the wrong description:
a single allele where read truth shows two non-reference modes, and a 5 bp
deletion where truth is +34 against +38.

That relocates the remaining work. It is not in the transfer and not in the
catalog -- it is upstream, in how the alignment channel constructs a site in a
repeat tract:

- `55,890,331` needs one locus with two alleles of opposite sign, where it
  currently produces an insertion record and a deletion record.
- `5,339,364` needs the alleles the reads actually carry to be proposed at all.
- `12,735,895` needs its het/hom verdict re-derived from counts this session
  already corrected, which as a blanket rule costs 380 concordant-to-discordant
  reads and so waits on a bridge gate.
