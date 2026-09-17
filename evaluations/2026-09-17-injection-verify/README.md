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

## Result over the six panel windows, 1,576 candidates

| window | candidates | dropped | duplicated | missing | alleleSet | depth | verdict |
|---|---:|---:|---:|---:|---:|---:|---:|
| 48,183,976 | 154 | 0 | 0 | 0 | 0 | 1 | 0 |
| 55,843,827 | 246 | 0 | 1 | 0 | 2 | 2 | 0 |
| 24,105,188 | 261 | 0 | 0 | 0 | 0 | 0 | 0 |
| 5,309,406 | 308 | 0 | 0 | 0 | 1 | 3 | 0 |
| 12,717,796 | 362 | 0 | 0 | 0 | 0 | 1 | 1 |
| 39,838,293 | 245 | 0 | 0 | 0 | 0 | 0 | 0 |
| **total** | **1,576** | **0** | **1** | **0** | **3** | **7** | **1** |

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
