# Two haplotype alleles at one position, and an undefined-behaviour read

Bug hunt with three instruments: a cross-record invariant scan over the emitted
VCF, an AddressSanitizer/UndefinedBehaviorSanitizer build, and read-level
scoring against the parental truth map.

## 1. One haplotype carrying two different alleles (106 positions)

Per-record auditing found nothing -- 62,485 records, no GT/AD/DP/AF violation.
The defect is only visible ACROSS records at a position: 508 positions carry
more than one record, 399 of them complementary (one allele per haplotype), and
**106 put two different ALTs on the SAME haplotype**, 105 of those inside one
phase set. 96 are a SNP together with the insertion that contains it --
`T>G` with `T>GAG` at 288,018, `T>A` with `T>AC` at 890,261 -- and 9 are an
insertion anchored on the reference base together with a SNP changing that base.
Applying both is incoherent, so one of each pair has to go.

**The reads decide, and they keep the longer.** Classifying every read at five
of these loci against both forms:

| position | reference | short form | long form |
|---|---:|---:|---:|
| 288,018 | 40 | **0** | 24 |
| 890,261 | 25 | **0** | 32 |
| 764,862 | 34 | 6 | 14 |
| 1,474,295 | 39 | 14 | 26 |

At the first two no read carries the bare SNP; at the other two the short form
appears on a minority of reads that share a parent with the long form, i.e. the
same event aligned two ways. Allele depth cannot arbitrate -- both records score
the same reads, 38,21 against 29,18 -- and the alignment arm is not a usable
oracle either: it emits the SNP alone at these loci (at 890,261 a `T>A` plus a
separate `T>C` at 890,262), splitting the event across records rather than
showing that a bare SNP is what the haplotype carries. Dropping the insertion
would report an allele no read has and lose the inserted bases.

So the rule is containment resolved towards the complete allele: the contained
record is dropped. For the 9 where neither allele contains the other, the two
records do not score the same reads and allele depth does separate them (at
1,907,006 the insertion carries 10 against the SNP's 5), so the better supported
one is kept, ties broken on the shorter allele for determinism.

Whole chr20: conflicting positions **106 -> 0**, records 62,485 -> 62,373, and
reads untouched -- 219,061 tagged, 323 read blocks, 2,543 misplaced, 1.161%,
333 phase sets. This is an output-representation fix and changes no phasing.

## 2. Undefined behaviour in the MSA refresh

An ASan+UBSan build over chr20:1-10,000,000 reported no memory error and two
UB sites, both the same one: `refresh_assigned_msa_observations`
(`collect_phase_noisy.cpp:1776`) binds `profiles[clu_read_ids[ci][ri]]` on an
EMPTY vector -- 18,579 times over that region.

`noisy_rvp` is allocated only on the single-consensus branch
(`collect_phase_noisy.cpp:1318`). On the two-consensus path the allocation lives
inside `update_cand_var_profile_from_cons_aln_str2`, which returns at its first
line when both haplotype variant lists are empty, leaving the vector empty while
the caller goes on to the merges and the refresh. The effect is nil -- with no
variants there is nothing to write, which is why it neither crashed nor changed
output -- but the reference binding is undefined and the guard is one line.

Fixed with a bounds check; chr20 output byte-identical with it in place.

## Gating it

The resolution is now `drop_conflicting_haplotype_alleles` in `collect_var.cpp`,
declared in `collect_var.hpp` and called from the graph writer, so it can be
tested directly. Six cases in `test_phase_predicates.cpp`; disabling the
function fails four of them.

Writing those tests corrected a misreading of the record shape. The first
version built the longer allele as a `VariantType::Insertion` at the same
`key.pos` as the SNP and failed: grouping is by `VariantKey::sort_pos()`, which
subtracts one for a non-SNP, so the two landed in different groups and nothing
was compared. Probing a real run at 288,018 shows what the pipeline actually
holds -- both records are `VariantType::Snp` at the same `key.pos`, differing in
`ref_len` (1 against 3) and alt (`G` against `GAG`). The tests are built that
way now, with one further case pinning the indel anchor convention: an
insertion recorded at pos+1 does group with the SNP at pos.

A panel-level assertion (`hap_allele_conflicts == 0`, asserted rather than
floored) was added to the window tests as well, but it is inert on this panel --
disabling the fix leaves it green, because none of the panel windows contains
one of the 106 positions, and a 100 kb region run around 288,018 does not
reproduce the pair either; the conflict needs whole-chromosome chunking. It is
kept as a cheap invariant, not counted as the gate.

Also recorded, since it cost a confusing half hour: building with sanitizers
leaves instrumented objects in `src/`, and a later plain `make` relinks them
WITHOUT the runtime, producing a binary that dies at startup with "ASan runtime
does not come first". Every run then fails and the window suite silently
collapses to 2 assertions. `make -B` after a sanitizer build.
