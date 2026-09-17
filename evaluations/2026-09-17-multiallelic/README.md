# Representing a locus whose haplotypes both differ from the reference

Two loci in the panel are the same shape: neither haplotype is the reference, so
no single biallelic record can describe them.

| locus | maternal | paternal | reference reads |
|---|---|---|---|
| `55,896,395` | 3 bp deletion, 32 reads | 16 bp deletion, 32 reads | **none** |
| `39,846,791` | ~21 bp insertion | a different ~21 bp insertion | **none** |

Emitted as two independent biallelic records, each one scores the other
haplotype's reads against its own allele. That is where the impossible records
come from: `0 reference / 30 alt` at allele fraction 1.000.

## Three fixes, all in the default path

**1. The multiallelic merge ran only under gap recovery.**
`merge_msa_insertion_alleles` -- which merges two co-located insertion records
into one record with two ALTs and a per-read vote -- was gated behind
`opts.recover_gaps && opts.private_msa_admit_all_in_region`, so the default
pipeline never reached it. A locus with two non-reference haplotypes needs that
treatment in every arm, so both merges now run unconditionally.

**2. There was no deletion equivalent.** `merge_msa_colocated_deletions` merges
two deletions of different length at one position the same way: the longer
deletion's span becomes REF, and each ALT is the anchor plus the bases that
allele retains, so `msa_insertion_alts` carries `{"", retained}` and the existing
multi-ALT machinery applies unchanged. The VCF writer's deletion branch ignored
`msa_insertion_alts` entirely and now honours it.

**3. The multi-allele site caller treated an empty query as reference.** In
`msa_site_event_allele`, the multi-allele branch returned allele 0 when the read
carried no query sequence at the site. For an insertion that is correct -- no
inserted bases means no insertion. For a deletion it is exactly backwards: an
empty query means the **whole footprint is deleted**, which is an allele, not the
reference. That scored 34 reads carrying the 16 bp deletion as reference and left
the allele at zero support, and the record genotyped `1|0`. The reference test is
now event-aware: `query == ref` for a deletion, `query.empty()` for an insertion.

## Measured

The alignment channel now describes both loci correctly:

| locus | before | after |
|---|---|---|
| `55,896,395` | two records: `CTTT>C` at `1\|0`, 30/28 and `CTTT...T>C` at `0\|1`, **0/30** | **one record** `REF=CTTTTTTTTTTTTTTTT ALT=CTTTTTTTTTTTTT,C`, **`1\|2`**, AD **0,35,34**, AF 0.507/0.493 |
| `39,846,791` | two records, each `AF 1.000` with **no reference reads**, genotypes `1\|0` and `0\|1` | **one record**, **`1\|2`**, AD **0,29,23** |

Truth at `55,896,395` is 32 maternal 3 bp and 32 paternal 16 bp deletions with no
reference read, so `0, 35, 34` is the right shape, and `1|2` is the right
genotype. Hiphase emits one record at each locus too.

Panel: **0 concordant -> discordant, 0 concordant tags lost, 0 records lost, 0
records gained**, 2,794 tagged at 99.68% -- unchanged, because the hybrid arm does
not yet reach these records (below).

## What the hybrid still does not do, and why gate-by-gate is the wrong route

In the hybrid the merged record is built correctly -- it is in the candidate table
as `REFlen=16 ALT=TTTTTTTTTTTTT,` DP 69, 0 reference / 69 alt, `NOISY_CAND_HET` --
and then two things stop it.

First it loses the same-key merge to the graph-claimed candidate, whose counts
`backfill_graph_candidate_counts` filled in as 2/35 and which
`classify_graph_only_candidates` then read as `CleanHom` at AF 0.946.

Second, and the real blocker: the multiallelic machinery is reachable only
through gap recovery. A chain of `gap_link_supported` conditions --
`collect_phase.cpp:234` (read-to-consensus scoring returns 0), `:678` (het
seeding requires `recovery_graph && gap_link_supported`), `:711`, `:749`, `:904`,
`:944`, `:955`, `:1196`, `:1242`, and `collect_pipeline.cpp:808`/`:837` -- exclude
a record carrying two ALTs unless a gap link vouched for it.

Relaxing those one at a time was tried and measured: letting the merged call win
the same-key merge and relaxing the scoring, seeding and inheritance gates
produced **zero new records and lost one** -- the locus that used to emit a wrong
`0|1` emitted nothing at all, because the record still earned no haplotype
(`PS = 0`, `HAP_ALT = 0`). Reverted. Trading a wrong record for no record is not
progress, and six more gates of the same kind remain.

The coherent change is to make "a record carrying two ALTs with reads behind
each" a first-class heterozygote across the whole assignment path, rather than a
gap-recovery special case -- one change with one measurement, not eleven. The
alignment channel already produces that record correctly, so the hybrid's job is
to carry that representation through, not to re-derive a verdict from a
recomputed allele fraction.
