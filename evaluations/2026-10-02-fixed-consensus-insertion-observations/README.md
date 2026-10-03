# Recover missing homopolymer insertion observations without changing BAM gauges

## Defect and correction

The ordinary haplotype-aware MSA uses one selected BAM phase set. Reads outside
its selected clusters can physically cover a candidate but contribute no
observation there. At 24.121713–24.131707 Mb, the source has separate four-T
and eight-T insertion rows. A spanning read with a seven-T CIGAR event can be
recognized by the existing one-edit local caller against these fixed paths,
but never reaches that caller in the initial recovery solve. This is an
observation-admission problem, not justification to copy HiPhase's binary
alleles onto both complementary rows.

The initial targeted BAM solve now requests local unplaced-read recall:

1. A cheap ordered scan of the fixed reference/consensus alignments checks for
   two lengths of the same inserted homopolymer at one reference coordinate.
   Other contrasts do not trigger extra read alignment work in this pass.
2. Reuse the existing WFA composition and coverage handling to evaluate both
   paths, without adding a read to either whole-region MSA cluster.
3. Restrict candidate admission to the requested seam. Separate insertion rows
   remain separate. Their one-base error neighborhoods must be disjoint; both
   composed paths and physical flanks must support a local call. A read must
   call the two rows complementarily. REF to both is a third allele, not proof
   of a haplotype. The existing allele-fraction check still applies.
4. Defer observations until ordinary discovery and phasing finish. An ordered
   exact-key/read-ID pass rejects inconsistent duplicates and fills only
   unknown observations at surviving MSA-verified phased heterozygotes. A
   previously phased source read keeps its own source block and must agree
   with that block's allele orientation. Rebuild counts and the read index.
5. Leave BAM genotypes and source gauges fixed. The usual graph transfer and
   stitching gates evaluate the expanded matrix. This mode neither authorizes
   a broad retry nor changes ordinary BAM operation.

No truth labels, competitor calls, fixture coordinates or new alignment
method enter production. No row merging or stitch-threshold relaxation occurs.

## Targeted result

Compare native owning 24–25 Mb replays, using input BAM read names that
physically overlap the nominated gap. Score HP after resolving each output
block's arbitrary orientation by parental majority; independent blocks are
not treated as a single connected block. Saved same-BAM HiPhase uses DV calls.

| Local metric | Starting pgphase | Updated pgphase | HiPhase DV |
|---|---:|---:|---:|
| Phased truth-scored overlaps | 60 | 89 | 93 |
| Correct assignments | 60 | 89 | 93 |
| Errors | 0 | 0 | 0 |
| Local read blocks | 3 | 3 | 1 |
| Dominant-block correct reads | 29 | 49 | 93 |

Thus 29 additional local reads are correctly phased, but the boundary pair
remains split. Equal local purity does not imply equal continuity. The owning
replay improves 3,895 scored / 3,893 correct / two errors to 3,924 scored /
3,922 correct / two errors. Its seven required rows and opposite insertion
alleles remain. The local overlap set contains 115 input read names.

The permanent 24 Mb test now requires the increased owning and local counts,
zero local errors, the separate insertion rows and the existing rejection of
a wrong flank join. It fails the starting binary in four coverage assertions.
The new synthetic regression checks fixed cluster membership and abstention
for literal-REF, different-position and tandem-repeat contrasts. The starting
alignment implementation fails its fixed-membership checks.

## Rejected trials

- Adding all recalled observations and rephasing the source changes later
  discovery. The whole chromosome loses 1,023 variant keys and increases
  errors from 7,127 to 9,469, despite promising local coverage. Rejected.
- An observation-only pass over insertion contrasts throughout the full
  recovery context keeps keys but increases chromosome errors by 32 and fails
  five existing assertions in four cases. Rejected; do not weaken their gates.
- Merely requiring exact local alleles still loses owning-chunk coverage and
  discards most of the correctly recoverable homopolymer reads. Rejected.
- Restricting coordinates alone is insufficient. Source membership, allele
  separation and the homopolymer contrast must also be respected.

The unrestricted fixed-source comparison and failed panel are preserved in
`rejected-context-wide-parity.json` and `rejected-context-wide-windows.txt`.
The all-site rephase comparison is in `rejected-all-sites-parity.json`.

## Final chromosome validation

Starting binary SHA256:
`74218715db7412bb7d0e88633b14ec6c85e3220047b6e9ce1d6df6210c7d7df4`.
Updated binary SHA256:
`d3b5ad20b227cc051d1d5280a6aa4fc82cbc1c1bc9393082aa5533005456a1c8`.
Use the same BAM, GAF, reference, catalog and default options. Outputs:
`test_data/tmp_gap_next27/partial-coverage/` and
`test_data/tmp_gap_next27/run-length-final/`.

| Metric | Before | After |
|---|---:|---:|
| Truth-scored phased reads | 237,170 | 237,199 |
| Correct assignments | 230,043 | 230,072 |
| Discordant assignments | 7,127 | 7,127 |
| Conditional read accuracy | 96.994983% | 96.995350% |
| Read phase sets | 661 | 661 |
| VCF variant keys | 63,630 | 63,630 |
| VCF phase blocks | 333 | 333 |
| VCF span N50 | 774,189 bp | 774,189 bp |

Exactly 29 read tags change, all previously unphased reads now truth-correct.
Every previously scored read retains its correctness state, and every
previously phased HP/PS tag remains. No VCF key, genotype or PS changes. Only
the two separate insertion records change counts, fractions and derived
quality: each gains 49 observations (DP 11 -> 60). Existing SNP gauges remain
coherent, and no old VCF block is merged into another.

The coordinate panel remains 86 spanned windows out of 100: no additional gap
closure is claimed. Four open coordinates are intentional controls; ten are
historical competitor targets. The newer 7.28/24.12 Mb boundary replays outside
that panel also remain split. The updated 24 Mb owning regression protects
coverage and rejects the wrong boundary orientation.

Build, all units, 981 predicate assertions in 42 cases, the 1,354 focused-window
assertions in three cases, and HiFi/ONT TSV/VCF goldens pass. HiFi one/four-thread
determinism holds. The only build warning is the existing unused-function
warning from abPOA's SIMD header. No panel expectation is lowered.
The complete window suite passes 7,038 assertions in 69 cases, preserving every
required connection, exact span expectation and parental-orientation guard.

The native chromosome run completes in 295.28 seconds at eight threads. The
105 fresh native panel requests complete in 395.26 seconds with four concurrent
workers. These are validation timings, not a controlled runtime comparison.
`native-cache.json` records the binary, argument vectors, input stats and hashes
of all native output files. `replay_cached_panel.py` verifies these identities
before reuse for scoring; uncached requests run natively. Starting-version or
rejected-version outputs are never scored as updated output.
