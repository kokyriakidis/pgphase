# Second-pass phasing of chr20:15,019,294–15,130,077

The existing hybrid MSA → k-means → graph-block merger can close this gap.
A local diagnostic recovers 211 additional truth-concordant reads and joins
both graph blocks, preserving every original graph assignment up to one
uniform orientation change per block. This is a regional experiment, not a
validated chromosome-wide policy. The MSA score default remains 24.

## Reproduce

From the repository root:

```bash
make -j"$(nproc)"
make unit-tests
OUT=/tmp/pgphase-gap-rephasing scripts/test_chr20_gap_rephasing.sh
```

Requires the existing `test_data/` chr20 BAM/reference/catalog/GAF and the
frozen chr20 evaluation data under `~/Downloads/pgphase-eval-data` (override
`DATA_ROOT`). `PGPHASE` can select another built executable. The script runs
only pgphase; it reads frozen graph and HiPhase outputs. It writes exact hybrid
commands, binary/whitelist hashes, candidate tables, allele matrices, BAMs,
read-truth evaluations, endpoint votes, and an assertion that graph blocks
retain their internal orientation.

The solve window is `CHM13#0#chr20:14,950,000–15,180,000`, providing flanking
anchors around the 110,783 bp target gap. `private_sites.vcf` contains the
previously audited 26 BAM-derived candidate keys from
`/tmp/chr20_competitor_chain_private.vcf`; it scopes noisy MSA intervals. It
is a fixed regional fixture, not an automatic genome-wide gap selector. The
four clean private SNPs are retained independently of the whitelist. Truth
is read only after phasing and stitching have finished.

## What sites are present?

The catalog query overlaps 1,450 records (including records starting before
its query boundary). The frozen graph run emits 22 candidates in the exact
gap: 18 repeat het indels, three clean het SNPs, and one clean het indel.
Catalog presence does not imply usable heterozygous observations.

The fresh native BAM run finds six clean het SNPs and 21 MSA het candidates
inside the gap. Four clean SNPs have no catalog record at their exact position:

| position | depth | alternate fraction |
|---:|---:|---:|
| 15,039,543 | 58 | 0.603 |
| 15,039,634 | 58 | 0.603 |
| 15,048,452 | 73 | 0.452 |
| 15,055,707 | 75 | 0.453 |

The scoped hybrid second pass admits ten MSA het candidates inside the gap,
including the boundary deletion at 15,095,642. The left boundary deletion at
15,019,256 is just outside the coordinate-defined gap. MSA and shared-call
indels use different representations; exact coordinate matching alone cannot
establish whether two indels are equivalent. All candidate rows and their
categories are in `gap_candidates.tsv`.

## Measured outcome

Counts cover the 230 kb solve window, using assembly read truth and a minimum
of five reads per evaluated phase set. These are read metrics, not NGC50 or
shared-VCF switch metrics.

| arm | evaluated reads | discordant | phase sets | read switch/flip events |
|---|---:|---:|---:|---:|
| Frozen graph | 495 | 0 | 2 | 0 |
| Frozen HiPhase | 970 | 8 | 1 | 8 |
| Fresh native BAM caller/phaser | 788 | 6 | 6 | 6 |
| Hybrid clean pass | 555 | 0 | 3 | 0 |
| Scoped MSA, margin 24 | 653 | 0 | 3 | 0 |
| Scoped MSA, diagnostic margin 1 | 653 | 0 | **1** | 0 |
| Graph + margin-24 proposal, stitched | 528 | 0 | 2 | 0 |
| Graph + margin-1 proposal, stitched | **706** | **0** | **1** | **0** |

Both MSA arms use SNP-first escalation, then eligible indels, with
`--link-by-alleles --block-link-window 8 --min-read-margin 2`. The same
candidate-category counts produce different connectivity because the actual
read observations differ:

* At 15,019,256 → 15,039,543, margin 24 leaves **zero** shared observed reads;
  margin 1 recovers **four**, with allele pairs 00=2 and 11=2.
* At 15,095,642 → 15,109,301, margin 24 leaves **zero**; margin 1 recovers
  **15**, with 10 opposite and five same-allele pairs. This second link is
  noisy; the diagnostic result is not evidence that each individual link is
  sufficiently reliable for general deployment.

These counts exclude missing/negative allele observations; physical read-span
coverage alone would overstate evidence. See `boundary_observations.tsv`.

The accepted graph join has 196 winning shared reads at its weaker endpoint,
passes the existing 10-read/5-margin/90%-purity/both-haplotypes gates, and keeps
all 495 original graph reads. Graph PS 14,729,749 is flipped uniformly into PS
15,115,387; the latter block keeps its orientation. Internal graph phasing is
preserved; literal HP integers and PS IDs necessarily change when joining
oppositely labelled blocks. `audit.json` verifies this invariant.

The merger's distance parameter measures **PS start-ID separation**, not gap
width. These IDs are 385,638 bp apart although the gap is 110,783 bp, so the
local margin-1 diagnostic explicitly uses a 400,000 bp cap. The margin-24 arm
uses the existing 300,000 bp cap. These are separate gates; successful MSA
phasing alone does not bypass the merger. No global distance default changed.

## Bugs fixed in the current working tree

1. The within-chunk orientation loop swapped diploid alleles twice, producing
   no change. It now swaps once.
2. Equal agree/conflict votes could join phase sets. Ties now abstain and the
   window search can continue to another preceding anchor.
3. `--private-sites` unconditionally set graph authority in CLI parsing,
   overriding additive-mode handling downstream. Authority is now decided
   by the existing pipeline logic, or explicitly requested by the user.
4. Region admission cleared the whitelist pointer, which also disabled
   replacement of an existing repeat candidate by its validated MSA call.
   Indel escalation now replaces such calls and transfers their observations.
5. Sorting MSA candidates did not remap original read-profile indices. The
   merge now preserves the site ownership of alleles and query positions.
6. The new read-only snapshot restore mixed old read HP/PS labels with new
   candidate orientations and could undo joins. It was removed; graph
   preservation uses the existing complete-block orientation/stitching step.
7. Noisy-region escalation treated PS 0 as invalid instead of negative PS.
8. Included dependency files could make an object the default `make` target,
   leaving the executable stale. `.DEFAULT_GOAL := all` now ensures rebuilds.

Regression tests first failed for the double swap, tied link, region-mode
repeat replacement, and unsorted profile mapping, then passed after fixes.
`make -j$(nproc)` and all five `make unit-tests` binaries pass without new
warnings. `make check` cannot reach its golden comparisons: its pre-existing
runner passes the unsupported `--phased-vcf-output` option and positional
FASTA/BAM arguments. The unrelated `.devcontainer` whitespace change was
also left untouched.

## Limits and interpretation

Changing HP1↔HP2 uniformly within a phase set is valid relabelling, not an
error. The previous claim that a percentage of changed integer tags alone
proved degradation was insufficient; read truth and block-relative orientation
are needed. The independent graph-preservation assertion here checks that.

This result does not require a replacement phasing core. It establishes a
working regional route through existing processes and identifies the missing
boundary observations. It does not establish that margin 1 is safe across
chr20, chr12, or chr18. The defaults remain conservative, and broader-panel
read/variant truth evaluation is still needed before promoting the diagnostic
admission policy. A block join is not, by itself, proof of a chromosome NGC50
gain. Production automation of regional site selection is also outside this
fixed-fixture test.
