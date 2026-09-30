# New BAM seam retry on chr20, 2026-09-29

The graph seam at 15.015–15.071 Mb acquires two separate BAM phase blocks in
the first recovery pass. The second pass previously omitted the new 15,056,025–
15,071,132 boundary because its interval overlaps the old graph seam. The
nearest verified BAM SNP pair has four qualified physical calls (three for one
relation, one against), so the new phase-set pair is eligible for a bounded
retry. The owning 15–16 Mb replay joins boundary rows at 15,056,025 and
15,071,132 in PS 15,039,543 and scores 4,257/4,279 truth-correct reads,
versus 4,259/4,279 before the retry.

A trial that admitted every new overlapping pair joined many unrelated blocks
incorrectly: full chr20 truth accuracy fell from 96.72% to 86.57%, with 25,283
previously correct reads becoming discordant. A second trial requiring only two
spanning molecules admitted the adjacent 17,839,398–17,839,399 pair, although
no reference base lies between those anchors; its local BAM sub-solve tagged
16 reads with poor truth purity. The final gate requires callable, oriented
SNP alleles on both sides and skips adjacent anchors.

| Full chr20 | Before | Final |
|---|---:|---:|
| Tracked gaps closed | 15/23 | 16/23 |
| VCF keys | 62,154 | 62,154 |
| Truth-scored phased reads | 236,847 | 236,847 |
| Correct / discordant | 229,090 / 7,757 | 229,088 / 7,759 |
| Accuracy | 96.7249% | 96.7240% |
| Read PS / VCF PS | 709 / 350 | 706 / 348 |

The full baseline is `/tmp/pgphase-gap60033-final-full/`; the final trial is
`/tmp/pgphase-second-pass-final-full/`. The broad rejected trial is
`/tmp/pgphase-second-pass-span-full/`. The 15 Mb and adjacent 17 Mb owning
chunk trials are `/tmp/pgphase-gap15-snp3/` and
`/tmp/pgphase-gap17-noadj/`.

Validation: `make -j4`, `make unit-tests`, `make check`, and
`make window-tests` all pass; the window suite reports 1,820 assertions in
36 cases. The 15 Mb target's graph-arm panel expectation now records a closed
span and a truth-separated floor of 0.89 (measured 0.897). The adjacent
17.839 Mb regression requires local truth concordance of at least 0.99.
