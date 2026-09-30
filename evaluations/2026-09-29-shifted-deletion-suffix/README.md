# Shifted deletion and disconnected graph suffix, 2026-09-29

At chr20:13,830,800–13,844,727, HiPhase connects a left deletion with the
right SNP. The pgphase BAM source presents the deletion at 13,830,801, but
four of the right-SNP REF reads place the equivalent two-base deletion 17
reference bases later in a `(TG)` repeat. The physical deletion bridge's
16-base search missed those calls. Expanding sequence-verified equivalence to
32 bases yields six qualified deletion REF and three ALT calls and a decisive
quality-weighted relation to the right SNP.

A whole-block stitch remained unsafe: the left graph phase set has unsupported
SNP edges at 13,752,640–13,773,452 and 13,773,452–13,804,846. No primary BAM
read spans the second SNP pair. Clean and repeat indels lie inside that edge,
so an empty-coordinate-cut assumption was wrong. The suffix split now checks
candidate physical spans and chooses the candidate boundary with the fewest
previously tagged reads that call both sides. It leaves the prefix independent,
unphases crossing reads, and transfers only the suffix into the right block.
The suffix's own graph SNP path and the physical deletion bridge are both
required. In the owning 13–14 Mb replay, SNP 13,773,452 remains in the left
phase set while SNP 13,804,846, deletion 13,830,800, and SNP 13,844,727 share
the right phase set.

| Full chr20 | Before | After |
|---|---:|---:|
| Tracked HiPhase-correct gaps connected | 16/23 | 17/23 |
| VCF keys | 62,154 | 62,154 |
| Truth-scored phased reads | 236,847 | 236,847 |
| Correct / discordant | 229,088 / 7,759 | 229,088 / 7,759 |
| Accuracy | 96.7240% | 96.7240% |
| Truth-scored read phase sets | 706 | 708 |

All VCF genotype dosages remain unchanged; five phased GT strings invert
together in the transferred suffix. Every truth-scored read keeps the same
correctness outcome in the full run. The phase-set count rises because the
unsupported graph prefix is separated instead of being carried into the
right-hand join.
The owning-chunk replay has 2,775 correct and 179 discordant of 2,954 scored
reads both before and after the stitch. Outputs are
`/tmp/pgphase-second-pass-final-full/` and `/tmp/pgphase-gap1383-full/`;
the local trial is `/tmp/pgphase-gap1383-trace/`. A dedicated owning-chunk
regression checks the right join, independent prefix, and local truth.

The graph window's measured truth-separated fraction rises from 0.54 to 0.73;
its panel span changes from open to closed, raising the graph panel total from
45 to 46. The target's existing required-site check retains deletion
13,830,801; the new owning-chunk test checks both boundary sites and the
independent upstream graph SNP.

Validation: `make -j4`, `make unit-tests`, `make check`, and the full
`make window-tests` suite pass. The window suite reports 1,834 assertions
in 37 cases. No new compiler warnings were emitted.
