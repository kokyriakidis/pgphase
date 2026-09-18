# Stage 2 keeps stage 1 and refines it

The solve runs in two rounds -- `collect_var.cpp:2162` over `kCandGermlineClean`,
then `collect_phase_noisy.cpp:2054` over `kCandGermlineVarCate` (clean plus the
noisy classes) after each noisy-region pass. As shipped, the second round
**discarded the first**: `collect_phase.cpp` cleared every read's haplotype and
phase set, re-initialised every participating site's consensus, and swept outward
from a pivot chosen over the NEW site set, so a noisy site could overturn the
parity the clean sites had established.

longcallD resets the same way (`assign_hap.c:491`) while its own comment
(`collect_var.c:2939`) says a round should use the previously obtained phasing as
initialization, so anchoring is a knowing divergence from upstream, in the
direction upstream documented.

## What anchored means, concretely

Three changes, all inside `assign_hap_based_on_germline_het_vars_kmeans`:

| | reset (previous) | anchored |
|---|---|---|
| read haplotypes and phase sets | cleared | kept |
| consensus of a site the previous round decided | re-initialised to -1/-1 | kept, and PINNED through the iteration |
| Phase 1 outward sweep from a fresh pivot | runs | skipped; the incoming consensuses are the gauge |

Two details that matter. The gauge is the per-site consensus, not the read
labels: Phase 3 recomputes every read's haplotype from the consensus, so
preserving read labels alone would anchor nothing. And the pin must SKIP the
write rather than revert it afterwards -- the iteration's convergence check
compares against a snapshot taken at entry, so a revert-after would keep
reporting "changed" and never converge.

The pin list is taken before initialisation, over sites with a decided
consensus. Measured over `chr20:1-2,000,000`: 18 stage-1 calls and 14 anchored
stage-2 calls, every one pinning most of its sites -- 24 of 31, 22 of 28, 26 of
29, 54 of 73, 19 of 24, 5 of 8. So the anchoring is active, not a no-op.

## Measurement, graph + recovery arm, chr20 at -t 16

| | wall | bridged | tagged | blocks | discordant | read hamming |
|---|---:|---:|---:|---:|---:|---:|
| stage 2 resets | 131 s | 127 | 203,751 | 281 | 2,851 | 1.399% |
| stage 2 anchored | 118 s | 126 | 203,751 | 281 | **2,732** | **1.341%** |

119 fewer misplaced reads at identical yield -- same tagged count, same 281
blocks, same 56,032 records -- for one fewer bridged block. On
`chr20:1-10,000,000` the two are byte-identical, so the gain is in windows the
slice does not contain.

Now the default. `--no-anchored-stage2` restores the resetting behaviour.

Scope of the evidence: only the graph + recovery arm was run, per the standing
instruction to measure one arm. The alignment and hybrid arms share this code and
were NOT measured, so their baselines are unverified against this change.

Unit 4/4, window 66/66, predicate 130/130.
