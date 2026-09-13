# Independent automatic recovery validation

Implementation `5316f9c`; binary fingerprint recorded alongside this report.
Selected before running: the two largest HiPhase-bridged gaps on chr12 and chr18,
and the next two largest chr20 gaps after the previously tested 15 Mb and 47 Mb cases,
from the frozen 2026-09-12 comparison. Exact source rows and solve regions are in
`manifest.json`. Each solve covers the audited gap plus 50 kb on either side.

Run `python3 evaluations/2026-09-13-auto-gap-validation/run.py`, then
`python3 evaluations/2026-09-13-auto-gap-validation/summarize.py` from this checkout.
Requires local test inputs, pysam, samtools, and read truth under
`~/Downloads/pgphase-eval-data`. Heavy outputs go to `/tmp/pgphase-auto-gap-validation`.
Commands are retained per case. No whitelist or threshold sweep was used.
Both arms use allele linking, block-link window 8, read margin 2, chunk size 500 kb,
and two threads. Recovery uses the default MSA margin 24.

| Audited gap | Clean → recovery phase blocks | Added phased reads | Internal gaps joined |
|---|---:|---:|---:|
| chr12:46725702–46838918 | 4 → 4 | 0 | 0/3 |
| chr12:63758735–63866834 | 2 → 1 | 40 | 1/1 |
| chr18:46185393–46302851 | 3 → 3 | 0 | 0/2 |
| chr18:33319321–33429051 | 4 → 3 | 0 | 1/3 |
| chr20:61737962–61810469 | 2 → 1 | 36 | 1/1 |
| chr20:17607780–17679082 | 2 → 2 | 0 | 0/1 |

Three of eleven residual hybrid gaps joined: chr20 61.8 Mb at the SNP tier,
and the other two joins at the indel tier. Two solve
regions reached one phase block; one improved partially. The other three did
not improve. The initial hybrid gaps differ from the frozen graph-only gaps
because the initial hybrid pass already incorporates clean BAM evidence.

All 2,293 initially tagged reads survived, with a single uniform PS/parity
transformation per original block. SNP evidence remained cumulative; indels
appeared only at tier 3; recovery stopped after successful joins. There were
76 additional phased reads. Read-truth discordance and switchflip counts did
not increase in any case. The chr18 46.2 Mb case retained one existing discordant
read and one switchflip; all other cases had zero in both arms.

Unresolved proposals often match both flanks but in separate proposal phase
sets: `LEFT_LINK=1, RIGHT_LINK=1, STATUS=partial` does not mean a bridge exists.
MSA evidence alone therefore does not guarantee connectivity. No preservation
or tier-ordering bug was detected by these checks.

These are targeted local runs, not a chromosome-wide NGC50 benchmark. Competitor
bridging is established by the frozen audit, not a new competitor run. Accuracy
here is read-truth accuracy, not a fresh variant-truth switch-error assessment.
Build and all unit tests passed. The existing `make check` runner still uses
unsupported CLI syntax, as documented in the implementation report.
