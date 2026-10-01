# Remaining chr20 gap coverage and allele evidence (2026-09-30)

The current same-callset HiPhase-correct panel contains 15 noncentromeric
windows whose pgphase graph plus targeted BAM result is still split. I checked
primary alignments in the indexed shared HiFi BAM at the exact panel endpoints.
Every spanning alignment in these windows has MAPQ at least 30; none is rescued
by reducing the recovery MAPQ threshold. Production targeted recovery already
uses `min(opts.min_mapq, opts.recovery_min_mapq)`, with recovery MAPQ default 1.
The span counts below are physical coverage, not callable allele pairs.

| chr20 endpoint interval | Primary reads spanning both |
| --- | ---: |
| 3,963,987–3,981,464 | 7 |
| 10,727,690–10,746,628 | 3 |
| 15,351,845–15,367,755 | 8 |
| 19,373,922–19,395,544 | 1 |
| 21,130,205–21,148,338 | 4 |
| 21,159,070–21,172,487 | 18 |
| 21,594,343–21,612,458 | 1 |
| 21,823,066–21,844,359 | 2 |
| 23,792,418–23,806,565 | 4 |
| 35,328,965–35,347,817 | 6 |
| 35,498,368–35,516,845 | 3 |
| 36,332,599–36,354,890 | 2 |
| 50,548,245–50,562,066 | 12 |
| 57,854,341–57,866,713 | 11 |
| 58,366,458–58,385,702 | 2 |

Focused owning-chunk replays under `/tmp/pgphase-57854-current/` and
`/tmp/pgphase-15351-current/` show why those relatively covered seams remain
split. At 57.854 Mb, all 11 MAPQ-60 spanning reads call the right clean SNP,
but none has a callable deletion allele in the selected BAM source matrix. The
left locus has complementary one- and two-base deletion rows. The physical
CIGAR events do not agree on one right-SNP parity: one-base deletion carriers
split 2 REF versus 3 ALT at the SNP, while two-base deletion carriers split
1 REF versus 2 ALT. Other reads carry different deletion lengths or no
deletion. Forcing a whole-block link would discard this conflict.

At 15.351 Mb, eight MAPQ-60 reads span the left deletion and right insertion.
Two carry the left one-base deletion, but all eight call the right insertion
REF and the nearby right clean SNP REF. The right ALT haplotype is absent among
these spanning molecules, so they cannot independently determine diploid
parity. Earlier exact-ALT backfill trials restored some missing deletion
observations but closed no tracked gaps and worsened full-chromosome accuracy;
see `evaluations/2026-09-30-sparse-gap-link-audit/README.md`.

At 19.373 Mb, one read spans the left insertion and right clean SNP. Its left
four-base insertion is shifted in the CIGAR and has one base quality 22; the
existing Q30 physical caller correctly abstains. A Q10 diagnostic obtains one
ALT-to-ALT pair, but no independent intermediate-to-right ALT pair. The BAM
source places these ALTs on opposite haplotypes while HiPhase places them
together, so importing that whole BAM phase set would create the wrong join.
The shifted-insertion repair trials are documented in
`evaluations/2026-09-30-shifted-insertion-backfill/README.md`.

At 36.332 Mb, two reads span the left SNP and right four-base deletion. Both
call the left SNP ALT; one has a six-base deletion at the right locus and the
other no deletion. Neither certifies the right candidate allele. A broad
shifted-insertion backfill trial joined this seam but fell below its protected
parental-read accuracy floor, so it was reverted.

The gap-window diagnostic also had a misleading label: a joined boundary with
zero interior phased heterozygotes was reported as `NO SITES`. It now reports
`CONNECTED` first. Adjacent boundaries, which have no interior reference
position to sample, now report `ADJACENT SITES` instead of an uninitialized
coverage minimum. The 12,256,072–12,269,535 integration regression exercises
the connected case (31 assertions), and a focused two-assertion unit regression
exercises both adjacent outcomes. These test-only corrections change no phasing
output or gap expectation. No gap was newly closed by this audit.
