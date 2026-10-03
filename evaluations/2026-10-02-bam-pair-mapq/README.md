# BAM candidate-pair MAPQ and remaining boundary evidence

## Fix

`local_run_boundary_flip` previously counted working-matrix observations using
GAF MAPQ for every candidate pair. Two BAM-injected sites now use their original
BAM-channel observations and BAM MAPQ. This prevents both rejecting a good BAM
alignment because its GAF mapping is weak, and accepting weak BAM observations
because its GAF mapping is strong. Missing BAM calls remain missing. Graph/mixed
pairs retain their existing working-matrix eligibility. Shared read eligibility,
allele-class support, disjoint read-half agreement and binomial cutoffs remain.

The helper moves to the graph/BAM adapter so the actual implementation can be
unit tested. Eight regression assertions cover same/cross polarity, weak/unknown
BAM MAPQ, absent BAM calls, mixed pairs, and missing allele-class support. The old
implementation fails five assertions; the final implementation passes.

## Verified result

Final binary SHA256:
`be30679a86113af54f4a885b29361f8914e9bee6ff3870c698023efc9fe2b5e8`.

A fresh final-build eight-thread chr20 run takes **256.93 s**. All **256,601 read
HP/PS pairs** and all **63,492 VCF rows** match the preceding accepted run exactly.
There are **237,127** phased/scored reads, **229,974** correct, **7,153** discordant
(**96.983473%**), **673** read phase sets and **337** VCF blocks; span N50 remains
**756,878 bp**. This fixes a channel eligibility bug; it closes **no additional
fixture gap**. The preceding audit's 42 competitor-supported coordinate
nominations remain. Do not treat this as phasing parity with HiPhase.

Build and all standalone units pass without new warnings. Predicate tests pass
822 assertions in 37 cases. Five fresh owning-window regressions pass 142
assertions, protecting the long/equivalent insertion joins, complementary-row
focused recovery, the 64 Mb suffix and aggregate BAM stitching. The full
95-coordinate panel was not rerun and no expectations were relaxed.

## Six additional gap audits

`boundary-evidence.json` measures primary BAM coverage, raw BAM-channel pairs
at the actual recovered boundaries, and parental read concordance. Both tools
use the same input BAM/reference. Native HiPhase uses DeepVariant calls; the
saved-pgphase-callset HiPhase arm uses the Oct 1 callset, not the current VCF.
Truth and competitors are evaluation-only inputs.

| Gap | MAPQ30 physical spanning reads | Callable BAM pairs | pgphase correct/phased | Native HiPhase correct/phased |
|---|---:|---|---:|---:|
| 4,866,153–4,874,129 | 23 | 0 | 40/40 | 50/51 |
| 24,121,713–24,131,707 | 29 | 0 at either insertion row | 60/60 | 93/93 |
| 34,094,604–34,102,867 | 39 | 1 at either right row | 75/75 | 76/76 |
| 41,879,449–41,880,908 | 59 | mostly double-REF; complementary ALT classes unresolved | 62/81 | 78/81 |
| 61,738,239–61,747,506 | 30 | 2, one allele class | 79/83 | 91/92 |
| 64,128,828–64,134,226 | 27 | 18 exclusive deletion-ALT pairs, unanimous SNP relation | 44/45 | 48/49 |

The 64 Mb pair is promising local evidence: 12 SNP-REF/first-deletion-ALT and
six SNP-ALT/second-deletion-ALT molecules agree, while the seven double-REF
calls at the two deletion rows are excluded from that certificate. It does not
certify the complete adjacent blocks: their graph observation paths have
missing edges and the downstream block has contradictory edges elsewhere.
The existing separation regression remains intact; no whole-block join is
forced. The next investigation is a certified local component transfer that
preserves the unsupported prefixes and all independent read-only gauges.

## Reproduce

```bash
make -j8
make unit-tests predicate-tests
PGPHASE_TEST_WORKDIR=/tmp/pgphase-gap-next22/window-tests ./test_gap_windows \
  'observed BAM insertion runs*,equivalent BAM and graph insertions*,a single seam admits*,phased right reads bridge*,complete BAM blocks use*'
```

The two Python evaluation scripts expose `--help`. The recorded full-run
outputs and six owning-chunk matrices reside in `/tmp/pgphase-gap-next22/`;
`measure_output_parity.py` accepts arbitrary before/after output directories.
`audit_gap_pairs.py --runs-root /tmp/pgphase-gap-next22` uses `baseline-<chunk>`
matrices. It decodes each competitor BAM once to keep the audit fast.
