# Contradictory allele edge and graph/BAM gauge in recovery stitch

## Defect and fix

In `stitch_recovery_phase_sets_left_to_right`, a graph-only edge computes a
shared-read graph/BAM gauge parity and `strongest_phase_set_edge` from paired
candidate alleles. The old branch marked the allele loci only when they
agreed, but its merge condition accepted **any** existing edge, even one with
opposite parity. The gauge could therefore flip an entire downstream block
against its directly observed molecule relation.

The merge now requires a present direct edge to agree with the gauge. On
conflict, the later ordinary allele-edge stitch can choose the direct parity.
No evidence threshold was changed. The synthetic `test_graph_bam_adapter`
case supplies 12 balanced, crossing clean SNP pairs and an opposing broad
gauge; it asserts the molecule-supported downstream orientation.

## Chr20 validation

Comparison against `/tmp/pgphase-8166-final-full/` using the same HiFi BAM,
GAF, graph catalog and parental truth (evaluation only):

| Metric | Accepted baseline | Contradiction veto |
| --- | ---: | ---: |
| VCF variant keys | 62,361 | 62,361 |
| Truth-scored tagged reads | 236,880 | 236,880 |
| Correct tagged reads | 229,143 | 229,590 |
| Discordant tagged reads | 7,737 | 7,290 |
| Read phase sets | 691 | 691 |

All 1,330 changed phased GT fields retain their PS labels. They are in two
peri-centromeric segments: 26,418,353–26,915,257 (1,078 records) and
29,270,077–29,303,608 (252 records). The large PS anchored at 25,971,072
falls from 547 to 142 discordant truth-scored reads; the PS anchored at
29,170,275 falls from 117 to 78. These per-PS changes explain 444 of the
447 fewer discordant reads; three other read assignments change. A temporary
vote trace confirmed that the direct edges behind those two vetoes had
9 same / 1 cross and 15 same / 0 cross paired observations, respectively;
they were not single-read nominations. The diagnostic logging was removed
from production. No target
noncentromeric gap was newly closed, and the centromeric region was not
used to choose or tune the fix.

Validation: `make -j8`, `make unit-tests`, `make window-tests` (3,788
assertions in 46 cases), `make check`, and `git diff --check` all pass.
Full-run artifacts: `/tmp/pgphase-gauge-veto-full/`.
