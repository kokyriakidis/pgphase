# Verified BAM SNP boundary and symmetric graph path checks

## New gap

Chr20:41,880,908–41,885,033 is outside the excluded centromeric region.
Both boundary keys have different PS labels in the accepted full-chromosome
baseline (`test_data/tmp_gap_fix50/full-final`). Saved native DeepVariant
HiPhase on the same BAM spans the interval with 87/91 locally phased reads
correct. These scores are evaluation only; production reads no truth or
competitor output.

## Defects and fix

- Recovery exposes a left MSA insertion opposite a verified noisy BAM SNP,
  but that new boundary was absent from the original stitch-seam snapshot.
  Revisit this supported boundary type within the existing recovery windows.
- The right A→T SNP at 41,885,033 is MSA and alignment verified, but its
  noise category excluded it from the clean-SNP bridge. Admit it only after
  its imported observations independently confirm the nearest following
  clean SNP's gauge. Here 20 MAPQ-30 paired calls, ten from each allele class,
  agree. Both deterministic read halves agree and the binomial tail passes
  the existing 0.001 bound. This does not change any verified allele call.
- The right block's weak graph SNP at 41,962,349 has both haplotypes on its
  incoming edge but loses one on its outgoing edge. The old bypass only
  skipped the right node of a failed edge. Add the symmetric left-node
  bypass with the same two-haplotype and significance requirements. The
  direct flanking pair at 41,958,183–41,965,652 has 17/7 supporting reads
  and zero opposing reads.
- A previously certified BAM join inside the right block has one opposing
  GAF vote. For this verified noisy-SNP boundary, reusing the certificate
  requires no statistically significant GAF reversal at the existing 0.01
  bound. Every original graph block still needs its internal path checked.
  Other bridge types retain the zero-reversal certificate rule.

The indel-to-SNP physical likelihood, allele-class requirements, source path
checks, exact inserted sequences, and transferred GT/observations remain in
place. No new realignment, coordinate-specific production rule, or truth
input is introduced.

## Owning 41–42 Mb replay

| Measure | Before | After |
| --- | ---: | ---: |
| Truth-scored phased reads | 4,154 | 4,156 |
| Correct | 4,126 | 4,129 |
| Discordant | 28 | 27 |
| Concordance | 99.32595% | 99.35034% |
| VCF blocks | 5 | 4 |
| VCF keys | 965 | 965 |
| Block span N50 | 341,860 bp | 460,682 bp |

All 4,126 previously correct reads remain correct; two new reads are correct
and one discordant read is corrected. Exact complementary +14/+18 CA-repeat
rows retain AD 22,31 and 32,21. Their ALTs stay opposite; the shorter insertion
ALT and both right SNP ALTs belong to the same parent.

For gap-overlapping reads, pgphase improves 87/89 to 90/91 correct;
saved native DV HiPhase has 87/91. Pgphase still has two read PS labels
among these overlaps despite one connected VCF block. The read grouping and
VCF block count are distinct measurements.

`transfer-audit.json` accounts for all 27,385 independent source calls and
27,385 checked calls, with no conservation errors. All 37,791 old overlay
calls, including 16 conflict abstentions, survive unchanged.

## Rejected trials

Applying certificate reuse globally joined 7.26 Mb but affected a protected
37 Mb replay. A broader exposed-insertion retry also joined the 36.62 Mb
control and lost seven previously phased reads. Merely narrowing that seam
route while leaving global certificate reuse enabled split three previously
closed chromosome intervals. None of these trial policies is retained. The
final certificate exception and added seam retry are restricted to the
verified noisy-SNP boundary; the general symmetric bypass retains its
original two-haplotype statistical certificate.

## Reproduction and tests

Use `test_data/tmp_gap_fix51/pgphase-before` for the saved baseline and current
`./pgphase` for the final implementation:

```bash
python3 test_data/tmp_gap_fix51/probe.py baseline 41
python3 test_data/tmp_gap_fix51/probe.py scoped41 41
make unit-tests
make predicate-tests
make check
```

The new committed owning-chunk test passes 34 assertions and fails eight
assertions on the baseline binary. The panel now has 105 coordinate cases;
this case asserts `spans=1`, read accuracy, parental orientation, and exact
complementary allele/depth retention. Its two insertion rows are also in the
required-sites list. The older 41.900 Mb regression now requires the intended
outer connection, retaining its downstream stitch and truth checks.

The broad trial's complete native suite ran 8,126 assertions in 78 cases and
identified three failures. Final focused reruns pass all 105 assertions in
the new case, the updated downstream 41 Mb case, and the 37 Mb counterexample;
the unchanged 36 Mb control passes 20 panel assertions. No accuracy floor or
control span expectation was weakened.

## Final default chromosome run

Final binary SHA256:
`69d998cab0b2319f6875e3514a7be4eaa1dc9b514419831614aeb012a13c4386`.
The checked-in-source rebuild and the binary used for the chromosome run are
identical. Outputs are in `test_data/tmp_gap_fix51/full-scoped/`.

| Measure | Accepted baseline | Final |
| --- | ---: | ---: |
| Output reads | 256,610 | 256,610 |
| Truth-scored phased | 237,328 | 237,330 |
| Correct | 230,232 | 230,235 |
| Discordant | 7,096 | 7,095 |
| Concordance | 97.010045% | 97.010492% |
| Scored read PS labels | 654 | 652 |
| VCF blocks | 331 | 330 |
| VCF keys | 64,188 | 64,188 |
| Block span N50 | 806,449 bp | 806,449 bp |
| Connected coordinate cases on expanded panel | 91/105 | 92/105 |

Every previously correct read stays correct. Two previously unphased reads
are correctly tagged and one previously discordant read is corrected. No
variant key or phased SNP is lost. Each original SNP block keeps one coherent
orientation; no old connected panel window opens. The sole new coordinate
closure is 41,880,908–41,885,033. The broad trial's 7.26 Mb join is not retained.

Build, all standalone units, 1,536 predicate assertions in 47 cases, HiFi/ONT
BAM golden outputs and HiFi thread determinism pass. Focused checks pass 105
assertions plus the 20-assertion control rerun.

Independent final verification also passes **8,106 assertions in all 78 native
cases**, using outputs from **110 fresh pipeline requests** to the exact final
binary. The prebuilt test executable still expected the reviewed 41 Mb outer
connection to remain split; that is the only failure in its fresh final-binary
run. The updated test source differs only in that assertion and its adjacent
comment. CLI construction and case ordering are unchanged. Rebuild the test,
then rescore the same native outputs, verifying their hashes and restoring the
hashed auxiliary matrices after tests remove old debug files. Every normalized
CLI request matches its saved request. This is documented in
`independent-native-cache.json`, `independent-native-requests.jsonl` and
`independent-replay-current.py`; the passing log is
`independent-window-tests.log`. No floor, ceiling or other assertion changes.

A separate full-chr20 run of the same binary confirms exactly one new gap,
91 to 92 connected coordinates among the 105 tracked cases, all 64,188 variant
keys and genotype allele sets retained, and every previously correct read
still correct. Two newly tagged reads are correct and one old error becomes
correct. See `independent-full-parity.json`, `independent-read-transitions.json`
and `independent-block-audit.json`. The isolated build introduces no warnings;
standalone units and HiFi/ONT golden checks also pass independently.
