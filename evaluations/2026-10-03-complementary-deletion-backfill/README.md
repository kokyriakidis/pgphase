# Complementary deletion evidence for source retry

Target: chr20:50,548,245–50,562,066 (13,821 bp), outside the centromere.

## Bug

Both complementary BAM/MSA deletion rows at the left boundary are homopolymer
indels. Ordinary post-solve backfill excludes homopolymer rows because literal
CIGAR REF cannot identify their ALT-versus-other MSA contrasts. Consequently,
all 12 physically spanning MAPQ60 reads retained unknown calls at both rows.
The source conflict detector never saw the missing diploid evidence, so it did
not request the existing focused full-block BAM/MSA retry. Retaining more reads
or forcing a read-HP vote did not repair this representation problem.

## Accepted design

Retry admission examines two separate, co-located MSA deletion rows only when
their biallelic genotypes are complementary in the same source phase set.
A known-MAPQ30 read with Q30 flanks can temporarily supply a jointly missing
pair if its original CIGAR deletion is exact or sequence-equivalent to exactly
one ALT. That row receives ALT; the other receives ALT absence. Third lengths,
literal REF, additional MSA alternatives, compound events, unsupported qualities,
and existing observations do not qualify. Existing edit-equivalence code is
used; this adds no realignment stage or merged candidate.

The extra calls select an existing focused source retry. Both read profiles and
their interval index are restored before retry validation or evidence transfer.
Only the accepted retry's validated BAM/MSA projection can enter the graph
matrix and read rescue. Candidate keys, discovery counts, genotypes and phase
labels are untouched by the temporary diagnostic. Production contains no
parental truth, competitor calls or coordinate-specific decisions.

The focused source covers the complete adjacent graph phase blocks within the
owning chunk. Its accepted block is transferred with the existing source gauge,
path and atomic stitch checks. The two deletion rows remain separate. Their
ALTs lie on opposite haplotypes; the two-base deletion ALT and right A→G SNP ALT
lie on the same haplotype.

## Rejected broader trial

Persisting those physical calls in the imported matrix also closed this gap,
but promoted homopolymer calls to new read-rescue evidence elsewhere. Full chr20
truth-scored reads rose 237,199 → 237,507 while discordance rose 7,127 → 7,254
and accuracy fell 96.995350% → 96.945774%. The 334 new assignments included 125
errors. That trial is rejected; `full-parity.json`, `block-audit.json` and
`read-transition-audit.json` preserve its measured failure.

## Owning-chunk regression

The short 114 kb replay remains disconnected: it lacks complete adjacent block
context. The committed panel now replays the complete 50–51 Mb owning chunk.
It asserts connection, preserved separate deletion rows, read coverage and
parental orientation, with at least 3,887 truth-scored / 3,876 correct reads and
at most 11 discordant reads. The accepted output has exactly those counts.
The previous binary has 3,877 scored / 3,867 correct / 10 discordant reads and
fails five of the regression's 27 checks (`owning-before.log`). The accepted
binary passes all 27 (`owning-after.log`). No accuracy floor or discordance
ceiling was weakened.

Unit tests cover exact and shifted ALTs, third alleles, compound edits, missing
quality, low MAPQ, separate source gauges, out-of-window calls, idempotence and
preservation of counts and observations. A control verifies ordinary backfill
does not promote these physical homopolymer calls to rescue evidence. The
original implementation fails six assertions in the new pair test; the accepted
predicate suite passes 1,325 assertions in 45 cases.

The complete retry also changes one SNP projection at 50,562,679 beside a newly
retained insertion, and removes a low-depth SNP at 50,719,983 beside two retained
BAM insertion alleles. These are recorded explicitly in the block audits. The
CIGAR audit at the first site finds most paternal reads carrying the insertion,
while maternal reads carry literal C; at the second, both parents predominantly
carry a CIGAR T. These representation changes are not treated as evidence of a
correct literal-SNP genotype. The full-block parental read regression remains
the gate against a haplotype switch.

## Reproduction

```bash
make -j4 pgphase
make unit-tests predicate-tests
make window-tests
```

Native chr20 evaluation uses the same reference, graph catalog, GAF and BAM as
the accepted baseline, with eight threads and the default graph recovery:

```bash
./pgphase collect-graph-variation \
  --ref test_data/chm13v2.0.chr20.renamed.fa \
  --bam test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam \
  --sites test_data/chr20.sites.striped.vcf.gz \
  --gaf test_data/HG002.chr20.annotated.coord.gaf.gz \
  -t 8 -r 'CHM13#0#chr20' -o candidates.tsv \
  --phased-vcf-out phased.vcf --phased-bam-out phased.bam
```

## Full chr20 result

| Metric | Before | Accepted fix |
|---|---:|---:|
| Truth-scored phased reads | 237,199 | 237,209 |
| Truth-correct reads | 230,072 | 230,081 |
| Discordant reads | 7,127 | 7,128 |
| Read accuracy | 96.995350% | 96.995055% |
| Truth-scored read phase sets | 659 | 656 |
| VCF phase blocks | 332 | 331 |
| Block span N50 | 774,189 bp | 790,093 bp |
| NG50 | 643,699 bp | 643,699 bp |
| Tracked connected intervals | 87/101 | 88/101 |

The only newly connected tracked interval is the target. No other old VCF
blocks merge. The formerly right-hand block retains 629 of its 630 SNP gauges;
the single changed projection and removed low-depth SNP are disclosed above.
Ten new truth-scored assignments are correct. Five formerly correct reads
become discordant and four formerly discordant reads become correct, for a
net gain of nine correct reads and one discordant read. This is a small measured
tradeoff, not an assertion of unchanged or improved accuracy. The broad trial's
125 erroneous new assignments are absent.

The native eight-thread full run took 289.57 seconds, concurrent with panel
precomputation; this is not a competitor runtime benchmark. Output parity and
orientation audits are in `full-narrow-parity.json`,
`full-narrow-block-audit.json`, `owner50-narrow-parity.json` and `ng50.json`.
Final binary SHA256:
`094186bd72b9a2b96cc888519d7e326ff205b8ef63d74d519554a0245beed6cc`.

Final window validation passes **7,155 assertions in 72 cases**. All standalone
unit binaries and 1,325 predicate assertions in 45 cases pass. The cache contains
106 normalized requests: 105 fresh native panel runs plus the fresh owning-chunk
regression run. Every cached result verifies the final binary hash, input file
metadata, exact normalized arguments and output hashes. Additional diagnostic
requests missing from the cache execute the native pipeline. No existing floor,
ceiling or connected span was weakened. Only the graph arm's newly validated
target span and aggregate connected count changed. `validation.json` and
`window-after.log` record the final gates.
