# Complete BAM block evidence and finalized stitching

## Design

Targeted recovery now snapshots every phased heterozygous site of each selected
BAM block, including sites outside the gap, with its callable allele observations.
The snapshot has independent site indices and source-scoped PS labels; graph
candidate ownership and allele representation are preserved. Duplicate molecule
names do not supply independent votes. The final stitch keeps both recovery
passes' matrices and considers every covering solve, even when a retry replaced
the active transfer gauges. A numeric label reused for a different source cannot
share a block-path certificate.

The new stitch scores complete source blocks and live graph blocks by molecule,
prioritizes clean SNPs, and translates the raw BAM allele basis to the current
live block basis. Both haplotypes must support a statistically decisive physical
MAPQ30/Q30 SNP relation. Each block must also independently agree with its current
read HP gauge; decisive conflicting observation channels veto the stitch.

All participating source and graph components must have supported internal paths.
Every edge is checked before a chain is mutated. A union changes only PS labels
and constant HP orientation, never genotypes or candidate representation. This
route runs after both recovery passes and source attachment, using the final live
block state. Production decisions use no truth or competitor data.

## Ordering regression

Applying complete evidence inside the earlier seam solver changed later genotype
recovery and source attachment. An initial chr20 trial lost 16 phased SNPs and
734 formerly correct read assignments. Path checks alone still lost 11 phased
SNPs. Moving the route after the shared solver preserved sites but still ran
before graph source attachment: the 37 Mb owner fell from 3,561 correct / 46
discordant reads to 3,415 correct / 191 discordant, despite a physical SNP edge
with 53 supporting molecules and no opposing votes.

The accepted ordering runs after source attachment and retries. On the 37 Mb
owner it preserves all 515 variant keys, every VCF row, and every read tag:
3,607 truth-scored reads, 3,561 correct, 46 discordant. The owning-chunk regression
allows a future correct closure but gates coverage, correctness, discordance and
parental orientation. Running it against the rejected binary fails all four of
those checks. The corrected binary passes its 11 assertions.

## Validation

The full chr20 run has exact read-tag and VCF-row parity with the accepted
baseline: 237,199 truth-scored reads, 230,072 correct, 7,127 discordant
(96.995350%), 63,630 variant keys and 332 VCF blocks. N50 remains 774,189 bp
and NG50 643,699 bp. The coordinate panel remains 87/101 connected; no additional
chr20 gap closes under these certificates. All previously accepted connections,
coverage and allele representations are preserved. The native eight-thread run
takes 291.01 s alongside panel validation; this is not a runtime comparison.

Final measurements are recorded in `validation.json`, `accepted-full-parity.json`
and `accepted-full-block-audit.json`. Window replays run the native candidate
binary. The replay manifest verifies the binary hash, normalized arguments, input
file metadata and every output hash; missing diagnostic requests run natively.
The manifest is only a test runtime optimization. No expectations are refreshed
or weakened.

Commands:

```bash
make -j4 pgphase
make unit-tests predicate-tests
make window-tests
./pgphase collect-graph-variation \
  --ref test_data/chm13v2.0.chr20.renamed.fa \
  --bam test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam \
  --sites test_data/chr20.sites.striped.vcf.gz \
  --gaf test_data/HG002.chr20.annotated.coord.gaf.gz \
  -t 8 -r 'CHM13#0#chr20' -o candidates.tsv \
  --phased-vcf-out phased.vcf --phased-bam-out phased.bam
```

All standalone unit binaries pass, with 1,193 predicate assertions in 44 cases.
The final gap suite passes 7,132 assertions in 71 cases, including the new
owning-chunk orientation regression. All 105 normalized native replays use the
final binary hash; fresh diagnostic requests also pass.

The physical certificate and internal-path gates deliberately leave unsupported
blocks independent. Complete evidence does not guarantee that every gap has an
identifiable, accurate join.
