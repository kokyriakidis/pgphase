# Independent gap check: chr20:47,671,540–47,762,233

Checked after pushing gap-repair commit `663b84a`. This was the next-largest
HiPhase-bridged chr20 gap in the frozen correct-bridge table, selected before
running the experiment. No additional C++ changes or threshold search were
made. The same clean/margin-24/margin-1 comparison was used.

**The 90,693 bp target gap closes with the existing repaired second pass and
stitcher.** Final output adds 70 evaluated reads to the frozen graph baseline
without losing or internally reorienting any original graph block. The final
220 kb window still has two blocks because a third, more distant graph phase
set is rejected by the unchanged 300 kb distance gate.

## Reproduction

```bash
bash evaluations/2026-09-13-gap-rephasing-47m/run.sh
```

Defaults: `OUT=/tmp/pgphase-gap-rephasing-47m`,
`DATA_ROOT=~/Downloads/pgphase-eval-data`, `PGPHASE=<repo>/pgphase`.
The script uses the existing chr20 test inputs and frozen graph/HiPhase BAMs.
It never reruns competitors. The solve window is
`CHM13#0#chr20:47,600,000–47,820,000`.

A fresh native BAM call selects 23 biallelic heterozygous keys within the gap
plus 5 kb padding, using GQ >=10 and allele fraction 0.30–0.70. Input GT
phasing and PS are cleared. These keys scope MSA; additive mode retains clean
BAM candidates independently. The generated whitelist is saved here as an
audit artifact; rerunning regenerates it from the native caller. Truth is
only loaded after all phasing and stitching runs.

## Sites and missing evidence

The catalog query overlaps 1,103 records. The frozen graph candidate table
contains six clean heterozygous SNPs and 19 repeat het indels in the exact gap.
All six hybrid clean SNP positions already have graph catalog records:
47,671,540; 47,689,418; 47,713,869; 47,726,592; 47,738,673; 47,762,233.
Unlike the first region, absent clean SNP positions are not the main problem.

The fresh native BAM run emits five clean heterozygous SNPs and 20 MSA het
candidates inside the gap. Hybrid obtains the sixth clean SNP from the graph.
The scoped second pass admits five MSA het candidates in the gap. Both MSA
margins produce the same candidate-category counts but different observations:

| boundary pair | shared observed reads, margin 24 | margin 1 | allele-pair counts at margin 1 |
|---|---:|---:|---|
| 47,694,120 → 47,713,869 | 0 | 5 | 00=2, 01=2, 10=1, 11=0 |
| 47,746,345 → 47,755,891 | 0 | 20 | 00=0, 01=5, 10=13, 11=2 |
| 47,751,481 → 47,755,891 | 0 | 38 | 00=12, 01=4, 10=12, 11=10 |

Negative/missing allele observations are excluded. The first link is only
3:2 in its preferred orientation; this is a weak local edge despite the
truth-correct aggregate result. This experiment does not establish that
margin 1 is a safe general policy. See `boundary_observations.tsv` and
`gap_candidates.tsv` for the measured evidence.

## Results

Read-truth results across the 220 kb solve window (minimum five reads per
evaluated PS), not chromosome-wide NGC50 or shared-VCF switch metrics:

| arm | evaluated reads | discordant | phase sets | read switch/flip events |
|---|---:|---:|---:|---:|
| Frozen graph | 462 | 0 | 3 | 0 |
| Frozen HiPhase | 825 | 5 | 1 | 5 |
| Fresh native BAM caller/phaser | 703 | 3 | 3 | 3 |
| Hybrid clean pass | 510 | 0 | 3 | 0 |
| Scoped MSA, default margin 24 | 510 | 0 | 3 | 0 |
| Scoped MSA, diagnostic margin 1 | 510 | 0 | **1** | 0 |
| Graph + margin-24 proposal | 532 | 0 | 3 | 0 |
| Graph + margin-1 proposal | **532** | **0** | **2** | **0** |

Both MSA arms use SNP-first escalation and
`--link-by-alleles --block-link-window 8 --min-read-margin 2`. Both stitch
runs use the same 300,000 bp PS-start-distance cap, plus the existing
10-shared-read, 5-vote-margin, 90%-purity, both-haplotypes gates.

The join from graph PS 47,666,163 to PS 47,762,233 has 35 winning endpoint
reads and distance 96,070 bp; it is accepted. The additional PS 47,121,694
has 175 winning reads but distance 640,539 bp; it is rejected. Therefore the
selected target closes while the broader window retains two graph blocks.
All 462 original graph assignments survive with consistent block orientation.
`summarize.py` asserts preservation and the target endpoint join.

The retained margin-24 and margin-1 stitched read counts are identical: in
this region the diagnostic improves connectivity, not the number of output
reads. Relative to the frozen graph, either proposal recovers 70 reads.
No new implementation bug was required to explain this second case; it
independently reproduces the loss of usable boundary observations at margin
24. Defaults remain unchanged, and broader-panel validation is still needed.
