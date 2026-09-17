# In the gap, discard the graph's sites and keep only verified alignment sites

`--gap-bam-only` (with `--retry-unphased-with-bam`, which is what detects and
re-solves an unphased window) discards every graph-catalog candidate inside that
window and keeps only the alignment channel's `msa_verified` ones, so the window
is re-solved on verified alignment evidence alone.

```sh
OUT=/tmp/panel-gapbam FLAGS="--retry-unphased-with-bam --gap-bam-only" \
  ./evaluations/2026-09-16-test-panel/run_panel.sh
```

## Result: a clean null

The two arms are **byte-identical in all six panel windows**.

| | retry | retry + `--gap-bam-only` |
|---|---:|---:|
| windows spanned | **3 of 6** | **3 of 6** |
| in-gap phased hets | 31 | 31 |
| reads tagged | 3,149 | 3,149 |
| concordance | 99.49% | 99.49% |
| gate: concordant -> discordant | **0** | **0** |
| new concordant / new discordant | 348 / 7 | 348 / 7 |

And it is not because there is nothing to discard. In the target window the
detector reports:

```
retry window 55,774,000-55,824,000 (50.0 kb): 421 candidates inside, 418 graph sites
retry window 55,843,827-55,889,114 (45.3 kb): 613 candidates inside, 600 graph sites
retry: 3 unphased window(s), 0 site(s) admitted
```

**600 of the 613 candidates in the gap are the graph's, and removing all 600
changes nothing.** They are in that window precisely because they failed to phase
it: they contribute no usable het or link evidence, so their presence is not
crowding anything out.

That rules out a hypothesis worth ruling out -- the gap is not unphased because
the catalog's sites outnumber and drown the alignment channel's. It is unphased
because the alignment channel's own sites are not usable there, which is the
`skip_noisy_kmeans` exclusion measured in the previous section and, beneath that,
the per-read allele call at indels.

## What did change: the retry, now that the representation is right

Worth recording because the same flag was measured before this work and behaved
differently. The retry arm on the panel:

| retry arm | spanned | in-gap hets | new concordant | new discordant | c -> d |
|---|---:|---:|---:|---:|---:|
| before the representation fixes | 0 of 6 | 19 | 333 | 9 | **1** |
| after | **3 of 6** | **31** | 348 | **7** | **0** |

Three windows now close, with no read flipped from correct to incorrect anywhere
in the panel. The fixes that moved it are the ones in
`evaluations/2026-09-17-multiallelic/`: merging a co-located pair into one record
with both alleles, the event-aware reference test, and the joint two-haplotype
orientation.


## Retired

The flag is removed (its measurement stands, above). It was an exact null: in all
six panel windows, discarding the graph catalog's candidates inside an unphased
window and keeping only the MSA-verified alignment ones produced output
byte-identical to the plain retry -- and not for want of material, since the
great majority of candidates inside those windows are graph sites and removing
all of them changed nothing.

That result is worth keeping and the code is not. What it cost to keep was four
pieces of CLI surface, an `Options` field, and two conditionals inside the
window-readmission loop, all of which a future reader has to reason about when
changing that loop. The loop now reads as what it does with the flag off, which
is the only way it ever ran: a graph site keeps whatever category the recovery
pass gave it, and only the alignment channel's candidates are readmitted.

Removal verified behaviour-neutral: panel stock defaults and
`--retry-unphased-with-bam` both emit records identical to the build before the
removal -- 0 lost, 0 gained -- suite 5/5, no references left in the tree, and
`--gap-bam-only` is now rejected as an unrecognized option.
