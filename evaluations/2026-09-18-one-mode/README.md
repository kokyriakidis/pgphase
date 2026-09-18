# One mode: the hybrid

The catalog's sites phase, the alignment recovers the gaps. There is no second
configuration to choose between.

## What was removed

| removed | why |
|---|---|
| `collect-graph-variation --bam` and its recovery | it was a second mode for the same job; measured equal to the hybrid on one window and 118 correctly-phaseable reads short on the other |
| `--graph-first` / `--no-graph-first` | graph-first IS the mode now, so a flag for it is a flag for "the mode", and its negation is the second mode |
| the `nographfirst` window-test arm | nothing left to select |

`recover_unphased_windows_from_bam` stays -- it is the hybrid's own recovery
entry, and it keeps the `contig_name` and `allow_import` parameters because the
tid translation and the index-safety they exist for are real constraints, now
documented at one caller rather than two.

Behaviour is unchanged: on both panel windows the emitted VCF is **identical**
to the pre-consolidation default, byte for byte over the records. Window tests
66 assertions, unit 4/4.

## The efficiency question, measured rather than assumed

The obvious target was the alignment discovery the first solve appears to throw
away: under graph-first every alignment-derived candidate is withheld before the
catalog's sites are injected, so the discovery pass looks like wasted work.

It is not. Skipping `collect_candidate_sites_from_records` and
`collect_allele_counts_from_records` outright on
`chr20:25,979,591-26,138,679`, single-threaded:

| arm | wall time | candidates | phased hets |
|---|---:|---:|---:|
| as shipped | 4.05 s | 448 | **447** |
| discovery skipped | **2.78 s** | 137 | **137** |

31% faster for a collapse from 447 phased heterozygotes to 137. The candidates
are discarded but their by-product is not: `pre_process_noisy_regs_pgphase` and
the post-classification passes build `chunk.noisy_regions` from them, that model
scopes the noisy-region MSA, and the MSA supplies most of the chunk's phased
sites. Discovery feeds the phasing even though its candidate rows do not reach
the solve.

So the redesign is a consolidation, not a restructuring: discovery runs, its
candidates are withheld from the first solve, the catalog's sites phase, and the
alignment returns inside the windows the solve could not phase. The comment at
the withholding site records the measurement so the next reader does not retry
the same optimisation.

## What a real efficiency gain would need

Separating the noisy-region model from candidate discovery -- deriving
`noisy_regions` from the read digars alone, without building and classifying a
candidate table that is then dropped. `pre_process_noisy_regs_pgphase` already
works from reads; the post-classification refinement is what needs candidates.
Whether the refinement is load-bearing is the measurement to take, and it is not
taken here.

Not measured: chromosome-wide.
