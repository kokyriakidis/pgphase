# Same-current-callset joint-gap controls

Question: can a joint solve over all eligible graph/BAM variants connect the
remaining gaps more naturally than separately solving and attaching blocks?

Export the accepted pgphase VCF (binary
`cb3b4ad62fb671abb1f375e02a48a44b1eac2553158b7ab39bc6c9755d94fc06`)
and original BAM over each gap with 100 kb context on each side. Preserve each
VCF record, genotype and allele representation; no sites or genotypes are
introduced or removed within that interval. HiPhase selects eligible input
variants with its own existing zygosity, GQ and variant-type rules. It ignores
incoming phase orientation and solves the input genotypes jointly. Run default
and `--disable-global-realignment` modes at MAPQ 5, one thread, with ignored
read groups. These are diagnostic competitor runs, not production changes.

| Gap | pgphase correct / phased (read PS) | Fresh HiPhase local on pgphase calls | Saved native-DV HiPhase |
|---|---:|---:|---:|
| 7,264,321–7,280,346 | 117 / 126 (3) | 120 / 122 (1) | 123 / 125 (1) |
| 24,121,713–24,131,707 | 60 / 60 (3) | 89 / 93 (1) | 93 / 93 (1) |

Both local-mode VCFs span the gap. For the spanning local phase sets, disjoint
10 kb flanks have consistent parental orientation and no discordant flank
reads: 42/42 left and 52/52 right at 7 Mb; 29/29 and 45/45 at 24 Mb.
`flank-votes.json` records these checks. Evaluation uses parental truth only
outside the solver. Local read correctness is oriented independently per PS,
and includes original primary read names overlapping the exact gap.

The default HiPhase mode on the same pgphase callset does **not** span either
pair: at 7 Mb it has 96/97 locally correct phased reads in two read blocks; at
24 Mb it has 59/79 in two blocks. Thus all-site access alone does not explain
success; representation-sensitive allele calling also changes the solve.
HiPhase local mode still performs its existing local allele alignment. This
control does not authorize or introduce realignment in pgphase and does not
isolate the solver on an identical read-by-site matrix.

The local runs take 1.02 s each after input export. Default runs take 1.16 s
and 3.34 s respectively during two concurrent window jobs. These are diagnostic
window timings, not chromosome runtime benchmarks.

The inspected HiPhase code uses `is_phasable_variant` in `src/block_gen.rs`,
then joint A* in `src/phaser.rs::solve_block` and
`src/astar_phaser.rs::astar_solver`. It filters variants and uses queue pruning;
it neither uses literally every catalog row nor promises an unrestricted exact
optimum. pgphase's graph solve, BAM subsolve and transfer preserve independent
gauges. Its joint MEC fallback is narrower: centered allele frequencies
(|AF−0.5| ≤0.12), category/provenance checks and a 20-variable resource bound.
These restrictions can remove intermediate informative sites from the solve.

The appropriate next experiment is a single gap observation matrix containing
all usable graph and BAM sites with allele identities and evidence channels
preserved, followed by a joint diploid solve. Gap variables may be rephased;
established flank blocks retain their internal orientations, and both possible
relative right-block gauges are evaluated. Connections should follow supported
read chains, rather than requiring one molecule across the entire gap. SNP
priority and uncertain/missing alleles must remain explicit. Compare methods
on the exact same observations before attributing the difference to a solver.

Reproduce the diagnostic export and HiPhase runs with
`compare_joint_calls.py` using the bench-phasers Python environment (pysam).
The script contains evaluation-only fixture paths and coordinates; outputs
are under `test_data/tmp_gap_next28/`. Score them with the existing
`evaluations/2026-10-02-padded-indel-context/compare_local_gap.py` and the
original input BAM/truth map. `metrics.json` retains the measurements.

No production behavior or test expectation changes in this investigation.
