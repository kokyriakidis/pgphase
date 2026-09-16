---
name: pgphase-gap-diagnosis
description: "Diagnose why pgphase leaves a phase-block gap unjoined, and whether a competitor closes it. Use this whenever the user asks about pgphase gaps, unresolved or unphased regions, `split`/`partial`/`rejected` gap-recovery verdicts, gap stitching or bridging, why a graph-only or hybrid run does not span a region, which gaps HiPhase/LongPhase/WhatsHap close that pgphase does not, or wants a gap-recovery change measured before it becomes a default. Also use it when validating any pgphase change that touches site admission, MSA verification, recovery tiers, or the MAPQ floor, since the read-level gate and truth-segregation controls here are what separate a real fix from a plausible one."
---

# Diagnosing a pgphase phase-block gap

A gap has exactly three possible causes, and they look identical in a tier
report. The whole point of this procedure is to tell them apart before anyone
proposes a fix:

1. **Nobody can close it** — no read crosses some point inside it. No
   site-admission policy can help; a competitor breaks there too.
2. **Admission policy** — the evidence is present, informative and
   read-chained, but a gate refuses the sites that would bridge it.
3. **Proposal structure** — sites are admitted and voting, but no single
   proposal block ever carries support on both flanks, so `joined` cannot
   become true however admission is widened.

Work in a dated directory under `evaluations/` (the repo convention) and keep
the runner script with the findings, so a result can be reproduced without
reconstructing the command.

## Step 0: is this gap even a deficit?

Never chase a gap no competitor closes. Of 311 chr20 pass-1 gaps (17.81 Mb),
only 187 gaps / 4.74 Mb are spanned by any of HiPhase, LongPhase or WhatsHap;
the other 13.07 Mb is spanned by nobody, and the largest gaps are
disproportionately in that class (mean 105 kb against 25 kb for the deficit
set). Two of the largest gaps absorbed a lot of investigation before this test
existed.

```sh
python3 evaluations/2026-09-16-competitor-deficit/find_deficit_gaps.py \
  --gaps evaluations/2026-09-15-graph-only-baseline/graph_only_gaps.tsv \
  --results-root "$EVAL_DATA/results/chr12-18-20-comparison/chr20" \
  --output <outdir>/deficit_gaps.tsv
```

A competitor "spanning" a gap means one of its phase sets starts at or before
the gap and ends at or after it. Note the weaker claim this test does *not*
make: a tool that spans nothing may still phase a large sub-block **inside** the
gap that we fragment further, so a nobody-spans gap can still hold recoverable
structure.

## Step 1: run the window with the audit export

```sh
./pgphase collect-hybrid-variation \
  --ref test_data/chm13v2.0.chr20.renamed.fa \
  --bam test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam \
  --graph-sites test_data/chr20.sites.striped.vcf.gz \
  --gaf test_data/HG002.chr20.annotated.coord.gaf.gz \
  -r "CHM13#0#chr20:$((GL-50000))-$((GR+50000))" -t 8 -q 1 \
  --link-by-alleles --block-link-window 8 --min-read-margin 2 \
  --recover-gaps --gap-decision-audit $D --gap-recovery-report $D/tiers.tsv \
  -o $D/candidates.tsv --phased-vcf-out $D/native.vcf -b $D/phased.bam
```

`--gap-decision-audit` is the instrument for every evidence question: it writes
`tid*.evidence.tsv` with `SITE` rows (position, type, category, `msa_verified`,
`is_homopolymer_indel`), `READ` rows, and per-read `OBS` rows carrying both the
BAM and graph alleles, plus `votes.tsv` and `reads.tsv` per tier.

**The audit and a normal run do not evaluate the same tiers.** The homopolymer
tier is reached in recovery pass 1 only when the audit is on, and
`stitch_gap_proposal` is called orientation-only whenever `audit`,
`recovery_pass == 1`, or `tier == kGapHomopolymerTier`. Identical final blocks
therefore do **not** prove the two runs agree — confirm any verdict you intend
to act on with a normal run.

## Step 2: diagnose

```sh
python3 evaluations/2026-09-16-competitor-deficit/diagnose_gap.py \
  --run-dir $D --gap-left $GL --gap-right $GR \
  --truth-map /tmp/truth_hap.tsv \
  --competitor-vcf "$HP/phased.vcf.gz" \
  --competitor-per-read "$EVAL/hiphase_reads/per_read.tsv.gz" \
  --output <outdir>/gap_${GL}_break_sites.tsv
```

It reports, in the order the reasoning needs them:

- **The residual break.** Our blocks rarely leave the whole gap open; the part
  that matters is the interval no single phase set covers. A 78 kb gap can carry
  a 19 kb break.
- **The competitor's site inventory against ours**, matched exactly and within
  15 bp. Indel anchor placement differs between callers, so exact-position
  matching alone overstates how many sites we miss.
- **Per-site truth segregation** for every site in the break, with its category,
  flags, and which recovery tier may admit it. This is the column that decides
  whether an exclusion is costing anything.
- **Read spanning across the break chain** — the test that separates cause 1
  from causes 2 and 3.
- **The competitor's own accuracy** on the reads overlapping our break. A join
  scoring 95% is a different proposition from one scoring 99.9%; do not chase a
  wrong join.

## Step 3: classify the blocker

`classify_blockers.py` (in `evaluations/2026-09-15-gap-link-site-gate/`) runs the
chain test across many gaps at once from the alignment alone, needing neither
truth nor the audit — useful for deciding which gaps deserve a drilldown. Its
`chain_present` verdict means only that coverage does not rule a join out; it is
a necessary condition, not evidence that a gate fix will produce one.

## Step 4: gate any change before believing it

Two arms, matched by read name against the read-level truth
(`compare_panel.py` in `evaluations/2026-09-15-gap-link-site-gate/` does this
for a whole window panel):

- **concordant → discordant must be 0.** A join that flips correct reads is
  worse than the gap.
- **Report concordant tags lost.** A join can be correct and still cost
  coverage: one verified fix closed its gap while 23 correct read tags went
  untagged, another closed a larger gap with zero reads changing label at all.
  Both are real outcomes; the difference only shows up if you count.
- **Run the panel, not just the motivating gap.** Deriving a change from the one
  gap that motivated it has overstated it more than once here. A fix that closes
  its own gap and is inert on eight others is a narrow fix, and worth saying so.

## Gotchas that have each cost a wrong conclusion

- **Phase-set extent is first-to-last phased variant.** Measuring from read
  start positions overstates a block's reach by up to a read length per end and
  understated a residual break by an order of magnitude.
- **A break interval is not centred on the gap.** Carry its real endpoints; a
  midpoint assumption counted sites in a window that did not overlap the break.
- **`msa_verified` is not informativeness.** In one 21 kb break all four
  verified sites carried `msa_verified = 1` and segregated at 0.509 (chance),
  0.607, 0.644 and 0.923.
- **Always include a trusted-class control.** Score `CLEAN_HET_INDEL` sites
  alongside whatever you are testing and require them to come out informative.
  This control caught two separate genotyping bugs that had produced a plausible
  negative result.
- **Genotype repeat-context indels from the net insert-minus-delete length over
  a window, never at the VCF anchor position** — the aligner places them
  arbitrarily within their repeat run.
- **Allele encoding in pgphase variant TSVs** (`collect_output.cpp`): a deletion
  is the deleted bases in `REF` with `ALT = .`; an insertion is the inserted
  sequence in `ALT` against a single anchor base. `len(ALT) - len(REF)` is not
  the length change and silently skips every deletion.
- **The emitted graph BAM is unaligned.** Region queries against it fail
  silently and every read looks untagged; read `phase_reads.tsv` instead. For
  the same reason the evaluator's `phase_block_n50_bp` / `phase_block_aun_bp`
  are invalid on that path — one reported auN of 15.47 Mb against a largest
  block of 1.28 Mb, which is impossible. Derive block spans from
  `phase_sites.tsv`.
- **`csv.writer` defaults to CRLF.** Pass `lineterminator='\n'`, or a shell
  runner silently matches zero windows because the last field reads as `gap\r`.
- **An absence claim is only as strong as the input inventory it was computed
  over.** Before concluding a region has no evidence, enumerate sites by
  disposition — retained by category *and* each filtered reason — not just what
  survived our filters.
