# Why gap-by-gap recovery is not generalizing

Investigation of the current executable
`c6bcf86a43b5301cd90960eeb8b525d14944494dfb8381f9b676c5641952e674`.
No production behavior, executable, test expectation or certified closure changes.
The finding is an integration and representation problem, not simply absence
of graph/BAM information. A joint binary solver already exists; the data it
receives and the decisions surrounding it remain fragmented.

## Current inventory

Deduplicate the latest regression contract by both gap endpoints. Audit every
window against `test_data/tmp_gap_fix91/final_current/0`, original primary BAM
overlaps and the existing globally oriented HiPhase tag cache. Assert all 141
HiPhase denominator/total/core measurements against the committed benchmarks.
The 141 windows contain 18,774 distinct truth-scorable primary reads.

Both native regressions and the full chromosome have 65 passing / 76 failing
windows under spanning, >=80% all-overlap correctness, and total/core HiPhase
parity. Their failure categories differ at six windows because owning context
and chromosome stitching differ; no short replay should certify a global join.
Of the 118 windows where HiPhase achieves 80%, 59 pass and 59 fail.

The remaining 59 on the full chromosome partition as follows:

| First failing condition | Windows |
| --- | ---: |
| No spanning variant block | 5 |
| Spans, but below 80% correct | 7 |
| Spans and >=80%, but fewer total correct than HiPhase | 30 |
| Total correctness passes, but dominant connected core is short | 17 |

Thus 54/59 are already spanned. Extending VCF blocks alone cannot address the
majority. In 50/59, every candidate strictly inside the gap is repeat/noisy;
six have no interior candidate and eleven have only unphased interior candidates.
These last three counts overlap and are not exclusive failure mechanisms.
`panel-audit.json` records every measurement and `native-gap-contract.tsv`
preserves the original regression report, including repeated identical rows.

## What the implementation actually does

1. `graph_collect.cpp:10693` first solves the graph/GAF matrix with
   `kCandGermlineClean`. Repeat indels have already been demoted.
2. `collect_pipeline.cpp:786` nominates seams from consecutive oriented sites
   whose phase-set lists have no identity in common. A spanning PS with poor
   read coverage is not itself a seam nomination.
3. `recover_phase_set_seams_in_place` runs separate BAM/MSA solves, selects
   source phase sets, translates graph-walk keys to sequence keys, imports
   candidates and records gauges/cuts. It does retain shared observations;
   this is not a claim that all BAM evidence is discarded.
4. `trusted_recovery_mec_path` already uses an exact diploid MEC kernel, but
   inside a guarded adjacent-block fallback. Its ordinary eligibility uses
   AF within 0.12 of 0.5, with several later exceptions, and a 20-variable
   bound. Source-path, existing block, and boundary certificate checks remain
   outside that optimization.
5. The whole-chunk BAM overlay adds exact matched observations to existing
   biallelic graph profiles. Private whole-chunk sites are not generally added
   by this overlay; independently phased BAM reads are staged as fallback.
6. After stitching, rescue fills only unassigned reads, graph first and BAM
   second. `graph_bam_adapter.cpp:2626` skips already tagged core/rescue reads;
   successful excluded-site rescues receive PS + 1,000,000,000 without changing
   candidates or core connectivity.
7. The batch then runs the separate insertion, source-rescue, masked deletion,
   complementary repeat, repeat-SNP, terminal, compound and tandem recovery
   passes at `graph_collect.cpp:10754`. Their checks are reusable code but
   specific to variant shapes, source topology and bounded neighborhoods.

Two consequences matter. A BAM source PS is a provisional clustering result,
yet several acceptance paths require its entire gauge/cut history to fit a
specific pattern. Separately, reads can have useful gap observations and a
correct rescue HP while their PS remains outside the core. Adding a new
variant pattern to one late pass does not create a common inference problem
for the other patterns.

## Fresh evidence from three owning chunks

Replayed 4–5, 5–6 and 6–7 Mb with the unchanged executable, recording complete
graph/recovery/source/BAM-overlay matrices. The native source state is in
`test_data/tmp_joint_evidence_investigation/{4,5,6}`.

At **5,309,406–5,345,085**, pgphase has **206 correct / 1 discordant / 19
unphased**, but only **140 correct in the dominant core**. HiPhase has
**220 / 1 / 5**, all 220 in one core. Among the 80 HiPhase-correct reads absent
from our correct core, all have a retained profile, none has a direct SNP
observation in the established core, and 75 lie in a callable read/site
component that reaches a core SNP. Repeat/noisy sites must carry their
connection. This is optimistic connectivity, not an independently validated
phase certificate: false, correlated or ambiguous observations can form edges.

At **6,513,891–6,516,221**, the two HiPhase-correct reads misassigned by pgphase
have a wrong graph SNP vote but a correct BAM/effective SNP vote in the current
matrix. Simply dropping graph/BAM disagreements abstains on both. Existing HP
labels survive because rescue is additive and the physical refresh requires
additional spaced/verified witnesses. This is evidence that the observations
can improve these reads, not a license to always override graph calls with BAM.

The **4,766,928–4,792,960** owning replay also exposes context sensitivity: it
has 154 correct core reads, already matching HiPhase, but 14 discordant reads.
Its short native contract has 156 total correct and only 112 core correct.
Fixing only its short-window representation would conceal the wider gauge
problem. `matrix-audit.json` lists individual read evidence for all three controls.

## Controlled HiPhase comparison: the missing information is partly present

HiPhase 1.6.0-ac3f399 uses exactly the same original alignments overlapping
5.25–5.40 Mb, with their full sequences/CIGARs retained. Score the same 226
original primary truth-scorable overlaps with 5,309,406–5,345,085. All input
GTs are unphased and PS is cleared. Existing pgphase VCF fields are preserved
in the surgical controls. Alternate competitor alleles are evaluation-only;
no truth or HiPhase output feeds production.

| HiPhase input/control | Correct | Discordant | Unphased | Dominant core correct |
| --- | ---: | ---: | ---: | ---: |
| Current emitted pgphase sites, global realignment | 199 | 2 | 25 | 126 |
| Same sites, global realignment disabled | 199 | 2 | 25 | 126 |
| Add only 5,315,591 C>CT | 208 | 8 | 10 | 208 |
| Replace only the compound repeat by one ALT/ALT contrast | 195 | 1 | 30 | 122 |
| Both marker and ALT/ALT contrast | 215 | 1 | 10 | 215 |
| Competitor DeepVariant sites | 220 | 1 | 5 | 220 |

The current graph catalog **already contains C>CT,CTT at 5,315,591**. The
selected +T candidate is phased internally into PS 5,239,261, classified
`REP_HET_INDEL`, and omitted from the emitted VCF. The working matrix already
has **43** left-SNP/+T observation pairs and **five** +T/downstream-deletion
pairs, with both diploid classes represented: three (REF,ALT), two (ALT,REF).
All five latter pairs use the BAM channel; the graph channel alone has zero.
HiPhase with the marker has 43 and six pairs respectively and connects the
blocks. This separates missing published site records, missing graph-channel
calls and an available combined read chain.

The compound repeat illustrates a different representation problem. Pgphase
publishes two independent insertion rows at 5,339,368, each contrasted against
genomic reference. DeepVariant supplies one multiallelic row at 5,339,363 with
GT 1/2 (written in the input's existing order), contrasting the two sample ALT
sequences. These are differently positioned descriptions, not interchangeable
binary observations. The combined control improves another seven correct
reads and removes seven discordant calls relative to the marker-only control.
It still falls five correct reads short of the complete competitor input;
this experiment does not attribute that remainder to either tested edit.

The extra marker alone supplies connectivity but does not supply accurate
allele assignment everywhere. The contrast alone improves repeat calling but
does not supply the missing upstream bridge. Both are required in this control.
Global realignment alone cannot recover a variant omitted from the input.
`competitor-controls.json` preserves counts, input hashes and named-segment
pair counts. HiPhase solves took 0.96–3.49 seconds on this bounded input.

A preliminary control normalized all existing pgphase VCF fields and gave a
marker-only 210/1/15 result. That was not the surgical comparison. The retained
table above uses preserved existing fields (208/8/10); only these final controls
are used for the conclusion.

## Does admitting every site fix the current solver?

Extract the current `solve_recovery_mec` implementation unchanged, compile a
small standalone driver and replay each current matrix. Keep emitted clean SNP
gauges in the dominant block fixed. Treat callable binary non-SNP sites between
the nearest clean flanks as variables. Compare AF-centered versus all such
variables, effective observations versus an agreement-only graph/BAM union,
and full reads versus deterministic FNV read halves. Compare uniform mismatch
costs with strict SNP priority (one fixed SNP mismatch outweighs every variable
mismatch in the problem). These are controlled kernel inputs, not invocations
of the production eligibility/stitch wrapper. Actual nonbinary sites abstain.
The binary model still treats decomposed rows as separate factors; it
does not provide a new complete-repeat allele likelihood or a production certificate.

| Current fixture | Variables, all sites | Current correct | Uniform all-site correct | SNP-priority all-site correct | HiPhase correct |
| --- | ---: | ---: | ---: | ---: | ---: |
| 4,766,928–4,792,960 | 11 | 154 | 153 | 153 | 154 |
| 5,309,406–5,345,085 | 7 | 206 | 201 | 202 | 220 |
| 6,513,891–6,516,221 | 3 | 86 | 93 | 88 | 88 |

The diagnostic predictions include no minimum-confidence read gate. They are
not emitted and do not establish accepted closures. The all-site 6.514 Mb
model gains net correctness but also loses five previously correct core
assignments. Strict SNP priority at 6.514 Mb instead reaches HiPhase's 88 correct
without losing any old correct core assignment. The centered two-variable
arm gives that same result in both halves; the all-site arm abstains on two
previously correct reads in half 0. At 5.309 Mb, the two read halves disagree
on one variable and neither cost scheme reaches parity. An
exact optimum under a lossy binary objective is not proof of an accurate phase.

The current kernel calls take approximately **2–5 ms**, including subprocess
startup, on these 2–11-variable problems. Replaying the same saved observation
state is a viable fast development loop. Adding observations/variables changes
the objective, so scores across centered/all arms are not directly comparable.
`joint-replay.json` and `joint-replay-snp-priority.json` record all states,
predictions and the extracted kernel hash.

The older `2026-10-02-joint-matrix-investigation` reached the same caution on
7/24 Mb observations: admission alone did not repair contradictory SNP calls
or missing complementary-repeat bridges. It is historical evidence, not a
measurement of this executable.

## Follow-up: duplicates and MSA representation

`audit_representation.py` checks the same full chromosome VCF and the three
saved overlay matrices. Reference normalization with bcftools 1.24 validates
every REF, realigns 504 of 64,483 records, and exposes **10 duplicate allele
pairs** despite zero exact duplicates in the unnormalized output. All ten
pairs agree on phased GT and PS; three differ in AD. For example, insertions
at 8,400,078 and 8,400,084 normalize to the identical allele at 8,400,078.
This establishes redundant output descriptions, not their upstream origin
or an independent-vote phasing defect. Do not sum their depths: shared original
molecules must count once. `representation-audit.json` records every pair.
The audit retains original-record tags because normalization can reorder rows.

The 4/5/6 Mb matrices have 25/47/24 raw candidate-key collision groups
(30/56/29 extra rows). These are not necessarily duplicate alleles. Graph
construction copies the site ID into `CandidateVariant.key.alt` and retains it
when decomposing alternatives; the selected allele lives in parallel metadata.
At 5,331,128, candidates 461 and 462 share an insertion key naming
`>115117600>115117608`, but that catalog site offers multiple distinct A-repeat
lengths. Deleting rows by this key would discard distinct alternatives.
Recovery's raw and translated sequence indexes use single-index `map::emplace`
(`collect_pipeline.cpp:3206–3214`), so collisions retain the first entry.
This is an identity/ambiguity hazard; the audit does not establish that these
particular raw collisions caused a wrong BAM match or independently weighted
vote. Exact minimal-VCF cross-snarl deduplication already exists, but it does
not left-align repeat-shifted descriptions.

MSA is a useful source of sample allele sequences, not blanket certification.
`make_noisy_candidate` sets `msa_verified=true` while constructing a candidate,
initially with zero allele counts and `NoisyCandHom` category
(`collect_phase_noisy.cpp:175–200`). The flag records MSA-derived provenance;
it does not independently certify genotype, per-read assignment, or phase.
`collect_noisy_reg_aln_strs` can also build consensuses using the current
read HP/PS grouping (`align.cpp:1970–2043`), so their orientation is not
independent evidence for that grouping.

The recommended first change is a common reference-anchored locus/allele
identity shared by graph and BAM/MSA, with repeat normalization and complete
sequence comparison for compound alleles. Merge equivalent descriptions and
retain provenance/path context; represent distinct sample alternatives as one
diploid locus. Preserve MSA consensus windows and per-read sequence confidence,
and rebuild counts from unique original molecules. Prefer a supported MSA
sequence when resolving an equivalent graph description, but infer genotype
and phase from joint evidence instead of importing the source HP as truth.
The two insertion rows at 5,339,368 describe distinct ALT alleles and should
become a joint ALT1/ALT2 contrast, not be deleted as duplicates. The controlled
HiPhase evidence above supports repairing that contrast as well as retaining
its missing upstream marker. A single-index collision must be resolved by
allele identity or treated as ambiguous, rather than selecting the first row.
Test duplicate-input and input-order invariance on cached states, then compare
all 59 HiPhase-correct targets under the existing parental/80%/total/core gates.
Production code, output and certified expectations remain unchanged.

Reproduce this follow-up independently:

```bash
/home/kokyriakidis/micromamba/envs/bench-phasers/bin/python \
  evaluations/2026-10-06-joint-evidence-investigation/audit_representation.py \
  > evaluations/2026-10-06-joint-evidence-investigation/representation-audit.json
```

## Proposed general change

Build a shared component inference stage over graph and gap evidence, rather
than routing each new allele shape through another late promotion pass:

1. **Keep complete locus alleles.** Recover the two sample sequence alleles
   using graph branches plus BAM/MSA consensuses. Retain parent/child snarl
   conditions and deletions masking SNPs. Overlapping/split records project
   from one locus variable; absence of one ALT is not automatically literal REF.
2. **Keep uncertainty per original molecule.** Store likelihood/cost for each
   compatible allele, with ambiguous calls remaining ambiguous. Graph and BAM
   observations of the same molecule/locus are dependent evidence, not two
   independent votes. Use physical sequence/quality to calibrate disagreements
   instead of blindly selecting a channel or trusting its early HP label.
3. **Solve sites, block orientation and read HP together in each component.**
   Verified flank relationships constrain the problem; conflicting calls and
   provisional BAM source PS do not become immutable haplotype facts. Use all
   usable intermediate sites with appropriate confidence. A site and a read
   should earn core membership from the same validated allele chain, rather
   than obtaining a correct HP but a permanent output-only rescue PS.
4. **Use connected evidence to nominate work.** Include unassigned/rescue reads
   and ambiguous loci inside an already spanning block. Physical allele
   disconnection, not just adjacent VCF PS identity, determines a boundary.
   Keep unsupported components separate.
5. **Cache evidence before optimization.** Persist immutable component
   observations and likelihoods with identities covering inputs, representations
   and extraction behavior. Invalidate solver results separately. Test generic
   mechanisms on saved states in milliseconds, then replay affected owning
   contexts and run the complete panel after a stable change. Frontier state or
   bounded search should scale by locally active uncertainty instead of total
   sites in a large block; its limits must be explicit and measured.

Start with the existing MEC kernel and matrix dumps. The first implementation
milestone should be a generic complete-locus evidence adapter and shared
component read assignment, evaluated together across **all 59** current
HiPhase-correct targets, not another dedicated marker rule. Contrast repaired
observations under the same optimizer before changing optimizer behavior.
Keep the >=80% all-primary and total/core HiPhase parity contract, parental
orientation checks, stronger existing floors and full-chromosome side-effect
audit. Parental truth remains evaluation-only. These results support this
direction; they do not prove the proposed model will close all targets.

## Reproduction

From the repository root, with pysam and the frozen competitor available:

```bash
BENCH_PYTHON=/home/kokyriakidis/micromamba/envs/bench-phasers/bin/python \
HIPHASE=/home/kokyriakidis/micromamba/envs/bench-phasers/bin/hiphase \
bash evaluations/2026-10-06-joint-evidence-investigation/replay.sh
```

For warm investigation, rerun just `audit.py`, `audit_matrices.py`,
`replay_joint.py` and `audit_competitor.py`. `manifest.json` fingerprints the
current executable, relevant outputs, controlled inputs and source identities.
The reference/BAM/GAF are identified by metadata to avoid rehashing gigabytes
for a diagnostic replay. All three owning runs, all six HiPhase controls and
all 72 exact-kernel calls completed successfully. Python/shell syntax checks
pass. Production code was not edited, so a production build/unit rerun is not
needed for this investigation. Existing uncommitted gap fixes are preserved.
