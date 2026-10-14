# AGENTS.md

Project-specific guidelines for AI agents working on pgphase.

## Project Overview

pgphase is a C++17/Rust tool for variant calling and haplotype phasing from long reads. It has three commands:

- `collect-bam-variation` — SNP/indel calling from BAM/CRAM alignments
- `collect-graph-variation` — variant calling from pangenome graphs (GBZ + GAF)
- `build-snarl-catalog` — preprocess GBZ graphs into phasing site catalogs

## Build System

- C++17 compiled with GCC, linked against htslib, zlib, pthreads.
- Third-party libraries: WFA2-lib, abPOA, edlib (git submodules under `third_party/`).
- Rust FFI component: `src/gbz_ffi/` (static library linked into the C++ binary).
- Rust tools: `third_party/gbz-base/` (query, gaf2db, gbz2db binaries).
- Build order: `make third-party-libs` → `make gbz-base` → `make -j$(nproc)`.

## Development loop

One entry point, `scripts/dev/pg` (Makefile aliases in brackets). Paths come
from `scripts/dev/env.sh`; put machine-specific ones (truth VCF, competitor
outputs, a Python with pysam, cache dir) in the untracked
`scripts/dev/local.env`.

- `pg test` [`make dev`] — parallel build of pgphase and all unit binaries,
  then every unit and Catch2 predicate test; one line per binary. A one-file
  edit to a union stage is ~7 s from save to green.
- `pg quick` [`make quick`] — `--bam` on the ten 1 Mb arm windows in
  `scripts/dev/quick_regions.bed` (~30 s), scored and tabled against the pinned
  baseline and the competitors (HiPhase, longphase) on the same reads.
- `pg full` [`make full`] — the whole of chr20 (~2 min), same report.
- `pg same [quick|full]` [`make same`] — output identity (candidates.tsv,
  phased.vcf body, every read's HP/PS) with the baseline. Required for any
  refactor; add `-- --graph-only` for the graph-only path.
- `pg diff [quick|full]` [`make diff`] — arms_phaseable reads whose class
  changed against the baseline, by transition and by 1 Mb bin, with a TSV of
  the reads for drill-down.
- `pg pin [quick|full]` [`make pin`] — make the current binary the baseline,
  once a change is accepted.

Runs are cached by binary hash and arguments under `PGDEV_CACHE`, so
re-scoring an unchanged binary is instant. The headline metric, checked every
iteration, is **arms_phaseable correct**: reads on the chromosome arms (outside
26-32 Mb) whose span crosses a truth heterozygote, phased and consistent with
their phase set's majority parental orientation. Watch its discordant and
unphased counts and arm switches alongside. Diagnostics, off by default:
`PGPHASE_TIMING=1` prints each chunk's seconds per stage (and the alignment
solve's load / classify / phase split) to stderr; `PGPHASE_EM_DUMP=DIR` writes
each chunk's final EM evidence for scripts/dev/why.py, locus.py and
proto_em.py; `PGPHASE_DIAG=DIR` writes the alignment solve's candidates
before pruning, its noisy regions, every injection decision, the graph rows,
the final EM's voting candidates, every locus window's fate (windows.*),
every read's profile as the final EM saw it (emreads.*), every EM block-cut
decision with its margins (joins.*), every noisy-region MSA variant with its
cluster consensus coverage and read calls before and after the refresh of
assigned reads (msa.*), and each read's calls in the alignment solve
(alnreads.*), every read an island block claims (islands.*), and every length
window with its read sides (lengthwin.*). `scripts/dev/pg diag` makes a cached
full-contig run with all of these; `scripts/dev/pg inspect REGION` prints
everything at a region by position (truth, catalog, graph, alignment,
injection, EM, our and every competitor's VCF, the reads of each tool, the
MSA calls and the EM's joins); `scripts/dev/pg reads REGION` shows every read
at a region against every truth het, rival site and EM site, its allele in the
alignment solve and in the EM against its truth haplotype;
`scripts/dev/pg gapsites` lists, per gap a competitor closes, the sites it
uses there (snp / indel / str / homopolymer, truth) and how far each gets
with us and why;
`scripts/dev/pg each [-j N] CMD` runs any pg command on every profile in
`PGDEV_EACH` (default chr20 chr12 chr18), each in its own environment;
`scripts/dev/pg sweep LABEL[:VAR=VAL,...] ...` runs the full tier per
configuration and prints each profile's change from its pinned baseline and
the total (a diag run doubles as the full run of its binary and environment,
so neither runs twice); `scripts/dev/pg islands` reports the island blocks
(classes, truth, purity and production features against truth);
`scripts/dev/pg calls` the EM's per-read call accuracy at rejected sites;
`scripts/dev/pg hiphase-exact` HiPhase's own allele assignment (a line-for-line
port, scripts/dev/hiphase_port.py, checked by test_hiphase_port.py) on its own
inputs, read by read against our EM's calls; `scripts/dev/pg missing` where
each real het without an EM site is lost; `scripts/dev/pg lengthwin` the length
windows against truth;
`scripts/dev/pg unphased` every arms_phaseable read we leave unphased in one
primary cause, crossed with the rivals' classes, with the loci it concentrates
in and how many of each locus's reads each rival phases correctly;
`scripts/dev/pg noise [--rival NAME]` the EM-reliable sites the parents do not
split and the discordant reads they label, by production feature;
`scripts/dev/pg switches [--all]` every arm switch (and flip) of the run against
the baseline, new, kept or gone, between which truth-matched hets;
`scripts/dev/pg audit [RIVAL] [--reads islands|FILE]` puts every read a competitor phases correctly
and we do not -- unphased or discordant -- in one primary cause, and splits
the discordant ones by their own calls judged by each parent's side at the
site (no allele representation compared); scripts/dev/gaptriage.py classes
gaps, scripts/dev/chancheck.py compares the EM's calls with the alignment
solve's against truth;
`scripts/dev/pg gaps` lists the gaps we leave that a competitor closes (by
locus: our flanking blocks, the competitor's phase set, our reason per read);
`scripts/dev/pg funnel [--rival NAME] [--support N]` traces where the reads a
competitor phases correctly and we do not are lost (by truth het, or by the
competitor's own phased sites). `PGDEV_PROFILE=chr12|chr18` evaluates another
chromosome (full tier). `scripts/dev/score.py` is the one
scorer (`run`, `table`, `diff`, `bins`); parental truth is evaluation-only.

## Testing

- `make check` — validation gates (requires `test_data/`).
- `make unit-tests` — builds and runs unit test binaries.
- `make gap-dev-check` — fast development loop over in-memory phasing fixtures
  (allele, adapter and Catch2 predicate tests); no production binary build.
  `PREDICATE="[deletion]"` optionally selects Catch2 predicates.
- `make predicate-tests`, `make parity-tests`, `make upstream-parity-tests` —
  Catch2 phasing predicates, port parity, and parity against a local longcallD
  checkout (`LONGCALLD_ROOT`).

## Key Files

| Path | Purpose |
|---|---|
| `CHECKPOINT.md` | Evaluation results, design decisions, noise filtering rationale |
| `src/noise_filter.hpp/cpp` | Shared noise detection (BAM + graph pipelines) |
| `src/graph_bam_adapter.hpp/cpp` | Graph-to-phasing bridge, noise filter wiring |
| `src/graph_collect.cpp` | Graph pipeline worker loops |
| `src/union_*.cpp` | `--bam` union gap phasing, one file per stage: `union_inject` (alignment sites into the chunk), `union_reads` (alignment-only reads, fills, realignment), `union_em` (global EM, block cuts), `union_windows` (haplotype and allele windows), `union_phase` (driver, weak-bridge sites); shared helpers in `union_internal.hpp` |
| `scripts/dev/` | Development loop: `pg`, `score.py`, `identity.py`, quick-tier regions |
| `src/collect_var.cpp` | BAM pipeline candidate classification |
| `src/collect_phase.cpp` | K-means phasing, category bitmasks |
| `src/phasing_types.hpp` | Core types: VariantKey, Interval, constants |
| `benchmark/` | Snakemake benchmark pipeline |
| `scripts/` | Evaluation and benchmark helper scripts |

## Architecture Notes

### Two pipelines, shared phasing core

Both pipelines produce a `PhasingChunk` with `CandidateVariant` entries, then share the same k-means phasing and chunk stitching code. The graph pipeline converts snarl site observations into candidates via `build_graph_chunk` in `graph_bam_adapter.cpp`.

### Noise filtering

Reference-based noise detection (homopolymer, tandem repeat, SDUST low-complexity) is shared between pipelines via `noise_filter.hpp/cpp`. The BAM pipeline has additional read-level noise filtering (XID clustering, noisy-reads-ratio, MSA recall) that depends on CIGAR data unavailable in the graph pipeline. See `CHECKPOINT.md` "Noise Filtering: BAM vs Graph Pipeline" for the full rationale.

### Category bitmask system

Candidates are classified with both a `VariantCategory` enum and a `lcd_var_i_to_cate` bitmask. The bitmask controls which candidates enter k-means phasing:

- `kCandGermlineClean` (SNP | Indel | Hom) — used for graph pipeline k-means.
- `kLongcalldRepHetVar` (0x010) — repeat/homopolymer indels, excluded from k-means.
- `kCandGermlineVarCate` — includes noisy-recalled candidates, used in BAM pipeline second k-means.

### Thread safety

`faidx_t` is not thread-safe. Each worker thread opens its own handle via `load_reference_index(opts.ref_fasta)`. The graph pipeline does this in both `process_graph_chunk_batch` and `process_graph_chunk_batch_indexed_gaf`.

---

## Behavioral Guidelines

Imported from [andrej-karpathy-skills](https://github.com/multica-ai/andrej-karpathy-skills/blob/main/CLAUDE.md).

### 1. Think Before Coding

Don't assume. Don't hide confusion. Surface tradeoffs.

- State assumptions explicitly. If uncertain, ask.
- If multiple interpretations exist, present them — don't pick silently.
- If a simpler approach exists, say so. Push back when warranted.
- If something is unclear, stop. Name what's confusing. Ask.

### 2. Simplicity First

Minimum code that solves the problem. Nothing speculative.

- No features beyond what was asked.
- No abstractions for single-use code.
- No "flexibility" or "configurability" that wasn't requested.
- No error handling for impossible scenarios.
- If you write 200 lines and it could be 50, rewrite it.

Ask yourself: "Would a senior engineer say this is overcomplicated?" If yes, simplify.

### 3. Surgical Changes

Touch only what you must. Clean up only your own mess.

When editing existing code:

- Don't "improve" adjacent code, comments, or formatting.
- Don't refactor things that aren't broken.
- Match existing style, even if you'd do it differently.
- If you notice unrelated dead code, mention it — don't delete it.

When your changes create orphans:

- Remove imports/variables/functions that YOUR changes made unused.
- Don't remove pre-existing dead code unless asked.

The test: Every changed line should trace directly to the user's request.

### 4. Goal-Driven Execution

Define success criteria. Loop until verified.

Transform tasks into verifiable goals:

- "Add validation" → "Write tests for invalid inputs, then make them pass"
- "Fix the bug" → "Write a test that reproduces it, then make it pass"
- "Refactor X" → "Ensure tests pass before and after"

For multi-step tasks, state a brief plan:

```
1. [Step] → verify: [check]
2. [Step] → verify: [check]
3. [Step] → verify: [check]
```

---

## C++ Conventions

Adapted from [XOOS C++ rules](https://github.com/Roche-DIA-RDS-CSI/XOOS).

### Naming

| Element | Convention | Example |
|---|---|---|
| Files | snake_case | `noise_filter.cpp`, `graph_sites.hpp` |
| Variables | snake_case | `ref_beg`, `chunk_rows` |
| Private members | trailing underscore | `loaded_tid_`, `seq_` |
| Functions | snake_case | `find_low_complexity_intervals()` |
| Classes/structs | PascalCase | `GraphSiteMeta`, `PhasingChunk` |
| Constants | kCamelCase | `kSdustThreshold`, `kDefaultNoisyRegMaxXgaps` |
| Namespaces | snake_case | `pgphase_collect` |
| Bitmask flags | kCamelCase | `kCandCleanHetSnp`, `kLongcalldRepHetVar` |

### Header rules

- Use `#ifndef` / `#define` include guards (not `#pragma once`).
- Use `"file.h"` for project-local files.
- All code under `pgphase_collect` namespace.
- Include order: (1) corresponding header, (2) project headers, (3) third-party headers (`sdust.h`, `cgranges.h`), (4) system/htslib headers, (5) standard library.

### Implementation rules

- Implement functions in `.cpp`, not `.hpp` (except templates).
- Const correctness: const parameters, const references, const locals.
- Use RAII for resource management (`std::unique_ptr` with custom deleters for htslib handles).
- No magic numbers — use named `constexpr` constants.
- Prefer `static` file-scope functions over anonymous namespaces for small helpers.
- Match existing style in the file you're editing.

### Error handling

- Use `std::runtime_error` for fatal errors.
- Use `std::cerr` for warnings and verbose output (gated by `opts.verbose`).
- Make error messages actionable — include the file path or parameter that failed.

### Comments

- Explain "why", not "what".
- Use `///` or `/** */` for public API documentation.
- Use `//` for inline implementation notes.
- Do not comment obvious code.

### Testing

- Unit test binaries: `test_allele_identity`, `test_allele_context`,
  `test_allele_genotype`, `test_graph_sites`, `test_graph_bam_adapter`,
  `test_noise_filter`, `test_union_phase`. Standalone `.cpp` files in `src/`,
  hand-rolled `check()` assertions, no framework. Run all: `make unit-tests`.
- Catch2 predicates: `src/test_phase_predicates.cpp`, built against the vendored
  single header in `third_party/catch2/`. Run: `make predicate-tests`.
- Every change to the hybrid (`--bam`) must match or beat HiPhase on the same
  truth-scorable overlapping reads and input alignments, and correctly phase at least 80%
  of those reads. Require both total correctly phased reads and correctly
  phased reads in the dominant connected phase set to be at least HiPhase
  counts. Output-only rescue phase sets do not count as connected core
  coverage. Count unphased reads in the denominator and report correct,
  discordant, unphased and dominant connected correct counts. Perfect
  accuracy or zero newly discordant reads is not required; measure and report
  changes outside the targeted region and preserve unrelated regression checks.
- Compare HiPhase on the same overlapping reads and input alignments. When
  HiPhase phases reads better, investigate the differing variant calls,
  observations, filtering and phasing decisions, record the cause, and fix
  the deficit before accepting the change. A diagnosis alone does not meet
  the acceptance rule. Do not assume pgphase lacks the required information.
  Parental truth is evaluation evidence only and must not enter production
  phasing decisions.
- When adding a new `.o` dependency, update both the main `pgphase` target and any test targets that link the dependent object.

---

## Shell Script Conventions

- Start every script with `set -euo pipefail`.
- Use `#!/usr/bin/env bash`.
- Quote all variable expansions: `"${var}"`.
- Use `readonly` for constants.
- Provide `--help` / usage output.
- Use `die()` or similar for fatal errors with exit code 1.

---

## Git Conventions

- Commit messages: imperative mood, first line ≤72 chars, blank line before body.
- Add co-author: `Co-authored-by: Ona <no-reply@ona.com>`.
- Only stage files relevant to the current task.
- Do not commit files modified before the task began unless directly related.

---

## Design Decisions Log

All significant design decisions, evaluation results, and filtering rationale are recorded in `CHECKPOINT.md`. Update it when:

- Adding or removing a filtering stage.
- Changing how a pipeline classifies or handles candidates.
- Running evaluation benchmarks with new results.
- Discovering pitfalls or gotchas worth preserving.

---

## Code Review Checklist

Before submitting changes:

- [ ] Every changed line traces to the user's request.
- [ ] `make -j$(nproc)` compiles with zero errors.
- [ ] `make unit-tests` passes.
- [ ] No new warnings introduced.
- [ ] Makefile updated if new `.cpp` files added (both main target and test targets).
- [ ] `CHECKPOINT.md` updated if design decisions were made.
- [ ] No magic numbers — named constants used.
- [ ] Const correctness maintained.
- [ ] Thread safety considered (no shared mutable state without synchronization).
- [ ] htslib handles wrapped in RAII (`std::unique_ptr` with custom deleter).

## Documentation

`docs/IMPLEMENTATION.md` is the single living description of how the code works.
**A commit that changes behaviour updates it in the same commit.** Part I is
current behaviour and is authoritative; Part II is per-component detail carrying
its own revision date.

Measurements, findings and history do NOT go there -- they go in `evaluations/`
by date and in `CHECKPOINT.md`. Keeping them out is what lets the implementation
doc be trusted without checking its date. `docs/HANDOFF.md`,
`docs/hybrid_design.md` and `docs/phase_graph_implementation.md` are history and
say so in their own banners.
