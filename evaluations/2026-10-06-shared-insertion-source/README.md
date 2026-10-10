# Verified shared insertion and local deletion read connection

Close **chr20:4,874,129–4,884,130**, the remaining rescue-to-core connection
reported by the complementary-deletion evaluation. The owning 4–5 Mb replay
has **53 correct connected-core reads**, matching HiPhase's 53. Pgphase has
**54 correct / 1 discordant / 9 unphased** of 64 original truth-scorable primary
overlaps (84.375% overall correctness), versus HiPhase's **53 / 1 / 10**.
Connected-core correctness is 53/64, or 82.8125%, with abstentions included.
`audit_block.py` independently verifies the original/HiPhase alignment geometry
and scores the actual HiPhase BAM on the same primary read names.

The canonical catalog A>AC insertion at 4,878,943 preserves its verified BAM
source edit while retaining the graph representation and omitting the MSA
flag. The old promoter therefore rejects it. Recognizing its exact retained
MSA insertion and checking the original complete C-repeat sequence promotes
17 correct rescues into the existing core without changing their haplotypes.

Five more correct rescues end before that insertion. The verified upstream
single-base deletion at 4,866,154 distinguishes them. The old source-component
check rejects unrelated unadoptable rows far upstream in the native owning
context. Its local deletion certificate now validates two shared clean SNPs
and **every intervening anchor** in one gauge, including weak/quality-cut and
unique-provenance checks. It also requires a separately anchored, verified
shared insertion from a different source within 20 kb in the same core.
Reads need their own complete deletion/REF sequence, bounded repeat length
and edit-distance agreement, known flank qualities, at most 5% physical plus
mapping error, an agreeing existing rescue haplotype and no contrary clean SNP.
The existing repeat caller can handle a neighboring sequence difference that
causes strict deletion equivalence to abstain. Other unshared deletions retain
their whole-cohort guards. No new haplotype or variant-block join is inferred.

The owning replay promotes 23 rescue phase-set tags: 22 correct and one already
discordant. Every previous haplotype and ordinary core tag survives. The correct
MAT read `m84031_231217_062403_s3/109381842/ccs` disagrees with both physical
marker calls and retains its original rescue tag; the native regression checks
that abstention explicitly. Whole-owner counts remain **3,862 correct / 26
discordant / 238 unphased**. One/four-thread tags, all candidate bytes and complete
VCF records match. The closure is added to the certified panel, measured HiPhase
reference, native replay map, required insertion evidence, parental-flank checks
and measured expectations. Existing expectations retain their floors.

Reproduction from the repository root (Python audits require pysam):

```bash
LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu make gap-dev-check
LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu make gap-owner-check GAP=4.874
LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu make window-tests
bash evaluations/2026-10-06-shared-insertion-source/replay.sh
python evaluations/2026-10-06-shared-insertion-source/audit_owner.py
python evaluations/2026-10-06-shared-insertion-source/audit_block.py
python evaluations/2026-10-06-shared-insertion-source/audit_preservation.py
python evaluations/2026-10-06-shared-insertion-source/audit_panel.py
python evaluations/2026-10-06-shared-insertion-source/audit_validation.py
```

Final executable SHA256:
`c6bcf86a43b5301cd90960eeb8b525d14944494dfb8381f9b676c5641952e674`.
Build, all units, 13 new in-memory checks, 1,551 phasing-predicate assertions
in 47 cases, HiFi/ONT goldens and HiFi thread determinism pass. The unified
suite passes **18,138 assertions in five cases**, plus four cache-helper tests.
The new native fixture passes 757 assertions in about 1.1 seconds warm; the
starting binary fails three assertions, including 31 versus 53 core parity.

The final native-panel audit checks **250 output labels / 128 independent
requests / 128 comparison pairs** against immutable task-start output.
Twelve labels change read tags, preserving every previous correct/phased
assignment and every complete VCF record. The 53 Mb control retains all its
read tags and VCF records exactly.

Full chr20 has **230,965 correct / 6,620 discordant / 19,027 unphased** primary
reads, identical to the baseline. All 47 changed tags promote existing rescues
with the same haplotype and truth status: 22 correct plus one already discordant
in the target neighborhood, 23 correct near 13 Mb, and one correct near 48 Mb.
The latter two effects are reported without claiming new gap closures.
`full-preservation.json` lists every changed read, original alignment span and
before/after tags. Every previous correct/phased read, ordinary core tag,
all 64,483 complete VCF records, every byte of the 82,287-site candidate TSV
and every block extent survives. **N50 remains 955,496 bp**, largest block
3,012,193 bp, with 250 blocks.

The entire 625,922 bp HiPhase block at 4,570,122–5,196,043 still has a deficit
outside this repaired interval: pgphase has 2,511 total / 2,448 connected-core
correct reads, versus HiPhase's 2,515 / 2,515 on 2,622 original primary overlaps.
Whole-block parity is not claimed. The repaired 10,001 bp interval meets both
HiPhase counts and the 80% contract. `results.json`, `owner-results.json`,
`panel-audit.json`, `full-preservation.json` and `validation.json` record the
completed measurements and checks.
