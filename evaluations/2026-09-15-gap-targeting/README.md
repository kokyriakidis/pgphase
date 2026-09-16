# Pinpointing the BAM region to attack, per unresolved gap

A gap does not have one informative region — it has a **linkage bottleneck**. A
gap is bridgeable only where some read crosses a position while carrying an
informative het site on *both* sides of that crossing. Scanning cut points across
the gap and counting such reads gives a linkage profile; the minimum of that
profile is the position where phasing actually breaks, and that position — not the
gap interval — is what a targeted recovery subprocess should be pointed at.

Reads do not need to span the whole gap for this to work: each cut only needs one
read reaching a site on either side of it, so a long gap is bridged as a chain of
overlapping reads. This is why the bottleneck is usually *interior* rather than at
a junction (87 of 164 gaps here), and why widening the recovery window around the
gap edges does not help.

```sh
python3 evaluations/2026-09-15-gap-targeting/locate_gap_evidence.py \
  --tiers /tmp/pgphase-private-snp-bridge/off/tiers.tsv \
  --candidates /tmp/pgphase-private-snp-bridge/off/candidates.tsv \
  --potential-vcf <het calls made using all reads> \
  --bam test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam \
  --output gap_targets.tsv
```

Two site sets are compared at each cut: the sites *we* called (`CLEAN_HET_SNP`,
`CLEAN_HET_INDEL` from the pipeline's own output) and the sites a caller finds
using every read regardless of mapping quality. The difference between the two
profiles is what separates "we failed to call the evidence" from "the evidence is
not there".

## Result on chr20: 164 unresolved gaps, only 34 are actionable

| cause | gaps | gap span | what the bottleneck lacks | remedy |
|---|---:|---:|---|---|
| `thin_linkage` | 130 | 11.96 Mb | no read carries a het site on both sides, even using all reads | none — information-limited |
| `sites_not_called` | 20 | 0.64 Mb | sites exist and a read chain exists; we did not call them | MSA discovery in the target window |
| `mapq_starved` | 12 | 2.59 Mb | reads exist but none pass the MAPQ floor at the cut | gated admission, then MSA |
| `linkage_present` | 2 | 0.02 Mb | nothing — reads and our own sites link across | solver or stitch defect |

**The `thin_linkage` verdict above was wrong, and a second pass overturned it.**
Both site sets in the first pass were BAM-derived — our own calls and a pileup
caller's het calls — and *both* omit repeat-context indels. The snarl catalog was
never consulted. It turns out to be dense at exactly those bottlenecks: a median
of **671 catalog sites within 25 kb**, of which ~96% are homozygous in this sample
(median 614 `ref_only` plus 31 `high_af`, unobserved sites negligible at 2%), but
with a median of **7 repeat-demoted het indels** (`REP_HET_INDEL`) and 4 low-AF
sites per window that are heterozygous and currently excluded from phasing.

Re-testing every bottleneck with those classes included
(`add_demoted_site_linkage.py`, `gap_targets_revised.tsv`):

| revised cause | gaps | gap span | meaning |
|---|---:|---:|---|
| `repeat_indels_would_bridge` | 88 | 6.45 Mb | admitting repeat-demoted het indels links the bottleneck |
| `low_af_sites_would_bridge` | 24 | 3.43 Mb | links once low-AF sites are added too |
| `no_linkage_from_any_site_class` | 52 | 5.32 Mb | no linkage at MAPQ >= 30 from any class |

So **112 of 164 unresolved gaps (9.89 Mb) have candidate linkage at their
bottleneck from graph sites the pipeline throws away** — median 17 linking reads
once repeat indels are included, in a median 50 kb target window. Of the 130 gaps
first written off as heterozygosity deserts, 99 fall in this group. The deserts
were an artifact of the site set, not a property of the sample.

Two honest qualifications. First, this is *candidate* linkage: those sites were
demoted because per-read genotypes at homopolymer and STR indels are unreliable,
so linkage computed from them may be phantom. They have to be MSA-verified inside
the target window before being trusted as anchors — which is exactly what the
existing `populate_gap_msa_cache` / `run_gap_msa_tier` path is for, and why this
result argues for that path rather than against it. Second, 12 of the 52
remaining gaps are the `mapq_starved` class, where no read passes the MAPQ floor
at the bottleneck at all, so for those this pass answers nothing — their
linkability is untested, not refuted. The genuinely unlinkable residue is at most
40 gaps, roughly 2.7 Mb. This pass also used the graph-only run's site table
rather than the hybrid candidates, so the `linkage_present` pair is not directly
comparable between passes.

## Validation

The scanner independently reproduces the diagnosis reached by hand in
`../2026-09-15-mapq-starved-gaps`: `chr20:25,834,662-25,883,079` comes back
`mapq_starved`, bottleneck at 25,870,298 (interior, not a junction), target
window 25,845,298-25,895,298, 20 of our sites against 98 potential ones nearby.
That diagnosis took a long manual session; the scan derives it in one pass and
does the same for the other 163 gaps.

## The per-gap subprocess this enables

`gap_targets.tsv` gives, for every unresolved gap, a `target_beg`/`target_end`
window plus a `cause` and `remedy`. A recovery subprocess can therefore be
dispatched per gap instead of running one policy chromosome-wide:

1. `sites_not_called` (20 gaps) — run gap MSA discovery restricted to the target
   window and admit verified hets as anchors. This is the case the existing
   `populate_gap_msa_cache` / `run_gap_msa_tier` machinery is built for, pointed
   at a 50 kb window rather than the whole gap.
2. `mapq_starved` (12 gaps) — admit sub-floor reads inside the target window only,
   and require MSA verification before any admitted site becomes an anchor. The
   whole-chromosome version of this corrupted 12 confidently mapped reads
   (`../2026-09-15-split-mapq-floor`); confining it to 50 kb windows where the
   scan shows the evidence is otherwise absent is what makes it defensible.
3. `linkage_present` (2 gaps) — do not touch evidence; these are solver or stitch
   defects and belong with the `split` investigation.
4. `repeat_indels_would_bridge` / `low_af_sites_would_bridge` (112 gaps, 9.89 Mb)
   — run MSA verification over the target window on the demoted het indels and
   low-AF sites already in the catalog, and admit only those the consensus
   resolves into two consistent haplotypes. This is the largest group and the
   one the graph catalog already has the sites for.
5. `no_linkage_from_any_site_class` (52 gaps, of which 12 are MAPQ-starved and so
   untested) — the only candidates for "information-limited", and even that
   should be re-checked with the low-MAPQ pass before anything is written off.

## One dependency worth removing

The "potential site" set here comes from an external caller's VCF, which is fine
for an audit but wrong for a production subprocess. Replacing it with a direct
pileup scan of the BAM over the target window — allele-balance test, no MAPQ
floor — makes the targeting self-contained, and that scan is work the MSA
discovery step would do anyway.

## Is the demoted-indel linkage real? Mostly not — 0.43 Mb of the 6.45 Mb

`verify_demoted_sites.py` genotypes every repeat-demoted het indel inside the 88
target windows straight from the alignment and scores how well its allele
partition follows the read-level truth (`diplinator` haplotypes, matched by read
name). A site whose alleles segregate with truth is real signal that screening
removed; one that splits reads at random is the noise the screen exists to catch.

| class | n | median segregation | informative (>= 0.90) | phantom (< 0.70) |
|---|---:|---:|---:|---:|
| `CLEAN_HET_INDEL` (control, trusted today) | 167 | 1.000 | 96% | 3% |
| noisy-candidate het (emitted today) | 138 | 0.728 | 25% | 48% |
| repeat-demoted het indel (excluded today) | 744 | 0.586 | 23% | 62% |

By type: demoted insertions are 33% informative, demoted deletions 18%; the
noisy-candidate class sits at 21-27% across DEL/INS/SNP.

**Recomputing each bottleneck with only the truth-validated sites admitted: 10 of
the 88 gaps retain linkage, spanning 0.43 Mb** — against 88 gaps and 6.45 Mb when
every demoted site is admitted. The 171 informative sites found across the 88
windows (about two per window) are real, but they are not positioned to span the
bottlenecks. So the apparent 6.45 Mb opportunity is carried mostly by phantom
sites, and even a *perfect* oracle gate — one admitting exactly the real sites —
would recover 0.43 Mb from this class. `apply_graph_noise_filter` is discarding
some genuine signal (23%), but it is right about the majority (62%).

### What this says about the verification gate

The noisy-candidate het class the pipeline emits today is only 25% informative
and 48% phantom — statistically indistinguishable from the demoted pile it would
be gating. That is *not* the same as saying MSA verification fails, and the
distinction matters: `NOISY_CAND_HET` in the output is the candidate **category**
(`collect_phase.cpp:72`), while `msa_verified` is a separate per-variant flag set
in `collect_phase_noisy.cpp:213` and consumed at `collect_phase.cpp:580,1188`.
The TSV does not expose it — only the `--gap-decision-audit` export writes
`msa_verified=` per variant (`collect_phase.cpp:303`). So the precision of MSA
verification itself is still unmeasured, and measuring it is the load-bearing
next step: the whole graph-blocks-plus-BAM-recovered-gap-blocks design rests on
MSA-verified sites being trustworthy as anchors, and right now the class that
carries them into phasing is three-quarters noise.

### Method corrections made while getting here

Two of my own measurement bugs, both caught by the same control — scoring
`CLEAN_HET_INDEL`, which the pipeline already treats as clean, and requiring it
to come out informative:

1. Genotyping indels at the VCF anchor position (`pileupread.indel != 0`) scored
   clean het indels at **2% scorable**. Repeat-context indels are placed
   arbitrarily within their repeat run by the aligner, so the allele has to be
   read from the net insert-minus-delete length over a window.
2. Deriving the expected length change from `len(ALT) - len(REF)` silently
   skipped every deletion, because this TSV writes deletions as the deleted bases
   in REF with `ALT` = `.` (`collect_output.cpp:110-145`). That inflated
   "unscorable" to 62% of sites and understated the real linkage.

After both fixes the control reads 96% informative with median segregation 1.000,
and 744 of 795 sites in the windows are scorable.
