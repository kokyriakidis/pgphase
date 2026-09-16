# Audit: the graph-only path (`collect-graph-variation`) on chr20

Whole-chr20 run against the snarl catalog and coordinate-indexed GAF, compared
site-by-site with the hybrid path's default arm and scored against the
diplinator read truth. Reproduce with:

```sh
./pgphase collect-graph-variation \
  --ref test_data/chm13v2.0.chr20.renamed.fa \
  --sites test_data/chr20.sites.striped.vcf.gz \
  --gaf test_data/HG002.chr20.annotated.coord.gaf.gz \
  -r 'CHM13#0#chr20' -t 8 \
  -o variants.tsv --phased-bam-out phased.bam \
  --phase-sites-out phase_sites.tsv --filtered-sites-out filtered_sites.tsv
```

## Where it stands

73,627 candidates retained, 914,975 catalog sites filtered (826,336 ref-only,
42,750 high-AF, 35,481 with no reads in chunk, 7,044 low-AF, 3,364 low depth).
Read-level accuracy against truth: 203,751 reads phased, 27,631 unphased, 371
phase sets, 1,744 discordant, read Hamming 0.856%, block N50 917,428 bp. The
hybrid default arm on the same chromosome phases 207,426 reads at 0.675% with
N50 1,064,616 in 177 phase sets, so the graph-only path is somewhat more
fragmented and roughly a quarter less accurate, which is the expected cost of
phasing catalog sites alone.

## Clean vs noisy classification: correct, and deliberately stricter

The noise filter does fire — 17,350 het indels are demoted to `REP_HET_INDEL`
(9,318 DEL, 8,032 INS) against 2,036 surviving as `CLEAN_HET_INDEL`, so 89.5% of
graph het indels are screened out. `batch_contig` comes from the FASTA-derived
synthetic header, so the reference fetch that guards `apply_graph_noise_filter`
resolves; the filter is not silently skipped in this configuration (it would be,
without warning, if that fetch ever failed — the `else` branch just frees and
continues).

At the 58,394 sites shared with the hybrid path, categories agree on 57,033
(97.7%). Every systematic disagreement is in the conservative direction: 397
sites the graph calls `REP_HET_INDEL` that the BAM classifier calls
`CLEAN_HET_INDEL`, plus 426 it demotes where the BAM path files them as
`NOISY_CAND_HET/HOM`. Nothing enters graph k-means that the BAM path would have
screened as repeat noise.

## Defect: equal-length multi-base substitutions are typed SNP and escape screening

`build_graph_chunk` derives the variant type from allele lengths alone
(`src/graph_bam_adapter.cpp:944-950`): shorter ref means insertion, longer ref
means deletion, **equal lengths mean `VariantType::Snp`**. `VariantType` has no
MNP member (`Snp`/`Insertion`/`Deletion`, mapped to CIGAR ops), so every
multi-base substitution becomes a SNP carrying `key.ref_len = ref.size()`.

Consequences on chr20:

- **781 sites** are equal-length multi-base substitutions typed `SNP`. All 781
  are classified `CLEAN_HET_SNP` — the highest-confidence anchor class, scoring 2
  in `pick_pivot`/`var_score` — and all 781 carry a `PHASE_SET`, so every one of
  them voted in k-means.
- They can never be demoted. `apply_graph_noise_filter` reconsiders only
  candidates whose category is `CleanHetIndel`, so a multi-base substitution in
  homopolymer or STR context is structurally unreachable by repeat screening.
- The context is frequently repetitive: of the first 400, 115 sit in a
  homopolymer or dinucleotide-repeat window — `GTG>CTC` inside `GTGTGTGTGTGTG`
  (chr20:145,051), `TAT>GAG` inside `TATATATATAGAG` (288,018), `CT>TC` inside
  `CCCTCTCTCTCCC` (863,876). The BAM path demotes the equivalent sites: at
  47,301,789 an `AT>TA` swap is `REP_HET_INDEL` there and `CLEAN_HET_SNP` here.
- Allele lengths run 2-12 bp (571 dinucleotide, 63 trinucleotide, the rest
  longer) and include at least one 367 bp equal-length allele pair at
  chr20:764,862 typed `Snp` with `ref_len = 367`, which also makes the key's span
  disagree with its type.

The BAM path is not right either — it types these as `Deletion`, which is also a
mislabel — but it at least routes them through repeat screening. Neither path
represents MNPs.

`mnp_typed_as_snp.tsv` lists all 781 with depth, AF, category and phase set.

## MSA validation: absent from this path

No MSA runs in the graph-only path. `abpoa` is reached only through
`src/align.cpp`, whose verification consumer is `src/collect_phase_noisy.cpp`
(the only place `msa_verified` is set true); `src/graph_collect.cpp` does not
include it and never sets `msa_verified`. Consistently, the graph-only output
contains no `NOISY_CAND_HET`/`NOISY_CAND_HOM` calls at all, against 4,003 and
1,964 in the hybrid arm.

So the graph-only path phases clean catalog sites and has no noisy-site rescue:
there is no MSA confirmation, and no mechanism to recover phase in regions where
the catalog offers only repeat-context indels. That is a capability boundary
rather than a bug, but it is the reason its unphased count (27,631) is 3.3x the
hybrid arm's.

## Two further observations

- **No homozygous calls are emitted.** 42,750 catalog sites are filtered with
  reason `high_af` and the output contains zero `CLEAN_HOM` records, while the
  hybrid arm emits 34,722. Harmless for phasing, an asymmetry for anything
  consuming the VCF as a call set.
- **The multi-allelic branch of `classify_graph_candidates` is unreachable as
  configured.** It classifies `n_uniq_alles > 2` sites as `CleanHetIndel` and, in
  doing so, skips the AF-centering anchor gate that the surrounding comment
  introduces specifically to stop paralog-collapse sites from voting. All 73,634
  retained sites are biallelic by then (decomposition sets `n_uniq_alles = 2`),
  so nothing currently takes that branch — dead logic that would bypass the
  paralog guard if decomposition order ever changed.
- The evaluator reports `switch_errors = 0` and `flip_errors = 0` for this path
  because the emitted BAM is unaligned; only the read-level concordance and
  Hamming figures above are meaningful for graph-only.

## Recommended fix for the defect

Smallest change that closes the screening hole without touching `VariantType`:
extend `apply_graph_noise_filter` to also examine `CleanHetSnp` candidates whose
catalog alleles are multi-base, and move those in low-complexity context out of
the anchor mask. The low-complexity intervals are already computed there
(`find_low_complexity_intervals` / `pos_in_low_complexity`), and
`classify_graph_candidates` already has the right destination mask for
"emit but do not vote" in `kCandNonAnchorHet` / `kLongcalldLowAfVar` under
`--emit-nonanchor-hets`. Note `is_noisy_site` derives indel length from allele
size difference, which is zero for these sites, so the demotion must key on
low-complexity overlap rather than on that helper.
