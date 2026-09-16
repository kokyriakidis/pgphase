# Phase the gap in windows and stitch, using the machinery that already exists

## How the pipeline chunks today

`split_region` (`collect_pipeline.cpp:239`) cuts a region into **non-overlapping
fixed windows**, `chunk_beg += chunk_size`, default `kDefaultChunkSize = 500000`
and settable with `--chunk-size`. Adjacent chunks are stitched by
`flip_chunk_hap` (`collect_phase.cpp`), which votes on the reads present in both
chunks (`up_ovlp_read_i` / `down_ovlp_read_i`) and feeds
`select_stitch_orientation`.

The property that matters here is its refusal: `if (n_cur_ovlp_reads <= 0) return
false`. **No read spanning the boundary means no merge.** That is the safety
guarantee a bespoke gap link does not have.

## Windows alone change nothing

Swept `--chunk-size` over `chr20:48,126,830-48,279,446` with the shipped
configuration (recovery on):

| `--chunk-size` | blocks | phased sites INSIDE the gap | spans gap |
|---|---:|---:|---|
| 500 kb | 2 | **0** | no |
| 50 kb | 2 | **0** | no |
| 20 kb | 2 | **0** | no |
| 10 kb | 3 | **0** | no |

At 10 kb it is worse -- the left two-site block fragments into two single-site
blocks. The window size is irrelevant because nothing inside the gap is admitted
to any window's solve: with `--recover-gaps`, every candidate with
`!graph_site` has its category zeroed before phasing
(`evaluations/2026-09-16-bam-sites-initial-solve/`).

## With the BAM's own sites admitted, the window size decides correctness

Same region, dropping `--recover-gaps` (which keeps the BAM categories) and
adding `--keep-noisy-kmeans` (which orients the noisy-candidate class):

| arm | blocks | in-gap phased sites | read accuracy |
|---|---:|---:|---:|
| 500 kb, recovery on (shipped) | 2 | 0 | 100.00% |
| 500 kb, sites admitted | 2 | 79 | 100.00% |
| **20 kb, sites admitted** | 4 | **77** | **100.00%** |

The 20 kb arm produces `48,147,227-48,183,976` (10 sites, holding the left
flank's own site 48,162,480, so it is anchored there by construction) and
`48,224,939-48,279,445` (50 sites, likewise anchored to the right flank). The
500 kb arm produces a single 82 kb block, `48,147,227-48,229,226`.

That difference is not cosmetic. Between the gap's two site clusters,
`48,183,976 -> 48,225,786` is **41.8 kb with zero reads covering both sites**
(the neighbouring spacings carry 41 and 51; the longest read in the region is
29.9 kb). Scoring each block's sites against read truth:

| arm | block | sites left of the hole | sites right of it | crosses it |
|---|---|---|---|---|
| 20 kb | `48147227` | 7 (6 agreeing, 1 weak) | 0 | no |
| 20 kb | `48224939` | 0 | 48 agreeing | no |
| 500 kb | `48147227` | 7 (6 agreeing) | **2, opposite orientation** | **yes, switched** |

The wide window asserts a relative phase across a hole no read supports, and gets
it wrong here. The narrow windows stop at the hole and emit two blocks, each
anchored to one flank, everything they claim being correct. Read-level accuracy
cannot distinguish the two -- both read 100.00% -- because no read spans the hole;
only the site-level truth check sees it.

## This retracts the closure claimed for this gap

`evaluations/2026-09-16-gap-lab/` reported this gap CLOSED with +107 concordant
reads and no flips. The frame that closed it spans the same 41.8 kb hole, and its
orientation across it is a coin flip: the gap arm run with 5 kb of read context
and the same arm run with 30 kb produce **opposite** genotypes at 48,225,786 and
48,229,226, and the 30 kb one links the right flank backwards (159 reads
discordant, caught only by the truth validation). The 5 kb run was right by luck.

The honest outcome for this gap is two blocks, not one: the left cluster anchored
to the left flank, the right cluster to the right flank, and 41.8 kb between them
that no read-based method can phase.

## What this implies for the fix

The mechanism to build is not a new link. It is the user's: phase the gap in
windows with the BAM's own sites admitted, and let `flip_chunk_hap` stitch them,
because its `n_cur_ovlp_reads <= 0` refusal is exactly the guard that makes an
unsupported join impossible.

Two things block adopting it as it stands, both already measured:

- `--recover-gaps` and the BAM's in-gap sites are mutually exclusive today, since
  recovery zeroes non-graph categories before phasing. The scoped change is to
  admit them inside gap intervals only.
- `--keep-noisy-kmeans` is chromosome-wide and documented at the call site as
  having phased ~8k extra reads at ~65% error, so it also needs the same interval
  scoping.

Reproduce:

```sh
./pgphase collect-hybrid-variation \
  --ref test_data/chm13v2.0.chr20.renamed.fa \
  --bam test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam \
  --graph-sites test_data/chr20.sites.striped.vcf.gz \
  --gaf test_data/HG002.chr20.annotated.coord.gaf.gz \
  -r 'CHM13#0#chr20:48126830-48279446' -t 8 -q 1 \
  --chunk-size 20000 --keep-noisy-kmeans \
  --link-by-alleles --block-link-window 8 --min-read-margin 2 \
  -o cand.tsv --phased-vcf-out native.vcf -b phased.bam
```
