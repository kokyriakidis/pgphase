# Remaining chr18 33.3 Mb boundary: singleton support

The remaining boundary after site-observation recovery is
chr18:33284672–33308109, a 23,437 bp link in the upstream flank of the originally
audited gap (33319321–33429051). The original target was already connected;
previous reports counted blocks across the entire solve window
chr18:33269321–33479051. HiPhase PS33179887 covers this upstream boundary too.

Both endpoint clean SNPs and their allele observations are already in our
phasing input. HiPhase has no heterozygous variant strictly inside this interval.
There is exactly one primary, nonsupplementary BAM read spanning both endpoints:

- `m84031_231217_062403_s3/253169506/ccs`
- Reference alignment [33281712,33309001), MAPQ60.
- Position 33284122: G, base quality 40; 33284672: C, quality 35;
  33308109: C, quality 40.
- Our allele matrix has one 0/1 observation linking 33284672 to 33308109.
  The same read also supplies the 33284122→33308109 link; those two site pairs
  are not two independent reads.
- HiPhase tags that read HP1/PS33179887. This is compatible with a singleton
  bridge; HiPhase's internal edge-selection trace was not instrumented.

The pgphase default `--min-block-link-reads 2` rejects this evidence. Retrieving
more copies of the existing sites or lowering the MSA margin does not create a
second molecule. The existing option `--min-block-link-reads 1` is sufficient:
the entire regional output becomes one block, with all 339 previously tagged
reads retained, zero read-truth discordance and zero switchflips. The two old
blocks both map uniformly to PS33284122 without internal haplotype changes.

No source-code change or implicit lowering of the default is needed for this
regional fix. The solved BAM is at
`/tmp/pgphase-singleton-link-validation/chr18_33319321/auto.bam`.
Use the exact command in `chr18_33319321/command.txt` to reproduce it.
The MSA margin 24, allele-based linking window 8, read margin 2, and normal
stitching settings are unchanged.

## Broader check

`run.py` reruns the same six fixed regions with min-block-link-reads 1 on the
recovery arm. `summarize.py` compares against the preceding site-observation
recovery at min-block-link-reads 2, and records read/block transformations.
It reports failures as well as successes; it does not assume a global singleton
policy is safe. This is targeted read-truth evaluation, not a chromosome-wide
variant-truth or NGC50 assessment. Production defaults remain unchanged.

| Region | Blocks with link minimum 2 → 1 | Discordant reads before → after |
|---|---:|---:|
| chr12 46.8 Mb | 4 → 2 | 0 → 0 |
| chr12 63.8 Mb | 1 → 1 | 0 → 0 |
| chr18 46.2 Mb | 1 → 1 | 1 → 1 |
| chr18 33.3 Mb | 2 → 1 | 0 → 0 |
| chr20 61.8 Mb | 1 → 1 | 0 → 0 |
| chr20 17.6 Mb | 2 → 1 | 0 → 0 |

All 2,386 tagged reads from the minimum-2 recovery outputs are retained.
Each previous block undergoes one uniform PS/parity transformation; no HP
label changes occur in these cases. Read-switchflip counts also remain unchanged.

### Clarification of the two remaining chr12 blocks

The original chr12 target is already one block, matching HiPhase. The remaining
boundary at 46842224–46874196 is in the downstream test-window flank; HiPhase,
WhatsHap, WhatsHap-opt and LongPhase all break there too. Thus the regional
4→2 count does not indicate a remaining failure to bridge the audited target.
Exact endpoint phase sets and the zero-spanning-read check are recorded in
`../2026-09-13-chr12-remaining-boundary/`.
