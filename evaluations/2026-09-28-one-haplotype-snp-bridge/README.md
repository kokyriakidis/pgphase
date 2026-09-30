# One-haplotype physical SNP bridge on chr20

At chr20:32,215,055–32,233,534, both flanks are clean graph SNPs, but
recovery leaves separate phase sets. The original BAM has two distinct primary
reads with MAPQ 38/40 and Q40 calls at both sites. Both read pairs are C/A:
ALT at the left SNP and REF at the right SNP. Their combined quality-weighted
likelihood exceeds the stitcher's 0.001 wrong-parity bound. The previous rule
rejected any bridge with more than one read unless the reads covered both
haplotypes. That discarded this valid low-coverage diploid parity observation.
The existing graph-path and reversing-vote checks still guard whole-block joins.

The owning 32–33 Mb chunk joins the boundary as 1|0 and 0|1 in one phase set;
local truth counts remain 1,311 correct and 208 discordant among 1,519 scored
reads. The regression runs the owning chunk with 50 kb of context and checks
this exact boundary orientation, span and local truth score. The existing
62.6 Mb unsafe whole-block control remains a separate test.

Full chr20, with the same graph, GAF and surjected BAM, changes tracked
HiPhase-correct gap closures from 10/23 to 11/23 and read phase sets from 717
to 716. It leaves 62,154 VCF variant keys and all 236,847 truth-scored reads
unchanged: 229,090 correct, 7,757 discordant (96.7249%). The newly closed
gap is the only change to the tracked target closure list.

Reproduction: run `collect-graph-variation` on
`CHM13#0#chr20:32000001-33000000` and on `CHM13#0#chr20:1-66210255`
with `test_data/chm13v2.0.chr20.renamed.fa`,
`test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam`,
`test_data/chr20.sites.striped.vcf.gz`, and
`test_data/HG002.chr20.annotated.coord.gaf.gz`. Score the two outputs with
`/tmp/score_pgphase_weakcut.py`, or run the committed window regression.

Validation: `make -j8`, `make unit-tests`, `make check`, and `make window-tests`
pass. The full window suite reports 1,741 assertions in 34 cases. Comparing
the two full-chromosome BAM outputs after orienting each phase set to parental
truth finds the same 236,847 read names and zero changed correctness outcomes.
