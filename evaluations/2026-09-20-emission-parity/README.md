# Parity bug: the VCF writer re-applied detection thresholds at emission

Found by asking which gate rejects the records upstream emits and we do not.

## The measurement that located it

Of the records upstream's own chr20 output carries and ours does not, 910 sit at
positions we never emit -- and we HOLD a candidate at 748 of them, all with a
phase set already assigned:

| category of the held candidate | count | phase state |
|---|---:|---|
| `NOISY_CAND_HET` | 674 | all phased |
| `CLEAN_HET_SNP` | 64 | all phased |
| `NOISY_CAND_HOM` | 9 | all phased |
| `CLEAN_HET_INDEL` | 1 | phased |

So the sites exist, they are phased, and something drops them at emission.
Instrumenting the four gates in `project_vcf_candidates` over whole chr20, on
phased candidates only:

| gate | phased candidates rejected |
|---|---:|
| `passes_vcf_depth_gates` | **1,316** (1,271 of them `NoisyCandHet`) |
| `is_alt_genotype` | 1,149 (763 noisy het, 358 clean het SNP) |
| `is_germline_output_category` | 0 |
| `lcd_make_variants_region_pass` | 0 |

## The bug

`passes_vcf_depth_gates` required `total_cov >= min_depth && alt_cov >=
min_alt_depth` at emission -- the DETECTION thresholds, re-applied to counts
that the MSA and the co-located merge have since rewritten.

longcallD does not do this. `make_variants` (`collect_var.c:1465-1562`) filters
on exactly two things: the candidate's category
(`var_i_to_cate & target_var_cate`) and whether the position falls inside the
active region. Its only reference to depth in that whole function is the
assignment `var->vars[i].DP = cand_vars[cand_i].total_cov`. min_dp and
min_alt_dp are applied where they belong, at candidate detection
(`collect_var.c:910-911`).

The thresholds themselves are already identical on both sides -- `min_dp 5`,
`min_alt_dp 2` -- so this was purely a duplicated gate.

## Fixed, and what it changes

| alignment arm, whole chr20 | records | identical to upstream | upstream-only | ours-only |
|---|---:|---:|---:|---:|
| before | 116,170 | 114,111 | 4,159 | 2,059 |
| **after** | **116,896** | **114,407** | **3,863** | 2,489 |
| upstream | 118,270 | -- | -- | -- |

296 of the recovered records are ones upstream emits; the remaining 430 fall in
the separate "ours-only" bucket, which is a different parity question (our
candidate set differs from upstream's at those loci) and is not addressed here.

Read placement is untouched, as it must be for a writer change: 216,962 reads
tagged, 349 read blocks, 1,516 discordant, 0.699% hamming -- identical before
and after. Positions putting two ALTs on one haplotype stay at 0. The graph arm
is unaffected (62,458 records either way) because it emits through
`graph_chunks_to_candidate_table`, not this writer.

Unit 3/3, predicate 151/151, window 125/125.

## Still open, with its measurement

`is_alt_genotype` rejects 1,149 phased candidates, and that one is NOT a
divergence to delete: upstream emits no `0|0` records either (chr20 GT
distribution is 41,705 `1|0`, 41,308 `0|1`, 35,257 `1|1`, nothing else), and its
emission loop reaches the same outcome by appending an ALT only for a
non-reference haplotype allele (`collect_var.c:1556-1558`).

The real divergence at those sites is the consensus VALUE. At `26,591,362` and
four neighbours, upstream emits `0|1:8:5,3` while our candidate carries the
same DP 8 and the same 5/3 split, `CLEAN_HET_SNP`, phased -- but
`hap_to_cons_alle = [0,0]`, both haplotypes called reference, so no ALT is
written. Both tools compute that consensus the same way: our
`update_var_hap_to_cons_alle` (`collect_phase.cpp:223-236`) is a faithful port
of upstream's (`assign_hap.c:244-268`), same argmax, same prefer-reference
tie-break, same ONT guard; and our complement inference
(`collect_phase.cpp:263-268`) ports `assign_hap.c:139-142` exactly. So the
divergence is in the read-to-haplotype assignment that fills the profile, in
low-coverage regions -- upstream also puts those five sites in ONE phase set
(26,591,362) where we split them across two (26,549,599 and 26,591,520).

That is the next parity step, and it is a solver-level question, not a gate.

## The `is_alt_genotype` class is NOT a bug on our side

Chased to the read level at `26,591,362`, one of the 358 `CLEAN_HET_SNP`
rejections, where upstream emits `0|1:8:5,3` and we emit nothing.

**Both tools are working from the same 8 sampled reads.** The site has 445
covering reads (330 reference, 115 alternate) yet both report DP 8 -- the
noisy-region read sampling -- and both label exactly the same 8 read names.

**Seven of the eight labels agree. The one that differs decides the record:**

| read | base | ours | upstream | truth |
|---|---|---|---|---|
| `215485304` | T | HP2 | HP2 | PAT |
| `36700407` | A | HP1 | HP1 | MAT |
| `64292608` | A | HP2 | HP2 | PAT |
| `85791316` | T | HP2 | HP2 | PAT |
| `143853999` | A | HP2 | HP2 | PAT |
| `54662206` | T | HP2 | HP2 | PAT |
| `73531961` | A | HP1 | HP1 | MAT |
| **`77401560`** | **A** | **HP2** | **HP1** | **PAT** |

Our labelling is 8/8 consistent with parental truth -- every HP2 read is
paternal, both HP1 reads are maternal. Upstream places `77401560`, a paternal
read, on HP1 beside the two maternal ones.

That one read is the whole difference in the output. Our hap2 profile becomes
`[3,3]`, a tie, which the prefer-reference argmax resolves to reference, so
`hap_to_cons_alle = [0,0]` and no ALT is written. Upstream's hap2 is `[2,3]`,
alternate wins, and it emits `0|1`. **The record upstream emits rests on its
misplaced read.**

Note also what the site looks like with full coverage: 115 of 445 reads carry
the alternate, and among our own 8 reads the paternal haplotype carries BOTH
bases. That is a collapsed-duplication signature, not a clean heterozygote.

### How much of the class is real

Sampling 30 of the dropped positions, 21 with a single-base upstream REF, 13
with at least 5 truth-attributed alternate reads:

| the alternate allele's parental purity | positions |
|---|---:|
| >= 0.90 -- a real heterozygote we suppress | **8** |
| mixed parents -- an artifact, suppression correct | **5** |
| fewer than 5 alt reads, not evaluable | 8 |

So roughly 60% of the class are real sites we lose and 40% are sites upstream
should not be calling.

### Why this is not fixed here

Three candidate fixes, all rejected on the measurement:

- **Change the tie-break toward ALT on het-classified sites.** Produces
  upstream's output but diverges from its rule, and would emit the 40%
  artifact fraction as heterozygous calls.
- **Change the read assignment to match upstream.** Ours is the more accurate
  labelling at this locus (8/8 against 7/8), so this would trade truth for
  parity.
- **Emit the tie as an unphased het.** Upstream emits no unphased records at
  all on chr20 (41,705 `1|0`, 41,308 `0|1`, 35,257 `1|1` and nothing else), so
  it does not close the gap either.

The class is therefore recorded as a known, quantified difference and not a
defect. The remaining parity buckets that ARE candidate defects: 162 positions
where we hold no candidate at all, and 172 loci where we hold the position but
a different allele (148 indels, 24 SNPs).

## Parameter audit: every threshold matches upstream

Checked because a clustered discovery difference (162 positions in 12 windows
of 10 kb, 146 of them SNPs, largest clusters 50 at 30.80 Mb and 89 across
26.68-26.71 Mb) would be explained by any one of these differing.

| parameter | upstream | ours |
|---|---|---|
| min mapping quality | 30 (`LONGCALLD_MIN_CAND_MQ`) | 30 (`kDefaultMinMapq`) |
| min candidate depth | 5 | 5 |
| min alternate depth | 2 | 2 |
| min allele fraction | 0.20 | 0.20 |
| max allele fraction | 0.80 | 0.80 |
| max noisy region length | 50,000 | 50,000 |
| max noisy region coverage | 1,000 | 1,000 |
| noisy region max X/gaps per window | 5 | 5 |
| slide window, HiFi / ONT | 100 / 25 | 100 / 25 |
| max noisy fraction per read | 0.5 | 0.5 |
| long end clip / clip flank | 30 / 100 | 30 / 100 (`kLongClipLength`, `kClipFlank`) |
| noisy region merge distance | 500 | 500 |
| noisy region flank length | 10 | 10 |
| sample reads above region size | 10,000 | 10,000 |
| min reads supporting a noisy region | `min_alt_dp` (`collect_var.c:603`) | `min_alt_depth` (`collect_var.cpp:1001`) |

All seventeen match, including the one upstream leaves as a commented-out
option and derives from `min_alt_dp` instead -- we derive it the same way.

So no remaining parity difference is a threshold. In the 30.80-30.81 Mb window
we hold 264 candidates and emit all 264 where upstream emits 315, from reads
that are mostly MAPQ 1-19 (21 of 28 over the first 2 kb, only 7 at or above the
shared floor of 30). At that depth the difference is which reads survive
digar-level filtering and how the noisy region is cut, not a parameter.
