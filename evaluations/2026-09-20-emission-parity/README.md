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
