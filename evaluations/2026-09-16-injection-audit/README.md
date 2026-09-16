# Auditing one window's BAM-site injection, in four stages

`audit_injection.py` takes what `collect-bam-variation` itself calls in the
interval as the reference set -- same reads, same reference, no graph channel --
and asks four questions about it, each returning a finding rather than a number
to be interpreted:

1. **PRESENT** -- does the hybrid chunk hold every het site the BAM channel calls?
2. **FIELDS** -- are depth, allele counts, allele fraction, type and category identical?
3. **ADMITTED/USED** -- is a row of that locus emitted as a phased het with a phase set?
4. **INFORMATIVE** -- do its alleles segregate with the read-level truth?

Stage 4 carries the clean-het control this project requires. Run on
`chr20:48,176,830-48,229,446` with `--retry-unphased-with-bam`:

```
window 48176830-48229446   BAM channel: 79 sites, 9 het
  1 PRESENT   9 of 9 het sites in the hybrid chunk (0 absent)
  2 FIELDS    5 het sites differ between channels
  3/4 USED    8 used as a phased het, 0 emitted homozygous, 1 unscorable
      control: 3 clean het sites, informative: YES (0.986, 1.000, 1.000)
      informative among used: 5
```

## What passes

Injection is not the problem. **Every het site the BAM channel calls is present
in the hybrid chunk, and 8 of 9 loci are emitted as phased hets with a phase
set.** No locus is emitted homozygous any more. That was the state this audit
was built to check, and it holds.

## What the audit found, and what it retracted

**The BAM channel emits both nested forms of a repeat deletion at one position as
independent contradictory hets.** At `48,177,780` it writes
`GAGAAAGAA>G 1|0` (29/45) *and* `GAGAAAGAAAGAAAGAAAGAA>G 0|1` (46/28); at
`48,225,786`, `CA>C 0|1` (35/16) and `CAAAAAA>C 1|0` (33/18). One locus cannot
carry both. This is exactly what `split_nested_msa_deletions` exists to prevent,
and it is off in the BAM channel, so the channel used as the reference set is
itself the one in error here. With the splitter on, the hybrid decomposes each
into the common deletion, homozygous because both haplotypes carry it, plus a
phased residual het a few bases away.

**Retracted: "the split destroyed the information".** The residual at
`48,177,789` scores 0.500 against truth, against 0.946 for the unsplit 20 bp
allele, which looked like the split flattening an informative site. It is not:
the reads there carry two deletion lengths, -8 (27 reads) and -20 (21 reads), and
a net-length test for a 12 bp residual has a window that necessarily contains the
common 8 bp deletion -- so a -8 read shows -8 there, is closer to -12 than to 0,
and is called alt as well (66 alt against 8 ref). **A split residual is
unscorable by this method by construction**, and the audit now says so and falls
back to the unsplit allele instead of reporting chance. At `48,225,787` the same
split is scorable and reads 1.000, so the mechanism is sound where it can be
measured.

## Open findings in this window

| locus | finding |
|---|---|
| `48,202,057` | used as a phased het but **segregates 0.507** over 69 reads -- and 0.509/0.507/0.515 tested as a 2, 3 or 4 bp deletion, so it is not a representation artifact. Not a split residual. A site with no haplotype information is being used to phase. |
| `48,177,726` | the BAM channel calls `NOISY_CAND_HET` (34/40, AF 0.541); the hybrid reclassifies it `REP_HET_INDEL` (35/39) and emits no record, so the locus links nothing. |
| `48,225,787` | 0.733 on the common 1 bp allele; its residual at `48,225,788` reads 1.000, so the locus is fine but the common part is not informative on its own. |
| `48,202,057`, `48,229,227`, `48,177,726` | allele counts differ from the BAM channel by 1-3 reads (`DP` 69 vs 68, 26/34 vs 23/37, 34/40 vs 35/39). Small, but the same reads and reference should give the same counts. |

`48,202,057` is the one to fix next: a site at chance segregation carrying a phase
set is the failure mode the whole gap effort has been chasing, and it is not
explained by representation.
