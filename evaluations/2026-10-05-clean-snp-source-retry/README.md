# Clean SNP retry with a complete independent BAM suffix

Target chr20:57,764,235–57,785,224 is a genuine open gap. Both tools have
128 correct / 7 discordant / 15 unphased among the same 150 original overlaps.
Before the fix, pgphase has only 65 correct reads in its largest core block;
the owning replay now places all 128 into one core, matching HiPhase (85.33%).

The exposed right boundary is a complementary MSA deletion family, followed
by a clean BAM SNP at 57,785,772. The deletion-family helper has one qualifying
Q20 molecule and correctly rejects its insufficient likelihood. A separate
primary MAPQ-60 molecule calls the left SNP REF and the clean right SNP ALT,
both at BQ40 (`physical-pairs.json`). The current gauges agree, with a
0.000202 two-base/two-mapping error bound and log odds −8.50704. The retry
nevertheless required observations of both left haplotypes before reaching
its ordinary physical likelihood and path checks.
Native HiPhase phases both clean SNPs under PS 57655293 even though its
18 bp boundary deletion call is missing (`hiphase-endpoints.json`). The
advantage is the available clean-SNP connection, not extra repeat information.

The retry now permits one-haplotype physical SNP pairs only when the right
block is a complete independent BAM-only run. Every anchor must have one
adoptable claim in the same source, a consistent current gauge, and no weak
or quality cut within the local run. Mixed source blocks retain the old
retry requirement. Unanimity and the 0.001 physical wrong-parity bound remain.
A complete graph path or BAM run certifies the left flank. The existing
deferred bridge retains every anchor and rechecks final block gauges after
rescue; the union changes only core labels. No read admission, genotype,
confidence threshold or production use of parental truth is introduced.

The owner audit preserves all 4,054 scored outcomes (3,860 correct / 194
incorrect), all rescue tags, variant keys, genotype alleles and nonphase
fields. It closes exactly this gap. Disjoint flanks agree in parental
orientation (100/0 and 113/7). See `owner-audit.json`.

The committed panel adds the gap with `spans=1`, a 0.85 core fraction floor,
its two deletion descriptions and both clean SNP witnesses. The owning
regression uses 57,000,001–58,000,000 and checks the existing owner read floors
and opposite SNP genotypes in the same phase set. Existing tests and required
sites are retained.

Run evaluation scripts from the repository root using bench-phasers Python
and `LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu`:

```bash
python evaluations/2026-10-05-clean-snp-source-retry/evaluate_gap.py \
  --before test_data/tmp_gap_fix63/full-deferred-source-runs \
  --after test_data/tmp_gap_fix64/full-source-retry \
  --output evaluations/2026-10-05-clean-snp-source-retry/full-audit.json
python evaluations/2026-10-05-clean-snp-source-retry/compare_hiphase.py
```

Full chr20 closes exactly this gap, reopens none, and preserves every read
correctness status and rescue tag. Only 129 core read phase labels and five
variant phase labels change. All 64,188 keys, genotype alleles and nonphase
fields are retained. The 237,334 scored reads remain 230,566 correct / 6,768
incorrect. Read phase sets fall 644→643, VCF blocks 324→323, and N50 remains
856,770; 99/112 tracked windows span. See `full-audit.json`,
`full-parity.json` and `hiphase-comparison.json`.

Final checks pass: build without compiler warnings, all unit tests, 47 phase
predicate cases / 1,536 assertions, deterministic HiFi/ONT validation gates,
and the complete native suite of 86 cases / 12,264 assertions over all 112
panel windows. The new regression fails three assertions on the preceding
binary and passes 747 after the fix. Frozen final production binary SHA256:
`f31a9166e0c8e871b28c3174ee95286e8bdc8cbd42c184dafc0e6eb7ce0019e6`.
See `validation.txt`, `native-validation.json` and `manifest.json`.
