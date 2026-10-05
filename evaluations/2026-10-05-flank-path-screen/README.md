# Further flank-path search: no accepted gap closure

All production trials were removed. The accepted starting and restored binary
is SHA256 `e635c76de38cc00ab0554529a278ca297b321a3a728634a40bc7652fb8b78c71`.
The baseline outputs are `test_data/tmp_gap_fix58/full-projected-final`.
No panel row, expectation, implementation description or production behavior
changes in this investigation.

Three separate full-chromosome experiments produce unchanged read tags,
genotype alleles, phased keys, phase-block gauges and parental truth statuses:

1. Check left as well as right graph flanks for the existing independently
   supported internal switch repair.
2. Also certify a failed graph edge when unanimous Q30/MAPQ30 primary BAM
   pairs confirm its existing orientation. Keep both haplotypes, two molecules
   per haplotype, count p<=0.01, quality-weighted wrong-parity <=0.001 and the
   complete remaining graph-path checks.
3. Accept a complete, consistently oriented, cut-free original BAM source run
   as the right path in clean-SNP stitching, symmetrically with the left path.

Each leaves 327 VCF blocks, 230,563 correct / 6,767 discordant scored reads,
and 95/108 tracked spans. None connects any newly uncovered gap or reopens a
tracked span. The adjacent audit JSON files contain exact comparisons. Native
outputs, frozen trial binaries and build logs are in `test_data/tmp_gap_fix59/`.
Full runs use the standard chr20 reference, annotated BAM, striped catalog and
coordinate GAF, with eight threads.

## Rejected actual join at 32,330,932–32,333,068

The imported left BAM source ends in a noisy five-base deletion at 32,330,932.
Its preceding clean G>A SNP at 32,327,334 pairs with the right C>G SNP at
32,333,071 on 17 primary Q30/MAPQ30 molecules. Both physical parities are
unanimously cross relative to the existing block genotypes (1/16 by left
haplotype), with log odds 139.506. The source-run predicate accepts the left
block and the right graph path passes.

An isolated trial lets a validated left BAM source nominate its nearest clean
SNP within 10 kb before an indel boundary. A fresh 32,300,001–32,400,000 replay
closes exactly this previously uncovered 2,136-base gap and retains all variant
keys, genotype alleles and uniform old block gauges. It is nevertheless
**unsafe**: independently oriented parental scores fall from 52 correct / 16
wrong to 43 correct / 25 wrong. Nine previously correct reads become wrong,
with no read improving. Reject the join and remove the fallback.

A complete original-source identity and unanimous physical boundary parity
are therefore insufficient to guarantee the biological orientation of a
repeat-rich source block. Preserve the independent parental audit; do not
accept this gap merely because the source-path predicate passes.

## Other rejected leads

The failing graph edge at 63,210,717–63,210,874 has 26 unanimous same-parity
primary Q30 pairs, on both haplotypes. Certifying it exposes another blocker
at 63,733,935–63,734,075. The latter right SNP sits immediately after a
100-base insertion and all 36 paired primary calls are reference C at that
right position, while the left SNP is heterozygous (19 T / 17 C). This is an
alignment/graph representation discrepancy; it does not justify declaring
the left SNP homozygous. Diagnostic omissions of either graph row produce no
owning-chunk output change and were removed.

Four actual gaps have a compound SNP/indel left boundary that the simple
SNP/indel boundary screens do not cover: 40,633,644–40,636,354,
45,416,750–45,439,920, 49,800,887–49,845,376 and 57,584,041–57,602,614.
The 45 and 49 Mb pairs have no primary MAPQ30 spanners; the 57 Mb pair has
one, with Q22 at its right insertion. The 40 Mb pair has 30 spanners but no
canonical right REF G call: aligned bases are T or the locus is deleted.
The left SNP and insertion classes both pair with T, so these calls cannot
supply a diploid boundary orientation. Full read-level observations are in
`compound-evidence.json`.

The wider screen checks the nearest three clean or noisy SNPs within 10 kb
of each side of every actual gap <=50 kb. None nominates a unanimous Q30 pair
with at least two molecules on each haplotype and two-sided count p<=0.01.
The corrected gap search advances the maximum covered endpoint across all
preceding blocks, rather than treating adjacent sorted-start blocks as gaps.

The Python screens run from the repository root using pysam, the above
accepted baseline VCF and the standard annotated chr20 BAM. They write their
reproduced JSON to `test_data/tmp_gap_fix59/`. The nomination screens do not
use parental truth. Parental labels are used only by the output audits.

Restoration verification: source matches its start-of-task snapshot byte for
byte, the rebuilt binary matches the accepted SHA256 exactly, and unit tests
pass. The restoration build introduces no warnings. Saved logs and the binary
manifest are adjacent to this document.
