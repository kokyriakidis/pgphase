# Correction: chr12's remaining boundary is also a competitor break

The original audited gap is chr12:46725702–46838918. With the explicitly tested
minimum-one-read setting, pgphase phases both endpoints in PS46709288, matching
HiPhase. The two reported pgphase blocks refer to the larger solve window
chr12:46675702–46888918, which includes a downstream competitor break.

At the remaining boundary (46842224–46874196), HiPhase changes from PS46709288
to PS46874196. WhatsHap, WhatsHap-opt and LongPhase change from PS46810372 to
PS46874196. None bridges this boundary. There are zero primary, nonsupplementary
input BAM reads covering both endpoint SNPs. This does not establish the absence
of every possible future source of evidence, but it rules out treating this as
an existing competitor bridge that pgphase has missed.

The earlier regional block-count report was accurate numerically but lacked
this target/flank distinction. No further merge or threshold relaxation is
justified by the comparison. No phasing code was changed in this audit.

Run `python3 evaluations/2026-09-13-chr12-remaining-boundary/audit.py` to recreate
`endpoint_phase_sets.tsv` and `spanning_reads.json` from the frozen competitor
VCFs, current regional output, and input BAM.
