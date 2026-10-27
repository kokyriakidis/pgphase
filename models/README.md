# Homopolymer models

plegma, pgphase's noisy-region MSA, decides each homopolymer run's length with a model of how reads report run
lengths: P(m | n), a read's length m for a true length n, per base, read strand and flanking letters, in shasta2's
csv (`base,strand,left,right,n,m,count,probability,logProbabilityDb`). Both are built into pgphase (make writes
them into src/msa_models.cpp); `--msa-model FILE` takes another, `--msa-model none` runs without one.

- `homopolymer-hifi-hg002.csv` (the default for HiFi reads): learnt from HG002 PacBio HiFi reads (Revio, CHM13
  alignments) on chr20, at reference runs of 2 to 50 letters with no HG002 truth variant within 20 letters,
  complementary strands pooled, lengths seen too rarely taken from the nearest well-seen shorter one
  (pgphase-dev scripts/dev/hifi_model.py). On held-out chr12 the length it finds from 10 reads is exact for 95% of
  runs of 20 to 24 letters.
- `homopolymer-ont-r10.4.1-sup5.0.0.csv` (the default with --ont): shasta2's model for ONT R10.4.1 reads basecalled
  with dorado sup 5.0.0 (github.com/paoloshasta/shasta2, conf/HomopolymerModel-ONT-R10.4.1-sup5.0.0-2026-09-28.csv,
  MIT licence; see THIRD_PARTY_NOTICES.md). Not yet evaluated in pgphase against a truth.
