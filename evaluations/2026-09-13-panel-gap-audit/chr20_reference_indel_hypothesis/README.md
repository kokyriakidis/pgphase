# Rejected reference-indel hypothesis trial

Six cases retain baseline target outcomes, but the 34.1 Mb window introduces
an original-block reversal (2 to 63 discordant reads). The 12.74 Mb deletion
recovers 26 reference / 48 alternate observations but still fails to join.
This trial restores the reference hypothesis when both MSA consensuses carry
the same alternate indel, reassesses provisional homozygous indels, handles
partial read coverage, and retains the original HP search bounds. It excludes
the earlier phase-scoped-read and generalized genotype-forcing changes, so
those changes are not required for the 34.1 Mb failure. The cumulative source
snapshot is remaining_reference_indel_experiment.patch; it is diagnostic,
not accepted production code. Local windows only; baseline multi_final.
