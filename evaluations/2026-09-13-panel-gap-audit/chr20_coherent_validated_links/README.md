# Coherent per-site repeat validation: rejected

Fourteen-case screen: thirteen joined targets and one split. Allowing consistent allele-group margins and a pure clean SNP with Q30 support on each allele restores the 20.5 and 38.2 Mb joins. However, 23.4 Mb now joins with an incorrect original-right-block orientation (32 assessed reads reversed). This is not accepted.

The full source experiment is archived as remaining_coherent_validated_experiment.patch. It includes phase-scoped updates in all rounds, homopolymer classification, original gap bounds, independently validated MSA anchors, joint repeat-genotype updates and ordinary link scoring for eligible repeats. Production excludes this combined experiment.
