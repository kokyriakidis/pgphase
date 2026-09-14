# Scoped k-means with partial-read refresh and original gap bounds

Rejected nine-case diagnostic. Preserving the original gap interval restores the existing 4.85 Mb join, and the 1.9 Mb endpoint repair remains. Rechecking partial MSA read observations does not prevent the wrong 3.95 Mb join: it still reverses the original right block (119 assessed reads).

The next diagnostic isolates the misclassified homopolymer insertion at 3,971,337. It has poor allele separation against both clean flanks but bypasses repeat validation because raw FASTA ASCII was compared to an nt4 code. This report does not establish that partial-read refresh is safe to deploy; that change is excluded from the subsequent implementation.
