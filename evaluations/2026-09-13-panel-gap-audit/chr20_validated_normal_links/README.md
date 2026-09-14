# Validated repeats with ordinary link scoring: rejected

Fourteen-case screen: ten joins and four splits, with no added assessed read errors. The 3.95 Mb target now joins correctly; the earlier wrong orientation is prevented by homopolymer classification and candidate validation. Two old joins (20.5 and 38.2 Mb) remain split.

Full 114-case evaluation: 89 joined, 23 split, two endpoint-unphased. Rejected because the 36.0 Mb window gains a wrong stitch of a five-read original block, also producing two variant hamming errors. Both 1.9 Mb endpoint repairs remain correct. This illustrates why endpoint join counts alone are insufficient.

This trial additionally relaxes repeat-edge scoring after candidate validation; its earlier strict-margin unit expectation is intentionally superseded only in the later experimental source, not production. Final production retains the established MSA update and repeat-link behavior and scopes only clean-candidate k-means updates.
