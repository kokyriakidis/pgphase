# Independent evidence for gap bridges

Two additional competitor-resolved chr20 gaps now join using existing cached
graph/BAM/MSA observations and the normal phasing/stitching core.

## 37.98 Mb: correlated observations vetoed the clean signal

Read `m84031_231217_034919_s2/158272567/ccs` observes a verified insertion at
VCF anchor 37,984,529, a verified deletion at 38,005,399, its neighboring SNP
at 38,005,401, and a separate SNP at 38,008,026. The neighboring deletion/SNP
observations disagree with the separate SNP in the right block's orientation.
The old component unanimity check discarded all evidence from that read.

For bridge anchors, SNPs within or immediately bordering an MSA-verified,
link-supported indel no longer count as independent clean SNPs. For each
read/component, independent clean SNPs are considered first; verified indels
are used when those SNP observations are absent. Conflicting independent SNPs
still reject the anchor. A confident indel cannot supply the confidence flag
for a weak SNP that displaced it. Exact BAM-validated biallelic insertions can
anchor singleton bridges under the same quality and conflict checks used for
deletions.

This joins the original blocks with opposite orientation, correctly against
parental truth. Local evaluation remains 458 concordant / 459 evaluated reads,
one discordant read and one flip. No truth labels or competitor calls enter
the phasing decision.

## 0.86 Mb: distinct insertion alleles were excluded as a group

The bridge uses `TC` versus `TCTC` at native insertion position 882,278.
Its singleton was excluded by the general multiallelic-insertion gate even
though its observed allele has an exact, high-quality BAM match. Singleton
eligibility now also accepts two alternate insertion sequences whose lengths
differ by at least two bases and neither consists of one repeated nucleotide.
Existing MSA verification, linkage support, sequence-equivalent BAM matching,
Q30 checks and the opposing-read veto still apply.

This joins the original blocks without flipping their relative orientation.
All 435 locally evaluated reads remain concordant, with no switch/flip errors.
The `C`/`CC` case at 23.48 Mb remains excluded by this singleton rule.

## Validation and remaining cases

`panel.tsv` compares the committed deletion-bridge baseline with
`/tmp/pgphase-gap-trials/runs/separated_insertion_v1`. Both changed regions
have parental read evaluations. In the other nine regions the ordered read
names, mapping starts, flags, HP and PS assignments match the previous run
exactly; their existing truth evaluations are reused and explicitly labeled.
The 37.98 Mb result is identical to its independently evaluated
`anchor_tiers_v1` result. The panel has 3 joined targets, 7 splits and 1
unphased endpoint, with no local read-error increases.

Full chr20 output is `/tmp/pgphase-anchor-final-check`. It joins 125 of 280
initial gaps (previously 123), with no previously recovered edge lost. The new
native gaps are 864,134–890,260 and 37,982,955–38,005,401. Cached gap solving
took 192.868 seconds; evidence-cache loading took 1.823 seconds. These timings
are a single run, not a statistically established speedup.

Full-chromosome parental read evaluation changes from 2,731 discordant /
189,001 evaluated to 2,734 / 189,041. All previously evaluated reads retain
their truth-concordance status; none are lost. The 40 additional reads comprise
37 concordant and 3 discordant assignments. Those three have truth HapQ 0;
they remain counted in the raw metrics. Switch counts change 353 → 355 and
flip counts 413 → 414 with the added reads. This is additional uncertain read
tagging, not a detected orientation regression in the existing blocks.
`common_read_regression.json` and `chr20.read_summary.json` preserve the checks.
No new shared-VCF NGC50 or variant Hamming metric is claimed.

The remaining cases are not claimed fixed. At 1.08 Mb the only observed
right-block SNP on the bridge read has Q10. At 19.37 Mb a candidate bridge
read has an exact four-base `TTCC` insertion, but an inserted base and the
left boundary have Q22 and fail the existing Q30 validator. That threshold
was not relaxed in this change.

Build and unit tests pass. New tests exercise correlated indel-boundary
observations, conflicting independent SNPs, unverified boundaries, low-quality
SNPs, singleton sequence-equivalent insertions, mononucleotide alternatives,
and low-quality insertion observations. The earlier shifted-insertion,
deletion and orientation tests remain enabled.
