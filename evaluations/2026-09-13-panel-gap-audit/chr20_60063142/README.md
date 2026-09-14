# chr20 60.1 Mb: unresolved multi-allelic insertion evidence

The remaining wrong join spans native anchors 60090439–60114911. HiPhase's
phased heterozygous evidence inside this interval includes VCF 60093416,
C→CTTTTTT / CTTTTTTT, a genotype with two alternate insertion lengths.
The local MSA likewise produces six-T and seven-T insertions at event position
60093417, but represents them as independent candidates. Exact sequence
observation leaves each candidate with alternate support and no reference
support (13 and 23 ALT observations respectively in the diagnostic run), and
both collapse to homozygous. The phasing core loses the 6-versus-7 distinction.

This representation issue is diagnosed, not fixed. A solution must retain two
alternate alleles in read observations and preserve their identity in output;
calling the shorter insertion genomic reference would be incorrect. The current
local recovery still reverses 39 original-block reads. Truth and competitor
alleles are used only for diagnosis, not runtime phasing.
