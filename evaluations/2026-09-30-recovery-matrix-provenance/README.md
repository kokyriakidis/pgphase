# Preserve each recovery matrix dump (2026-09-30)

With `--phase-matrix-dump`, targeted BAM groups all used
`RegionChunk::chunk_id == -1`. Every group's initial solve, MSA retry, and
selected source therefore wrote the same `PREFIX.chunk-1.*.tsv` paths.
Later groups silently replaced earlier evidence. The whole-chunk BAM read solve
also reused the graph chunk ID, replacing the graph's `flags*.tsv` matrices;
the full-block validation solve could replace the last targeted BAM matrix.

The dump prefix now names the graph chunk and recovery window. Initial, MSA,
focused, remainder, and selected-source dumps have separate names. Validation
and whole-chunk BAM solves have their own prefixes. This changes diagnostics
only: the saved 58-59 Mb and 57-58 Mb phased VCFs and candidate TSVs are
byte-identical before and after the first naming fix.

The corrected 58-59 Mb run writes three source matrices, covering
58.198-58.501, 58.510-58.630, and 58.654-58.885 Mb. The old single
`chunk-1.recovery-source.tsv` retained only the last group and did not
contain the 58.366-58.385 Mb gap. In its actual source, the left clean SNP
at 58,366,458 and right deletion at 58,385,703 are both represented, but
no read has callable alleles at both. Intervening graph rows provide no
strong diploid path between them. Thus the old "empty source" interpretation
was an artifact, while the conclusion to abstain from a join remains.

The 15-16 Mb replay writes three source matrices; the 15.351-15.367 Mb
seam is in the middle one. Its left deletion has only three ALT observations,
and paired reads to the right insertion call only the left REF. The 57.854
Mb source has no heterozygous interior BAM site between the two boundary
rows. Graph indels at 35.498 Mb have mixed or one-sided links to its
boundary deletion and SNP. None of these audits justify a new phase-set
connection.

At 21.159-21.172 Mb, 18 graph-matrix reads call the left SNP and the
right clean SNP, but only one calls right ALT; the others call right REF
(10 with left REF and seven with left ALT). The right deletion has a strong
local link to the right SNP (26 agreeing versus two disagreeing), yet no
read calls both the left SNP and that deletion. Thus the right block's
internal path does not supply the missing cross-gap orientation.
