"""The supplement, written once; render.py turns it into paper/supplement.tex and build/paper/supplement.html.
Same markup as main_text.py; its tables come from tables.STABLES."""
from tables import STABLES

TITLE = "Supplementary material for pgphase: phasing long reads at a pangenome's known alleles and the sample's own heterozygous variants"
TABLES = STABLES

BODY = [
('sec', 'Per-chromosome results'),
('p',
 'Table S1 gives, per chromosome and tool, the phaseable reads on the arms by class, the discordant reads whose '
 'parental label is confident, the switch and flip errors of a record-level comparison and the phased truth '
 'heterozygotes. The MHC lies inside chr6.'),
('table', 'perchrom'),
('sec', 'Switch errors by haplotype sequence'),
('p', 'Table S2 classes every arm switch of Table S1 by the haplotype sequences around it (main text Section 2.8).'),
('table', 'census'),
('sec', 'Running time and memory'),
('table', 'resources'),
('sec', 'Commands'),
('p',
 'Each tool was run on the reads of one chromosome, CHM13 contig names without the graph\'s `CHM13#0#` prefix. '
 '`REF` is the chromosome\'s CHM13 sequence, `BAM` the surjected reads, `VCF` the DeepVariant calls, `T` the '
 'threads.'),
('list',
 ['DeepVariant 1.10.0: `run_deepvariant --model_type=PACBIO --ref=REF --reads=BAM --regions=CHROM '
  '--output_vcf=VCF --num_shards=T`',
  'HiPhase 1.6.0: `hiphase --reference REF --bam BAM --vcf VCF --threads T --output-bam phased.bam '
  '--output-vcf phased.vcf.gz`',
  'LongPhase 2.0.2: `longphase phase -s VCF -b BAM -r REF -t T --pb -o phased`, then `longphase haplotag -s '
  'phased.vcf.gz -b BAM -r REF -t T -o phased`',
  'WhatsHap 2.8: `whatshap phase --reference REF --ignore-read-groups --indels [--distrust-genotypes] -o '
  'phased.vcf.gz VCF BAM`, then `whatshap haplotag --reference REF -o phased.bam phased.vcf.gz BAM`',
  'longcallD v0.0.11: `longcallD call --hifi -o phased.vcf -b phased.bam -t T REF BAM CHROM`',
  'pgphase: `pgphase collect-graph-variation --ref REF --bam BAM --sites CATALOG --gaf GAF -r CHROM --threads T '
  '-o candidates.tsv --phased-vcf-out phased.vcf --phased-bam-out phased.bam`',
  'Catalog and GAF index: `pgphase build-snarl-catalog` on the GBZ graph; `vg convert -G` of the reads\' GAM, then '
  '`pgphase index-gaf --graph GBZ --reference-sample CHM13`',
  'Truth labels: `minimap2` of each read to each parental haplotype of HG002 v1.1, then `hiphap` on the two '
  'alignments']),
]
