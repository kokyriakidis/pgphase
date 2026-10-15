# Graph pipeline test data — chr20:25,000,001–25,500,000

A 500 kb slice of CHM13 chr20 with HG002 HiFi reads, aligned both to the HPRC
pangenome (GAF) and to CHM13 (BAM). `make check` runs union gap phasing on it
and compares the output with `test_data/expected/union_chr20_expected.*`.

## Running

```bash
./pgphase collect-graph-variation \
  --ref   test_data/graph_chr20/ref.fa.gz \
  --sites test_data/graph_chr20/chr20_25M.sites.vcf.gz \
  --gaf   test_data/graph_chr20/HG002_chr20_25M.coord.gaf.gz \
  --bam   test_data/graph_chr20/HG002_chr20_25M.bam \
  -o output.tsv --phased-vcf-out phased.vcf --phased-bam-out phased.bam
```

Expected: 576 candidate sites; 2,006 of the 2,032 reads in the phased BAM carry
an `HP` tag. Without `--bam` the run phases from the graph alignments alone.

`collect-bam-variation` runs on the BAM alone:

```bash
./pgphase collect-bam-variation \
  --ref test_data/graph_chr20/ref.fa.gz \
  --bam test_data/graph_chr20/HG002_chr20_25M.bam \
  --hifi -o output.tsv
```

## Files

| File | Size | Description |
|------|------|-------------|
| `ref.fa.gz` | 261 KB | bgzipped reference FASTA; chr20 N-padded to genome coords |
| `ref.fa.gz.fai` | 23 B | FASTA index |
| `ref.fa.gz.gzi` | 6.2 KB | bgzip index |
| `HG002_chr20_25M.bam` | 12 MB | HiFi reads aligned to chr20 (contig = `chr20`, LN=25,500,000) |
| `HG002_chr20_25M.bam.bai` | 17 KB | BAM index |
| `chr20_25M.sites.vcf.gz` | 213 KB | Graph snarl sites VCF (tabix-indexed) |
| `chr20_25M.sites.vcf.gz.tbi` | 488 B | Tabix index |
| `HG002_chr20_25M.coord.gaf.gz` | 5.6 MB | Coordinate-sorted GAF, HG002 HiFi (tabix-indexed) |
| `HG002_chr20_25M.coord.gaf.gz.tbi` | 2.1 KB | Tabix index |

## Notes

**Reference FASTA**: The slice covers positions 25,000,001–25,500,000 on chr20.
The contig is N-padded from position 1 so that FASTA coordinates match genome
coordinates. This is required because pgphase tiles the reference into chunks
and uses absolute positions for tabix queries against the sites VCF and GAF.

**Sites VCF**: Subset from the full CHM13 chr20 snarl catalog with
`tabix ... CHM13#0#chr20:25000001-25500000`. The pipeline resolves
`CHM13#0#chr20` → `chr20` via suffix matching.

**GAF**: Coordinate-sorted and bgzipped, indexed with `tabix -p bed`.
The format has reference contig/start/end in columns 1–3.

## Regenerating

```bash
# BAM (subset from full chr20 BAM, rename contig, fix header length)
(samtools view -H full_chr20.bam \
     | awk '!/^@SQ/ || /SN:CHM13#0#chr20/' \
     | sed 's/CHM13#0#chr20/chr20/g; s/LN:66210255/LN:25500000/g'
 samtools view full_chr20.bam "CHM13#0#chr20:25000001-25500000" \
     | sed 's/\tCHM13#0#chr20\t/\tchr20\t/g') \
    | samtools sort -O bam -o HG002_chr20_25M.bam
samtools index HG002_chr20_25M.bam

# Reference
samtools faidx ref_full.fa "CHM13#0#chr20:25000001-25500000" \
    | grep -v "^>" | tr -d '\n' > seq.txt
python3 -c "
seq = open('seq.txt').read().strip()
print('>chr20')
padded = 'N'*25000000 + seq
[print(padded[i:i+80]) for i in range(0, len(padded), 80)]
" | bgzip > ref.fa.gz
samtools faidx ref.fa.gz

# Sites VCF
tabix sites_full.vcf.gz "CHM13#0#chr20:25000001-25500000" \
    | bgzip > chr20_25M.sites.vcf.gz
tabix -p vcf chr20_25M.sites.vcf.gz

# GAF
tabix reads_full.coord.gaf.gz "chr20:25000001-25500000" \
    | bgzip > HG002_chr20_25M.coord.gaf.gz
tabix -p bed HG002_chr20_25M.coord.gaf.gz
```
