# Graph pipeline test data

Small pangenome graph test case for `collect-graph-variation`, covering the
MICB (chr6) and KIR3DL1 (chr19) loci from the HPRC pangenome.  Based on the
haplotype sampling test data from
[gbz-base](https://github.com/jltsiren/gbz-base/tree/main/test-data).

## Running the graph pipeline

```bash
./pgphase collect-graph-variation \
  --ref test_data/graph/ref.fa.gz \
  --sites test_data/graph/micb-kir3dl1.sites.vcf.gz \
  --gaf test_data/graph/micb-kir3dl1_HG003.coord.gaf.gz \
  -o output.tsv
```

Expected output: **292 candidate variant sites** (53 on chr6, 239 on chr19),
the same with `--sites test_data/graph/micb-kir3dl1.catalog`.

To restrict to a single locus:

```bash
./pgphase collect-graph-variation \
  --ref test_data/graph/ref.fa.gz \
  --sites test_data/graph/micb-kir3dl1.sites.vcf.gz \
  --gaf test_data/graph/micb-kir3dl1_HG003.coord.gaf.gz \
  --region "CHM13#0#chr6:31350873-31363898" \
  -o output.tsv
```

## Files

### Source files (checked in)

| File | Description |
|------|-------------|
| `micb-kir3dl1.gbz` | GBZ pangenome graph (CHM13 + GRCh38 refs, 46 samples) |
| `micb-kir3dl1_HG003.gaf` | GAF alignments for sample HG003 (12,439 reads) |

### Derived files (checked in, can be regenerated)

| File | Description |
|------|-------------|
| `micb-kir3dl1.catalog/` | Site catalogs from `build-snarl-catalog`: `chr6.pgsites` (202 sites), `chr19.pgsites` (715), and the contigs' snarls (`snarls/`), from which `make check` rebuilds them without vg |
| `micb-kir3dl1.sites.vcf.gz` | The same 917 sites as a sites-only VCF from `vg deconstruct -a` (ALTs in vg's order) |
| `micb-kir3dl1.sites.vcf.gz.tbi` | Tabix index for sites VCF |
| `micb-kir3dl1_HG003.coord.gaf.gz` | The GAF indexed by reference position, from `index-gaf` (`make check` rebuilds it) |
| `micb-kir3dl1_HG003.coord.gaf.gz.tbi` | Its tabix index |
| `ref.fa.gz` | Reference FASTA (bgzipped, N-padded to genome coordinates) |
| `ref.fa.gz.fai` | FASTA index |
| `ref.fa.gz.gzi` | bgzip index |

### Notes on the reference FASTA

The GBZ subgraph has offset paths (e.g., CHM13 chr6 starts at genome position
31,350,872).  The reference FASTA is N-padded so that FASTA coordinates match
genome coordinates.  This is required because pgphase tiles the FASTA into
chunks and uses those coordinates for its region queries.

## Regenerating derived files

```bash
# Site catalogs (requires vg and pgcat: make pgcat)
./pgphase build-snarl-catalog --graph test_data/graph/micb-kir3dl1.gbz \
  --reference-sample CHM13 --out test_data/graph/micb-kir3dl1.catalog

# Reference-indexed GAF (requires pgcat)
./pgphase index-gaf --graph test_data/graph/micb-kir3dl1.gbz \
  --gaf test_data/graph/micb-kir3dl1_HG003.gaf --reference-sample CHM13 \
  --out test_data/graph/micb-kir3dl1_HG003.coord.gaf.gz

# Sites VCF (requires vg, bgzip and tabix)
vg deconstruct -a -P 'CHM13#0#' test_data/graph/micb-kir3dl1.gbz | cut -f1-8 \
  | bgzip > test_data/graph/micb-kir3dl1.sites.vcf.gz
tabix -p vcf test_data/graph/micb-kir3dl1.sites.vcf.gz

# Reference FASTA (requires vg + samtools)
# 1. Extract CHM13 paths from the GBZ
vg paths -x test_data/graph/micb-kir3dl1.gbz -S CHM13 -F > /tmp/chm13_paths.fa

# 2. N-pad each contig so FASTA coordinates match genome coordinates
python3 -c "
import re, sys
seqs = {}
with open('/tmp/chm13_paths.fa') as f:
    name = None
    for line in f:
        if line.startswith('>'):
            m = re.match(r'>(\S+)\[(\d+)\]', line)
            name, offset = m.group(1), int(m.group(2))
            seqs[name] = (offset, [])
        else:
            seqs[name][1].append(line.strip())
for name, (offset, parts) in seqs.items():
    seq = ''.join(parts)
    print(f'>{name}')
    padded = 'N' * offset + seq
    for i in range(0, len(padded), 80):
        print(padded[i:i+80])
" > /tmp/ref_padded.fa

# 3. Compress and index
bgzip -c /tmp/ref_padded.fa > test_data/graph/ref.fa.gz
samtools faidx test_data/graph/ref.fa.gz
```
