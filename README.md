# pgphase

pgphase phases long reads (PacBio HiFi or ONT) using their alignments to a
pangenome graph. Each read gets a haplotype (`HP`) and phase set (`PS`), and the
heterozygous variants are written to a phased VCF.

Variant sites come from the graph's snarls, so pgphase can phase on variation
that a linear reference represents badly. Given the reads' linear alignments as
well (`--bam`), it adds the sample's own heterozygous sites, those not in the
graph, and phases both kinds of site together ("union gap phasing"). This is the
recommended mode.

## Commands

| Command | What it does |
|---|---|
| `collect-graph-variation` | Phase reads from their graph alignments (GAF), optionally with their linear alignments (`--bam`) |
| `collect-bam-variation` | Phase reads from linear alignments alone (BAM/CRAM); a C++ port of [longcallD](https://github.com/yangao07/longcallD)'s phasing |
| `build-snarl-catalog` | Build the per-chromosome site catalog that `collect-graph-variation` reads |

## Installation

Requirements:

- a C++17 compiler (GCC 8+ or Clang 7+)
- [htslib](https://github.com/samtools/htslib) with headers, zlib, pthreads
- [vg](https://github.com/vgteam/vg), only for `build-snarl-catalog`
- [pggaf](https://github.com/kokyriakidis/pggaf), to prepare the GAF input
- [samtools](https://github.com/samtools/samtools), only for `collect-bam-variation --refine-aln`

```bash
git clone --recursive https://github.com/kokyriakidis/pgphase.git
cd pgphase
make third-party-libs   # WFA2-lib and abPOA
make -j"$(nproc)"
./pgphase --version
```

`make portable-bundle` writes a self-contained binary with its shared libraries
to `dist/pgphase-<os>-<arch>/`.

## Quick start

pgphase needs four inputs: the reference FASTA, a site catalog built from the
graph, the reads' graph alignments as a coordinate-indexed GAF, and (recommended)
the same reads' linear alignments as an indexed BAM.

### 1. Build the site catalog (once per graph)

```bash
for chr in chr{1..22} chrX; do
    pgphase build-snarl-catalog --ref-sample CHM13 --contig "$chr" -t 16 \
        -o "$chr.sites.vcf.gz" graph.gbz
done
```

For each chromosome, `build-snarl-catalog` extracts the chromosome's connected
component (`vg chunk -C`), computes its snarls (`vg snarls`), and runs
`vg deconstruct -a` with them on the reference path `<ref-sample>#0#<contig>`,
writing a bgzipped, tabix-indexed sites-only VCF. On the HPRC v2.1 CHM13 graph,
chr20 takes 28 minutes on 10 cores with a 25 GB peak. The chromosomes are
independent jobs.

Do not pass a whole-genome snarls file with `--contig`: `vg deconstruct`
enumerates the haplotypes through every snarl it is given before checking which
reference path it belongs to, so each chromosome would pay for the whole genome.

### 2. Index the GAF by reference coordinate

```bash
pggaf annotate-gaf --gaf reads.gaf --gbz graph.gbz --ref-sample CHM13 \
    --out-gaf reads.annotated.gaf --out-sets reads.pgs
pggaf index-gaf --in reads.annotated.gaf --out reads.coord.gaf.gz
```

### 3. Phase

```bash
pgphase collect-graph-variation \
    --ref chm13.fa \
    --sites chr20.sites.vcf.gz \
    --gaf reads.coord.gaf.gz \
    --bam reads.bam \
    -r 'CHM13#0#chr20' \
    -o candidates.tsv \
    --phased-vcf-out phased.vcf \
    --phased-bam-out phased.bam \
    -t 16
```

Leave out `--bam` to phase from the graph alignments alone.

## Outputs

| Option | Content |
|---|---|
| `-o FILE` | Every candidate site with its allele counts, category and phase |
| `--phased-vcf-out FILE` | Heterozygous calls, phased (`GT:DP:AD:VAF:GQ:PS`) |
| `--phased-bam-out FILE` | One unaligned record per read seen, with `HP` and `PS` tags when phased |

Read tags can be copied onto an aligned BAM by read name, for example for
haplotype-aware variant calling.

## `collect-graph-variation` options

| Option | Description | Default |
|---|---|---|
| `--ref FILE` | Reference FASTA (indexed) | required |
| `--sites FILE` | Site catalog from `build-snarl-catalog` | required |
| `--gaf FILE` | Coordinate-indexed GAF from pggaf | required |
| `--bam FILE` | Indexed BAM of the same reads (union gap phasing) | off |
| `-t INT` | Worker threads | 1 |
| `-q INT` | Minimum mapping quality | 5 (1 with `--bam`) |
| `-D INT` | Minimum total depth | 5 |
| `--min-alt-depth INT` | Minimum alt depth | 2 |
| `--min-af FLOAT` / `--max-af FLOAT` | Allele fraction band for a het | 0.20 / 0.80 |
| `--chunk-size INT` | Region chunk size | 500000 (1000000 with `--bam`) |
| `-r STR` / `--region-file FILE` | Restrict to regions | whole genome |
| `--hifi` / `--ont` | Read technology | `--hifi` |

`pgphase collect-graph-variation --help` lists everything, including
diagnostic outputs.

## `collect-bam-variation`

Linear-only phasing, as in longcallD, for when no graph alignment is available:

```bash
pgphase collect-bam-variation --ref ref.fa --bam reads.bam --hifi \
    -o candidates.tsv --phased-vcf-out phased.vcf -b phased.bam -t 16
```

Use `--ont` for ONT reads; `--help` lists all options.

## Testing

```bash
make unit-tests   # unit tests
make check        # end-to-end golden tests on the bundled data in test_data/
```

`make upstream-parity-tests LONGCALLD_ROOT=/path/to/longcallD` compiles
longcallD's original C and checks the ported phasing functions against it.

## Licence

MIT; see `LICENSE`. Third-party code and its licences are listed in
`THIRD_PARTY_NOTICES.md`.
