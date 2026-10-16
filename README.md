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
| `build-snarl-catalog` | Build the site catalogs that `collect-graph-variation` reads, from the GBZ graph (runs `pgcat build`) |
| `index-gaf` | Index the reads' graph alignments (GAF) by reference position, for `--gaf` (runs `pgcat index-gaf`) |

## Installation

Requirements:

- a C++17 compiler (GCC 8+ or Clang 7+)
- [htslib](https://github.com/samtools/htslib) with headers, zstd, zlib, pthreads
- [Rust](https://rustup.rs) (cargo), to build `pgcat`, which `build-snarl-catalog`
  and `index-gaf` run
- [vg](https://github.com/vgteam/vg), only for `build-snarl-catalog`
- [samtools](https://github.com/samtools/samtools), only for `collect-bam-variation --refine-aln`

```bash
git clone --recursive https://github.com/kokyriakidis/pgphase.git
cd pgphase
make third-party-libs   # WFA2-lib and abPOA
make -j"$(nproc)"
make pgcat              # catalogs and GAF indexing, placed beside pgphase
                        # (pgphase finds it there, then on PATH; PGCAT=PATH overrides)
./pgphase --version
```

`make portable-bundle` writes a self-contained binary with its shared libraries
to `dist/pgphase-<os>-<arch>/`.

## Quick start

pgphase needs four inputs: the reference FASTA, a site catalog built from the
graph, the reads' graph alignments as a coordinate-indexed GAF, and (recommended)
the same reads' linear alignments as an indexed BAM.

### 1. Build the site catalogs (once per graph)

```bash
pgphase build-snarl-catalog --graph graph.gbz --reference-sample CHM13
```

This writes `graph.catalog/<contig>.pgsites` (or, with `--out DIR`,
`DIR/<contig>.pgsites`) for every reference contig: every snarl
the reference path `CHM13#0#<contig>` passes through, at every nesting level,
with all its alleles, in position order. It is what `vg deconstruct -a` reports,
without the genotypes and allele counts, which need the haplotype names and
take most of its time. `--sites graph.catalog` then gives pgphase every contig.

- `--reference-sample` may be left out when the graph declares a single
  reference sample (`pgcat reference-samples graph.gbz` lists them); HPRC graphs
  declare both CHM13 and GRCh38.
- Contigs may be named either way: `chr20`, or the path name `CHM13#0#chr20`.
- `--region chr20` (repeatable) builds chosen contigs, and
  `--region chr20:10000000-12000000` part of one, as `chr20_10000000-12000000.pgsites`.
- Each contig's snarls come from its own connected component (`vg chunk -C`,
  then `vg snarls`, several contigs at a time) and are kept in the catalog's
  `snarls/` folder for later runs on the same graph.

On the HPRC v2.1 CHM13 graph (10 cores), all 25 contigs take 58 minutes: 12 to
extract the contigs' components (`vg chunk`, 24 GB), 11 for their snarls (up to
`--snarls-jobs` 6 `vg snarls` at once, about 10 GB each) and 34 for the catalogs. The catalogs total 1.2 GB with every site
and allele; chr20's is 30 MB, against 656 MB as a sites-only VCF. Records too
long to hold in memory (chr2 has a snarl spanning 11.4 Mb, a 6.5 GB line) are
streamed from their haplotype walks, so pgcat itself peaks at the graph plus one
contig's working set: 12 GB for chr20, 16 GB for chr2.

A `.pgsites` file holds the catalog's VCF lines in independently
zstd-compressed position blocks with an index, so a region query reads only the
blocks it needs. `pgcat view graph.catalog/chr20.pgsites [chr20:START-END]` prints it,
or the records overlapping a region as `tabix` would, as VCF; pgphase also reads
a bgzipped, tabix-indexed sites VCF (one file for every contig).

### 2. Index the alignments by reference position

```bash
pgphase index-gaf --graph graph.gbz --gaf reads.gaf.gz --reference-sample CHM13
```

This places every alignment on the reference by the nodes of its walk that lie
on the reference sample's paths, and writes `reads.coord.gaf.gz` (sorted,
bgzipped) with its tabix index. Unmapped alignments, and those without a unique
reference interval, are left out and counted. Header lines record the longest
alignment on each contig: pgphase uses them to skip sites no read can traverse
(every allele longer than any alignment), which changes nothing but memory and
time. The GAF may be plain, gzip or bgzip; sorting spills to disk beyond
`--memory` (a quarter of the machine's memory). chr20 of HG002 HiFi (5.7 GB of
GAF) takes 61 s.

### 3. Phase

```bash
pgphase collect-graph-variation \
    --ref chm13.fa \
    --sites graph.catalog \
    --gaf reads.coord.gaf.gz \
    --bam reads.bam \
    -r 'CHM13#0#chr20' \
    -o candidates.tsv \
    --phased-vcf-out phased.vcf \
    --phased-bam-out phased.bam
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
| `--sites PATH` | Site catalogs from `build-snarl-catalog` (the folder, or one `.pgsites` file), or a tabix-indexed sites VCF | required |
| `--gaf FILE` | Reference-indexed GAF from `index-gaf` | required |
| `--bam FILE` | Indexed BAM of the same reads (union gap phasing) | off |
| `-t INT` | Worker threads | all cores |
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
    -o candidates.tsv --phased-vcf-out phased.vcf -b phased.bam
```

Use `--ont` for ONT reads; `--help` lists all options.

## Testing

```bash
make unit-tests   # unit tests
make pgcat-tests  # pgcat unit tests
make check        # end-to-end golden tests on the bundled data in test_data/
```

`make upstream-parity-tests LONGCALLD_ROOT=/path/to/longcallD` compiles
longcallD's original C and checks the ported phasing functions against it.

## Licence

MIT; see `LICENSE`. Third-party code and its licences are listed in
`THIRD_PARTY_NOTICES.md`.
