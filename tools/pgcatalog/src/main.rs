//! pgcatalog: pgphase's inputs from a GBZ pangenome graph: the site catalogs
//! (`build`), and for dromos (github.com/kokyriakidis/dromos), which indexes the
//! reads' alignments by reference position, the reference map (`ref-map`) and
//! the node table (`node-table`).
//!
//! A catalog lists, for one reference contig, every snarl the reference path
//! passes through (at every nesting level), in position order, with all its
//! alleles: what `vg deconstruct -a` reports, minus what needs haplotype names
//! (genotypes, AC/AF/AN/NS), which is most of its run time.

mod build;
mod catalog;
mod nodetable;
mod graph;
mod pgsites;
mod pileup;
mod record;
mod refmap;
mod reference;
mod snarls;
mod surject;

use std::io::Write;

const USAGE: &str = "\
pgcatalog: pgphase's inputs from a GBZ pangenome graph

  pgcatalog build --graph GRAPH.gbz [options]          site catalogs
  pgcatalog ref-map --graph GRAPH.gbz --out FILE      the reference map dromos index
                                                   places alignments by
  pgcatalog node-table --graph GRAPH.gbz --out FILE   node lengths and letters, for
                                                   dromos gam-to-gaf
  pgcatalog surject --graph GRAPH.gbz --gaf READS.gaf --contig CONTIG
                                                   alignments to a reference contig (SAM;
                                                   experimental: pgphase --bam without vg surject)
  pgcatalog pileup --graph GRAPH.gbz --gaf READS.gaf --contig CONTIG
                                                   candidate variants from the reads' graph
                                                   edits (experimental)
  pgcatalog view FILE.pgsites [CONTIG[:START-END]]     a catalog as VCF
  pgcatalog reference-samples GRAPH.gbz                the graph's reference samples

Run `pgcatalog <command> --help` for a command's options.
";

const BUILD_USAGE: &str = "\
pgcatalog build: write site catalogs (DIR/<contig>.pgsites) from a GBZ graph

  pgcatalog build --graph GRAPH.gbz [options]

  --graph FILE              the pangenome graph (GBZ)
  --out DIR                 where catalogs are written: <contig>.pgsites, or
                            <contig>_<start>-<end>.pgsites for a region
                            [GRAPH.catalog, in the current folder]
  --reference-sample NAME   the sample whose paths are the reference: CHM13 for
                            CHM13#0#chr20 [the graph's declared reference
                            sample, if it declares exactly one; HPRC graphs
                            declare CHM13 and GRCh38 and need this]
  --region REGION           a contig (chr20) or part of one
                            (chr20:10000000-12000000, 1-based, inclusive),
                            repeatable [every reference contig]; a region holds
                            the records a tabix query of the whole contig's
                            catalog returns
  -t, --threads N           threads [all cores]
  --snarls-jobs N           vg snarls runs at a time [6, at most --threads];
                            fewer when memory is short (each needs about 20
                            times its component's size)
  --snarls-cache DIR        where each contig's snarls are kept and reused by
                            later runs on the same graph [DIR/snarls]
  --vg PATH                 the vg binary, for vg paths/chunk/snarls [vg]
  --level N                 zstd level, 1-22 [6]; higher is smaller and slower

Every snarl the reference path passes through is kept, at every nesting level,
with all its alleles. `pgcatalog view` prints a catalog as VCF; for a tabix-indexed
copy: pgcatalog view chr20.pgsites | bgzip > chr20.sites.vcf.gz && tabix -p vcf ...
";

const SURJECT_USAGE: &str = "\
pgcatalog surject (experimental): a GAF's alignments as alignments to one reference
contig, from the graph alone, as SAM on standard output

  pgcatalog surject --graph GRAPH.gbz --gaf READS.gaf --contig CONTIG [options]
      | samtools sort -o READS.bam && samtools index READS.bam

  --graph FILE              the GBZ graph the reads were aligned to
  --gaf FILE                their alignments (GAF, or dromos index's output; plain,
                            gzip or bgzip); they need cs, and bq for qualities
  --contig NAME             the reference contig (chr20)
  --reference-sample NAME   the sample whose path is the reference [CHM13]
  --sam-contig NAME         the contig's name in the SAM [<sample>#0#<contig>]
  --graph-placement         place read bases only where the graph puts them: a
                            stretch of more than 50 bases off the reference
                            (another haplotype's or copy's sequence) is not
                            aligned to it, and such read ends are clipped
  --catalog FILE            the contig's catalog (.pgsites): a read spanning one
                            of its snarls (windows as pgcatalog pileup's) is aligned
                            there as its allele, its own letters between the
                            window's anchors against the reference's, spelled as
                            the truth sets spell variants (minimap2 asm5 costs)
  -t, --threads N           threads [all cores]

Each read is rebuilt from its walk and cs string; a base on a node the
reference path visits once takes that node's position. Between placed
stretches, a substitution bubble is aligned base for base and an indel bubble
globally against the reference; read ends off the reference are soft-clipped.
Reads with no base on the contig are left out.
";

const PILEUP_USAGE: &str = "\
pgcatalog pileup (experimental): candidate small variants from a pileup of the reads'
own graph edits (cs mismatches, insertions, deletions), no alignment to the
reference; TSV on standard output

  pgcatalog pileup --graph GRAPH.gbz --gaf READS.gaf --contig CONTIG [options]

  --graph FILE              the GBZ graph the reads were aligned to
  --gaf FILE                their alignments (GAF with cs; plain, gzip or bgzip)
  --contig NAME             the reference contig (chr20)
  --reference-sample NAME   the sample whose path is the reference [CHM13]
  --out-contig NAME         the contig's name in the output [<sample>#0#<contig>]
  --min-alt N               reads an edit needs to be a candidate [2]
  --region START-END        only reads and candidates in this part of the contig
                            (1-based, inclusive)
  --catalog FILE            the contig's catalog (.pgsites): reads' alleles at its
                            snarls (of at most 500 reference bases) from their
                            own sequences, decomposed into SNVs and indels
  --reads FILE              each read's evidence (name, reverse, MAPQ, reference
                            ranges, placed edits, edits off the reference)
  -t, --threads N           threads [all cores]
";

fn pileup_main(args: &[String]) {
    let mut threads = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1);
    let mut o = pileup::PileupOptions {
        graph: String::new(),
        gaf: String::new(),
        contig: String::new(),
        reference_sample: "CHM13".into(),
        out_contig: None,
        min_alt: 2,
        region: None,
        reads_out: None,
        catalog: None,
    };
    let mut i = 0;
    let value = |i: &mut usize| -> String {
        *i += 1;
        args.get(*i).cloned().unwrap_or_else(|| fail(&format!("{} needs a value", args[*i - 1])))
    };
    while i < args.len() {
        match args[i].as_str() {
            "-h" | "--help" => {
                print!("{PILEUP_USAGE}");
                return;
            }
            "--graph" => o.graph = value(&mut i),
            "--gaf" => o.gaf = value(&mut i),
            "--contig" => o.contig = value(&mut i),
            "--reference-sample" => o.reference_sample = value(&mut i),
            "--out-contig" => o.out_contig = Some(value(&mut i)),
            "--region" => {
                let v = value(&mut i);
                let (a, b) = v.split_once('-').unwrap_or_else(|| fail("--region: START-END (1-based, inclusive)"));
                let num = |x: &str| x.replace(',', "").parse::<u64>().unwrap_or_else(|_| fail("--region: not a number"));
                o.region = Some((num(a) - 1, num(b)));
            }
            "--reads" => o.reads_out = Some(value(&mut i)),
            "--catalog" => o.catalog = Some(value(&mut i)),
            "--min-alt" => o.min_alt = value(&mut i).parse().unwrap_or_else(|_| fail("--min-alt: not a number")),
            "-t" | "--threads" => threads = value(&mut i).parse::<usize>().unwrap_or_else(|_| fail("--threads: not a number")).max(1),
            a => fail(&format!("unknown option {a}; see pgcatalog pileup --help")),
        }
        i += 1;
    }
    if o.graph.is_empty() || o.gaf.is_empty() || o.contig.is_empty() {
        eprint!("{PILEUP_USAGE}");
        std::process::exit(1);
    }
    rayon::ThreadPoolBuilder::new().num_threads(threads).build_global().ok();
    if let Err(e) = pileup::run(&o) {
        fail(&e);
    }
}

fn surject_main(args: &[String]) {
    let mut threads = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1);
    let mut o = surject::SurjectOptions {
        graph: String::new(),
        gaf: String::new(),
        contig: String::new(),
        reference_sample: "CHM13".into(),
        sam_contig: None,
        graph_placement: false,
        catalog: None,
    };
    let mut i = 0;
    let value = |i: &mut usize| -> String {
        *i += 1;
        args.get(*i).cloned().unwrap_or_else(|| fail(&format!("{} needs a value", args[*i - 1])))
    };
    while i < args.len() {
        match args[i].as_str() {
            "-h" | "--help" => {
                print!("{SURJECT_USAGE}");
                return;
            }
            "--graph" => o.graph = value(&mut i),
            "--gaf" => o.gaf = value(&mut i),
            "--contig" => o.contig = value(&mut i),
            "--reference-sample" => o.reference_sample = value(&mut i),
            "--sam-contig" => o.sam_contig = Some(value(&mut i)),
            "--graph-placement" => o.graph_placement = true,
            "--catalog" => o.catalog = Some(value(&mut i)),
            "-t" | "--threads" => threads = value(&mut i).parse::<usize>().unwrap_or_else(|_| fail("--threads: not a number")).max(1),
            a => fail(&format!("unknown option {a}; see pgcatalog surject --help")),
        }
        i += 1;
    }
    if o.graph.is_empty() || o.gaf.is_empty() || o.contig.is_empty() {
        eprint!("{SURJECT_USAGE}");
        std::process::exit(1);
    }
    rayon::ThreadPoolBuilder::new().num_threads(threads).build_global().ok();
    if let Err(e) = surject::run(&o) {
        fail(&e);
    }
}

/// A quarter of the machine's memory, from /proc/meminfo; 4 GB elsewhere.
const NODE_TABLE_USAGE: &str = "\
pgcatalog node-table: a GBZ graph's node lengths and letters as the node table dromos reads GAM with

usage: pgcatalog node-table --graph GRAPH.gbz --out FILE [--no-letters]

  --no-letters              lengths alone (enough for GAF's columns, not for its cs)
";

fn node_table_main(args: &[String]) {
    let (mut graph, mut out, mut letters) = (String::new(), String::new(), true);
    let mut i = 0;
    let value = |i: &mut usize| -> String {
        *i += 1;
        args.get(*i).cloned().unwrap_or_else(|| fail(&format!("{} needs a value", args[*i - 1])))
    };
    while i < args.len() {
        match args[i].as_str() {
            "-h" | "--help" => {
                print!("{NODE_TABLE_USAGE}");
                return;
            }
            "--graph" => graph = value(&mut i),
            "--out" => out = value(&mut i),
            "--no-letters" => letters = false,
            a => fail(&format!("unknown option {a}; see pgcatalog node-table --help")),
        }
        i += 1;
    }
    if graph.is_empty() || out.is_empty() {
        eprint!("{NODE_TABLE_USAGE}");
        std::process::exit(1);
    }
    let g: gbz::GBZ = simple_sds::serialize::load_from(&graph).unwrap_or_else(|e| fail(&format!("{graph}: {e}")));
    if let Err(e) = nodetable::write(&g, &out, letters) {
        fail(&e);
    }
}

const REF_MAP_USAGE: &str = "\
pgcatalog ref-map: where each node of a GBZ graph lies on its reference sample's paths, saved once per graph for
dromos index, which places alignments on the reference by it

usage: pgcatalog ref-map --graph GRAPH.gbz --out FILE [--reference-sample NAME]
";

fn ref_map_main(args: &[String]) {
    let (mut graph, mut out, mut sample) = (String::new(), String::new(), None);
    let mut i = 0;
    let value = |i: &mut usize| -> String {
        *i += 1;
        args.get(*i).cloned().unwrap_or_else(|| fail(&format!("{} needs a value", args[*i - 1])))
    };
    while i < args.len() {
        match args[i].as_str() {
            "-h" | "--help" => {
                print!("{REF_MAP_USAGE}");
                return;
            }
            "--graph" => graph = value(&mut i),
            "--out" => out = value(&mut i),
            "--reference-sample" => sample = Some(value(&mut i)),
            a => fail(&format!("unknown option {a}; see pgcatalog ref-map --help")),
        }
        i += 1;
    }
    if graph.is_empty() || out.is_empty() {
        eprint!("{REF_MAP_USAGE}");
        std::process::exit(1);
    }
    if let Err(e) = refmap::save_ref_map(&graph, sample.as_deref(), &out) {
        fail(&e);
    }
}

fn fail(msg: &str) -> ! {
    eprintln!("error: {msg}");
    std::process::exit(2);
}

fn build_main(args: &[String]) {
    let threads = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1);
    let mut o = build::BuildOptions {
        graph: String::new(),
        out: String::new(),
        reference_sample: None,
        regions: Vec::new(),
        snarls_jobs: None,
        snarls_cache: None,
        vg: "vg".into(),
        catalog: catalog::CatalogOptions { threads, level: 6 },
    };
    let mut i = 0;
    let value = |i: &mut usize| -> String {
        *i += 1;
        args.get(*i).cloned().unwrap_or_else(|| fail(&format!("{} needs a value", args[*i - 1])))
    };
    let number = |s: String, what: &str| -> usize { s.parse().unwrap_or_else(|_| fail(&format!("{what}: not a number: {s}"))) };
    while i < args.len() {
        match args[i].as_str() {
            "-h" | "--help" => {
                print!("{BUILD_USAGE}");
                return;
            }
            "--graph" => o.graph = value(&mut i),
            "--out" => o.out = value(&mut i),
            "--reference-sample" => o.reference_sample = Some(value(&mut i)),
            "--region" => o.regions.push(value(&mut i)),
            "-t" | "--threads" => o.catalog.threads = number(value(&mut i), "--threads").max(1),
            "--snarls-jobs" => o.snarls_jobs = Some(number(value(&mut i), "--snarls-jobs").max(1)),
            "--snarls-cache" => o.snarls_cache = Some(value(&mut i)),
            "--vg" => o.vg = value(&mut i),
            "--level" => o.catalog.level = (number(value(&mut i), "--level") as i32).clamp(1, 22),
            a => fail(&format!("unknown option {a}; see pgcatalog build --help")),
        }
        i += 1;
    }
    if o.out.is_empty() && !o.graph.is_empty() {
        let stem = std::path::Path::new(&o.graph).file_stem().map(|s| s.to_string_lossy().into_owned());
        o.out = format!("{}.catalog", stem.unwrap_or_else(|| "pgcatalog".into()));
    }
    if o.graph.is_empty() {
        eprint!("{BUILD_USAGE}");
        std::process::exit(1);
    }
    rayon::ThreadPoolBuilder::new().num_threads(o.catalog.threads).build_global().unwrap();
    if let Err(e) = build::run(&o) {
        fail(&e);
    }
}

/// `chr20` names the contig of any sample's path; `CHM13#0#chr20` only CHM13's.
fn same_contig(a: &str, b: &str) -> bool {
    let split = |s: &str| match (s.find('#'), s.rfind('#')) {
        (Some(i), Some(j)) => (Some(s[..i].to_string()), s[j + 1..].to_string()),
        _ => (None, s.to_string()),
    };
    let ((sa, ca), (sb, cb)) = (split(a), split(b));
    ca == cb && (sa.is_none() || sb.is_none() || sa == sb)
}

/// The header, then every record (or those overlapping the region), as VCF.
fn view_main(args: &[String]) {
    if args.is_empty() || args[0] == "-h" || args[0] == "--help" {
        println!("pgcatalog view FILE.pgsites [CONTIG[:START-END]]: print a catalog, or the records\n\
                  overlapping a region (as tabix would), as VCF");
        return;
    }
    let r = pgsites::Reader::open(&args[0]).unwrap_or_else(|e| fail(&e.to_string()));
    let (mut beg, mut end, mut contig) = (1usize, usize::MAX, None);
    if let Some(q) = args.get(1) {
        let (c, span) = match q.rsplit_once(':') {
            Some((c, s)) => (c, Some(s)),
            None => (q.as_str(), None),
        };
        contig = Some(c.to_string());
        if let Some(s) = span {
            let (a, b) = s.split_once('-').unwrap_or_else(|| fail("region: CONTIG:START-END"));
            beg = a.replace(',', "").parse().unwrap_or_else(|_| fail("region start"));
            end = b.replace(',', "").parse().unwrap_or_else(|_| fail("region end"));
        }
    }
    let out = std::io::stdout();
    let mut w = std::io::BufWriter::with_capacity(1 << 24, out.lock());
    let _ = w.write_all(r.header.as_bytes());
    if let Some(c) = contig {
        if !same_contig(&c, &r.chrom) {
            return;
        }
    }
    let res = r.for_each_record(beg, end, |line| {
        let _ = w.write_all(line);
        let _ = w.write_all(b"\n");
    });
    if let Err(e) = res {
        fail(&e.to_string());
    }
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match args.first().map(|s| s.as_str()) {
        Some("build") => build_main(&args[1..]),
        Some("view") => view_main(&args[1..]),
        Some("ref-map") => ref_map_main(&args[1..]),
        Some("node-table") => node_table_main(&args[1..]),
        Some("surject") => surject_main(&args[1..]),
        Some("pileup") => pileup_main(&args[1..]),
        Some("reference-samples") => {
            let g = args.get(1).unwrap_or_else(|| fail("usage: pgcatalog reference-samples GRAPH.gbz"));
            for s in graph::reference_samples(g).unwrap_or_else(|e| fail(&e)) {
                println!("{s}");
            }
        }
        Some("-h") | Some("--help") | None => print!("{USAGE}"),
        Some(c) => fail(&format!("unknown command {c}\n{USAGE}")),
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn same_contig() {
        use super::same_contig;
        assert!(same_contig("chr20", "CHM13#0#chr20"));
        assert!(same_contig("CHM13#0#chr20", "chr20"));
        assert!(same_contig("CHM13#0#chr20", "CHM13#0#chr20"));
        assert!(!same_contig("GRCh38#0#chr20", "CHM13#0#chr20"));
        assert!(!same_contig("chr2", "CHM13#0#chr20"));
    }
}
