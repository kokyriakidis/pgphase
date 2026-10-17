//! pgcat: pgphase's inputs from a GBZ pangenome graph: the site catalogs
//! (`build`) and the reads' alignments indexed by reference position
//! (`index-gaf`).
//!
//! A catalog lists, for one reference contig, every snarl the reference path
//! passes through (at every nesting level), in position order, with all its
//! alleles: what `vg deconstruct -a` reports, minus what needs haplotype names
//! (genotypes, AC/AF/AN/NS), which is most of its run time.

mod bgzf;
mod build;
mod catalog;
mod gaf;
mod graph;
mod pgsites;
mod record;
mod reference;
mod snarls;
mod tabix;

use std::io::Write;

const USAGE: &str = "\
pgcat: pgphase's inputs from a GBZ pangenome graph

  pgcat build --graph GRAPH.gbz [options]          site catalogs
  pgcat index-gaf --graph GRAPH.gbz --gaf READS.gaf [options]
                                                   alignments by reference position
  pgcat view FILE.pgsites [CONTIG[:START-END]]     a catalog as VCF
  pgcat reference-samples GRAPH.gbz                the graph's reference samples

Run `pgcat <command> --help` for a command's options.
";

const BUILD_USAGE: &str = "\
pgcat build: write site catalogs (DIR/<contig>.pgsites) from a GBZ graph

  pgcat build --graph GRAPH.gbz [options]

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
with all its alleles. `pgcat view` prints a catalog as VCF; for a tabix-indexed
copy: pgcat view chr20.pgsites | bgzip > chr20.sites.vcf.gz && tabix -p vcf ...
";

const INDEX_GAF_USAGE: &str = "\
pgcat index-gaf: index a GAF by reference position, for pgphase --gaf

  pgcat index-gaf --graph GRAPH.gbz --gaf READS.gaf [options]

  --graph FILE              the GBZ graph the reads were aligned to
  --gaf FILE                their alignments (GAF; plain, gzip or bgzip)
  --out FILE                [READS.coord.gaf.gz, in the current folder], with
                            FILE.tbi beside it
  --reference-sample NAME   the sample whose paths give the coordinates [the
                            graph's declared reference sample, if it declares
                            exactly one]
  -t, --threads N           threads [all cores]
  --memory GB               alignments sorted in memory before spilling sorted
                            runs to disk [a quarter of the machine's memory]
  --tmp-dir DIR             where sorted runs are spilled [the output's folder]
  --level N                 BGZF compression level, 0-9 [6]

Each alignment is placed on the reference by the walk's nodes on the reference
sample's paths and written as `contig  start  end  <GAF line>` (1-based,
inclusive), sorted, bgzipped and tabix-indexed. Unmapped alignments, and those
without a unique reference interval (all off the reference, or ambiguous), are
left out and counted.
";

/// A quarter of the machine's memory, from /proc/meminfo; 4 GB elsewhere.
fn default_memory() -> usize {
    std::fs::read_to_string("/proc/meminfo")
        .ok()
        .and_then(|m| m.lines().find(|l| l.starts_with("MemTotal:")).and_then(|l| l.split_whitespace().nth(1)?.parse::<usize>().ok()))
        .map(|kb| kb * 1024 / 4)
        .unwrap_or(4 << 30)
}

fn index_gaf_main(args: &[String]) {
    let mut threads = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1);
    let mut o = gaf::IndexGafOptions {
        graph: String::new(),
        gaf: String::new(),
        out: String::new(),
        reference_sample: None,
        memory: default_memory(),
        tmp_dir: None,
        level: 6,
    };
    let mut i = 0;
    let value = |i: &mut usize| -> String {
        *i += 1;
        args.get(*i).cloned().unwrap_or_else(|| fail(&format!("{} needs a value", args[*i - 1])))
    };
    let number = |s: String, what: &str| -> f64 { s.parse().unwrap_or_else(|_| fail(&format!("{what}: not a number: {s}"))) };
    while i < args.len() {
        match args[i].as_str() {
            "-h" | "--help" => {
                print!("{INDEX_GAF_USAGE}");
                return;
            }
            "--graph" => o.graph = value(&mut i),
            "--gaf" => o.gaf = value(&mut i),
            "--out" => o.out = value(&mut i),
            "--reference-sample" => o.reference_sample = Some(value(&mut i)),
            "-t" | "--threads" => threads = (number(value(&mut i), "--threads") as usize).max(1),
            "--memory" => o.memory = ((number(value(&mut i), "--memory") * (1u64 << 30) as f64) as usize).max(1 << 26),
            "--tmp-dir" => o.tmp_dir = Some(value(&mut i)),
            "--level" => o.level = (number(value(&mut i), "--level") as u32).min(9),
            a => fail(&format!("unknown option {a}; see pgcat index-gaf --help")),
        }
        i += 1;
    }
    if o.graph.is_empty() || o.gaf.is_empty() {
        eprint!("{INDEX_GAF_USAGE}");
        std::process::exit(1);
    }
    if o.out.is_empty() {
        let name = std::path::Path::new(&o.gaf).file_name().map(|s| s.to_string_lossy().into_owned()).unwrap_or_default();
        let stem = name.trim_end_matches(".gz").trim_end_matches(".bgz").trim_end_matches(".gaf");
        o.out = format!("{stem}.coord.gaf.gz");
    }
    rayon::ThreadPoolBuilder::new().num_threads(threads).build_global().unwrap();
    let sample = graph::resolve_reference_sample(&o.graph, o.reference_sample.as_deref()).unwrap_or_else(|e| fail(&e));
    if let Err(e) = gaf::run(&o, &sample) {
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
            a => fail(&format!("unknown option {a}; see pgcat build --help")),
        }
        i += 1;
    }
    if o.out.is_empty() && !o.graph.is_empty() {
        let stem = std::path::Path::new(&o.graph).file_stem().map(|s| s.to_string_lossy().into_owned());
        o.out = format!("{}.catalog", stem.unwrap_or_else(|| "pgcat".into()));
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
        println!("pgcat view FILE.pgsites [CONTIG[:START-END]]: print a catalog, or the records\n\
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
        Some("index-gaf") => index_gaf_main(&args[1..]),
        Some("reference-samples") => {
            let g = args.get(1).unwrap_or_else(|| fail("usage: pgcat reference-samples GRAPH.gbz"));
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
