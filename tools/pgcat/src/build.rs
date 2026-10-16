//! `pgcat build`: catalogs for every reference contig of a graph, for chosen
//! contigs, or for regions of them.
//!
//! 0. The reference paths come from the GBZ's path metadata, read alone.
//! 1. Snarls come from each contig's connected component alone: one `vg chunk -C`
//!    pass extracts the components of the contigs not yet in the snarls cache,
//!    and `vg snarls` runs on several at once (it uses about one core each).
//!    Whole-genome snarls would make every contig pay for the genome: vg
//!    enumerates a snarl's haplotypes before checking which contig it is on.
//! 2. One GBZ load, then each selection is built (see catalog.rs).
//!
//! The cache records the graph's fingerprint and refuses another graph's snarls.

use crate::catalog::{build_catalog, CatalogOptions};
use crate::graph::{fingerprint, reference_paths, resolve_reference_sample};
use crate::reference::ReferencePath;
use gbz::GBZ;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::process::{Child, Command, Stdio};
use std::time::Instant;

pub struct BuildOptions {
    pub graph: String,
    pub out: String,
    pub reference_sample: Option<String>,
    pub regions: Vec<String>,
    /// vg snarls runs at a time; None for the default.
    pub snarls_jobs: Option<usize>,
    pub snarls_cache: Option<String>,
    pub vg: String,
    pub catalog: CatalogOptions,
}

struct Selection {
    contig: String,
    region: Option<(usize, usize)>,
}

fn log(t0: &Instant, msg: &str) {
    eprintln!("[pgcat {:7.1}s] {msg}", t0.elapsed().as_secs_f64());
}

/// `chr20`, `chr20:10,000,000-12,000,000` or `CHM13#0#chr20:...`.
fn parse_region(r: &str, sample: &str) -> Result<Selection, String> {
    let prefix = format!("{sample}#0#");
    let r = r.strip_prefix(&prefix).unwrap_or(r);
    if r.contains('#') {
        return Err(format!("region {r} is not on the reference sample {sample}'s paths"));
    }
    match r.rsplit_once(':') {
        None => Ok(Selection { contig: r.to_string(), region: None }),
        Some((contig, span)) => {
            let num = |s: &str| s.replace(',', "").parse::<usize>().map_err(|_| format!("bad region {r}"));
            let (a, b) = span.split_once('-').ok_or_else(|| format!("bad region {r}: use contig:start-end"))?;
            let (a, b) = (num(a)?, num(b)?);
            if a == 0 || b < a {
                return Err(format!("bad region {r}: start-end are 1-based and inclusive"));
            }
            Ok(Selection { contig: contig.to_string(), region: Some((a, b)) })
        }
    }
}

/// chr1 < chr2 < ... < chr10 < chrX: digits compared as numbers.
fn natural_key(s: &str) -> Vec<(String, u64)> {
    let mut out = Vec::new();
    let mut text = String::new();
    let mut num = String::new();
    for c in s.chars() {
        if c.is_ascii_digit() {
            num.push(c);
        } else {
            if !num.is_empty() {
                out.push((std::mem::take(&mut text), num.parse().unwrap_or(0)));
                num.clear();
            }
            text.push(c);
        }
    }
    out.push((text, num.parse().unwrap_or(u64::MAX)));
    out
}

fn wait_one(running: &mut Vec<(Child, String)>) -> Result<(), String> {
    loop {
        for i in 0..running.len() {
            if let Some(status) = running[i].0.try_wait().map_err(|e| e.to_string())? {
                let (_, what) = running.swap_remove(i);
                return if status.success() { Ok(()) } else { Err(format!("{what} failed")) };
            }
        }
        std::thread::sleep(std::time::Duration::from_millis(200));
    }
}

/// Removes the `<prefix><pid>` work folders in `dir` of earlier runs that were
/// killed (their process is gone), which would otherwise stay indefinitely:
/// a build's contig graphs, an index-gaf's sorted runs. Needs /proc to tell;
/// elsewhere they are left alone.
pub fn remove_abandoned_work(dir: &str, prefix: &str) {
    if !Path::new("/proc/self").exists() {
        return;
    }
    let Ok(entries) = std::fs::read_dir(dir) else { return };
    for e in entries.flatten() {
        let name = e.file_name().to_string_lossy().into_owned();
        let Some(pid) = name.strip_prefix(prefix).and_then(|p| p.parse::<u32>().ok()) else { continue };
        if !Path::new(&format!("/proc/{pid}")).exists() {
            eprintln!("[pgcat] removing {}, left by an interrupted run", e.path().display());
            let _ = std::fs::remove_dir_all(e.path());
        }
    }
}

/// Computes the snarls of `contigs` into the cache, in a work folder under
/// --out that is removed afterwards, also when a step fails (it holds the
/// contigs' graphs: 24 GB for HPRC v2.1).
fn compute_snarls(
    t0: &Instant,
    o: &BuildOptions,
    sample: &str,
    paths: &BTreeMap<String, Vec<String>>,
    contigs: &[String],
    cache: &Path,
) -> Result<(), String> {
    remove_abandoned_work(&o.out, ".pgcat-work-");
    let work = PathBuf::from(&o.out).join(format!(".pgcat-work-{}", std::process::id()));
    std::fs::create_dir_all(&work).map_err(|e| format!("{}: {e}", work.display()))?;
    let result = compute_snarls_in(t0, o, sample, paths, contigs, cache, &work);
    let _ = std::fs::remove_dir_all(&work);
    result
}

fn compute_snarls_in(
    t0: &Instant,
    o: &BuildOptions,
    sample: &str,
    paths: &BTreeMap<String, Vec<String>>,
    contigs: &[String],
    cache: &Path,
    work: &Path,
) -> Result<(), String> {
    let path_list = work.join("paths.txt");
    let names: Vec<&String> = contigs.iter().flat_map(|c| paths[c].iter()).collect();
    std::fs::write(&path_list, names.iter().map(|n| format!("{n}\n")).collect::<String>()).map_err(|e| e.to_string())?;
    log(t0, &format!("extracting {} contig component(s) with vg chunk -C", contigs.len()));
    let prefix = work.join("comp");
    let status = Command::new(&o.vg)
        .args(["chunk", "-x", &o.graph, "-C", "-P"]).arg(&path_list)
        .args(["-O", "pg", "-t", &o.catalog.threads.to_string(), "-b"]).arg(&prefix)
        .stdout(Stdio::null()).status().map_err(|e| format!("running {} chunk: {e}", o.vg))?;
    if !status.success() {
        return Err(format!("{} chunk failed", o.vg));
    }
    // One component file per targeted path; fragments of a contig that share a
    // component may come out in one file.
    let mut jobs: Vec<(String, PathBuf, PathBuf)> = Vec::new(); // contig, component, snarls part
    for c in contigs {
        let parts: Vec<PathBuf> = paths[c].iter()
            .map(|n| PathBuf::from(format!("{}_{n}.vg", prefix.display())))
            .filter(|p| p.exists())
            .collect();
        if parts.is_empty() {
            return Err(format!("vg chunk wrote no component for {sample}#0#{c}"));
        }
        for (k, p) in parts.into_iter().enumerate() {
            jobs.push((c.clone(), p, work.join(format!("{c}.{k}.snarls"))));
        }
    }
    // vg snarls uses about one core and, on HPRC chr1, 10 GB: at most 6 at
    // once by default, within --threads.
    let snarls_jobs = o.snarls_jobs.unwrap_or(6).min(o.catalog.threads).max(1);
    log(t0, &format!("computing snarls of {} component(s), {} at a time", jobs.len(), snarls_jobs.min(jobs.len())));
    let per_job = if o.catalog.threads >= 2 * snarls_jobs { 2 } else { 1 };
    let mut running: Vec<(Child, String)> = Vec::new();
    for (c, comp, part) in &jobs {
        if running.len() >= snarls_jobs {
            wait_one(&mut running)?;
        }
        let out = std::fs::File::create(part).map_err(|e| e.to_string())?;
        let child = Command::new(&o.vg).args(["snarls", "-t", &per_job.to_string()]).arg(comp).stdout(out).spawn()
            .map_err(|e| format!("running {} snarls: {e}", o.vg))?;
        running.push((child, format!("vg snarls on {c}")));
    }
    while !running.is_empty() {
        wait_one(&mut running)?;
    }
    // A contig's snarls file is its components' files one after another (both
    // the stream format and gzip allow concatenation).
    for c in contigs {
        let tmp = cache.join(format!("{c}.snarls.tmp"));
        let mut out = std::fs::File::create(&tmp).map_err(|e| e.to_string())?;
        for (_, _, part) in jobs.iter().filter(|j| &j.0 == c) {
            std::io::copy(&mut std::fs::File::open(part).map_err(|e| e.to_string())?, &mut out).map_err(|e| e.to_string())?;
        }
        std::fs::rename(&tmp, cache.join(format!("{c}.snarls"))).map_err(|e| e.to_string())?;
    }
    Ok(())
}

pub fn run(o: &BuildOptions) -> Result<(), String> {
    let t0 = Instant::now();
    let sample = resolve_reference_sample(&o.graph, o.reference_sample.as_deref())?;
    if o.reference_sample.is_none() {
        log(&t0, &format!("reference sample {sample} (the only one the graph declares)"));
    }
    std::fs::create_dir_all(&o.out).map_err(|e| format!("{}: {e}", o.out))?;
    let cache = PathBuf::from(o.snarls_cache.clone().unwrap_or_else(|| format!("{}/snarls", o.out)));
    std::fs::create_dir_all(&cache).map_err(|e| format!("{}: {e}", cache.display()))?;
    let id = fingerprint(&o.graph)?;
    let id_file = cache.join("graph.id");
    match std::fs::read_to_string(&id_file) {
        Ok(cached) if cached.trim() != id => {
            return Err(format!("{} holds snarls of another graph ({}); use another --snarls-cache", cache.display(), cached.trim()));
        }
        Ok(_) => {}
        Err(_) => std::fs::write(&id_file, format!("{id}\n")).map_err(|e| e.to_string())?,
    }

    let mut selections: Vec<Selection> =
        o.regions.iter().map(|r| parse_region(r, &sample)).collect::<Result<_, _>>()?;
    let mut contigs: Vec<String> = selections.iter().map(|s| s.contig.clone()).collect();
    contigs.sort();
    contigs.dedup();
    let cached = |c: &str| cache.join(format!("{c}.snarls")).is_file();
    let need_paths = selections.is_empty() || contigs.iter().any(|c| !cached(c));
    if need_paths {
        let paths = reference_paths(&o.graph, &sample)?;
        if paths.is_empty() {
            return Err(format!("the graph has no {sample}#0#<contig> reference paths"));
        }
        if selections.is_empty() {
            let mut all: Vec<String> = paths.keys().cloned().collect();
            all.sort_by_key(|c| natural_key(c));
            selections = all.iter().map(|c| Selection { contig: c.clone(), region: None }).collect();
            contigs = all;
        }
        if let Some(c) = contigs.iter().find(|c| !paths.contains_key(*c)) {
            return Err(format!("no reference path {sample}#0#{c} in the graph"));
        }
        let missing: Vec<String> = contigs.iter().filter(|c| !cached(c)).cloned().collect();
        if !missing.is_empty() {
            compute_snarls(&t0, o, &sample, &paths, &missing, &cache)?;
        }
    }
    log(&t0, &format!("{} selection(s) on {} contig(s); snarls in {}", selections.len(), contigs.len(), cache.display()));

    let g: GBZ = simple_sds::serialize::load_from(&o.graph).map_err(|e| format!("{}: {e}", o.graph))?;
    log(&t0, "graph loaded");
    let mut reference: Option<ReferencePath> = None;
    for s in &selections {
        if reference.as_ref().map(|r| r.chrom != format!("{sample}#0#{}", s.contig)).unwrap_or(true) {
            reference = Some(ReferencePath::load(&g, &sample, &s.contig)?);
        }
        let reference = reference.as_ref().unwrap();
        let prefix = match s.region {
            Some((a, b)) => format!("{}/{}_{a}-{b}", o.out, s.contig),
            None => format!("{}/{}", o.out, s.contig),
        };
        let snarls = cache.join(format!("{}.snarls", s.contig));
        let sm = build_catalog(&g, reference, &snarls.to_string_lossy(), s.region, &prefix, &o.catalog)?;
        log(&t0, &format!(
            "{prefix}.pgsites: {} records from {} snarls ({} on the reference{}), {} blocks, {:.1}s",
            sm.records, sm.snarls, sm.on_reference, if s.region.is_some() { " and the region" } else { "" },
            sm.blocks, sm.seconds,
        ));
    }
    log(&t0, "done");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn regions() {
        let s = parse_region("chr20", "CHM13").unwrap();
        assert_eq!((s.contig.as_str(), s.region), ("chr20", None));
        let s = parse_region("CHM13#0#chr20:10,000,000-12,000,000", "CHM13").unwrap();
        assert_eq!((s.contig.as_str(), s.region), ("chr20", Some((10_000_000, 12_000_000))));
        assert!(parse_region("chr20:0-5", "CHM13").is_err());
        assert!(parse_region("chr20:9-5", "CHM13").is_err());
        assert!(parse_region("chr20:5", "CHM13").is_err());
        assert!(parse_region("GRCh38#0#chr20", "CHM13").is_err());
    }

    #[test]
    fn contigs_in_natural_order() {
        let mut v = vec!["chrX", "chr10", "chr2", "chr1", "chrM", "chr22", "chrY"];
        v.sort_by_key(|c| natural_key(c));
        assert_eq!(v, ["chr1", "chr2", "chr10", "chr22", "chrM", "chrX", "chrY"]);
    }
}
