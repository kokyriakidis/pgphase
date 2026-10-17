//! Building one contig's catalog, or one region of it, from its snarls.
//!
//! Pass 1 finds each snarl's reference traversal; snarls the reference does not
//! pass through have no position and no record (vg deconstruct skips them
//! too). The rest are sorted by position, ancestors first at a tie, and built
//! in batches: in parallel, then LV/PS from the ancestors already built, then
//! written in POS order by a writer thread while the next batch builds. Only one
//! or two batches are ever in memory.
//!
//! One snarl's walk search is serial and a batch waits for its slowest snarl,
//! so the few snarls spanning tens of kilobases to megabases (chr2 of HPRC
//! v2.1: 11.4 Mb, 30 s) are built ahead of their turn on helper threads, from
//! the start. A helper also writes and compresses a giant record's line once
//! its ancestors are known, so the writer only copies the block.

use crate::pgsites;
use crate::record::{build_record, RefTrav, Record};
use crate::reference::ReferencePath;
use crate::snarls::{canonical, read_snarls, snarl_name};
use gbz::GBZ;
use rayon::prelude::*;
use std::collections::HashMap;
use std::sync::atomic::{AtomicBool, AtomicU8, AtomicUsize, Ordering};
use std::sync::{Condvar, Mutex};
use std::time::Instant;

/// Snarls built per batch: enough to keep every thread busy, few enough that
/// one or two batches of records fit in memory.
const BATCH: usize = 50_000;

/// Records with longer lines are not formatted in memory: they are streamed
/// from their walks into a block of their own.
const STREAM_BYTES: usize = 64 << 20;

/// Snarls spanning at least this much of the reference are heavy: built ahead
/// of their turn on helper threads. A snarl's ancestors span at least as much,
/// so the ancestors of a heavy snarl are heavy too.
const HEAVY_SPAN: usize = 50_000;

/// What building a heavy snarl gives: its record, or, for a giant one, the
/// compressed block it becomes.
enum Built {
    Record(Option<Record>),
    Block(pgsites::Block),
}

/// A heavy snarl's result, built by a helper or by the batch reaching it
/// first, whichever claims it.
struct Heavy {
    /// 0 waiting, 1 being built, 2 built.
    state: AtomicU8,
    /// Whether it has a record (once built).
    has_record: AtomicBool,
    out: Mutex<Option<Built>>,
    built: Condvar,
}

impl Heavy {
    fn new() -> Self {
        Heavy { state: AtomicU8::new(0), has_record: AtomicBool::new(false), out: Mutex::new(None), built: Condvar::new() }
    }

    /// Builds the result if no one has claimed it.
    fn try_build(&self, build: impl FnOnce() -> Built) {
        if self.state.compare_exchange(0, 1, Ordering::AcqRel, Ordering::Acquire).is_ok() {
            let b = build();
            self.has_record.store(!matches!(b, Built::Record(None)), Ordering::Release);
            *self.out.lock().unwrap() = Some(b);
            self.state.store(2, Ordering::Release);
            self.built.notify_all();
        }
    }

    fn wait(&self) -> std::sync::MutexGuard<'_, Option<Built>> {
        let mut out = self.out.lock().unwrap();
        while self.state.load(Ordering::Acquire) != 2 {
            out = self.built.wait(out).unwrap();
        }
        out
    }

    /// Whether it has a record, waiting until it is built.
    fn emitted(&self) -> bool {
        drop(self.wait());
        self.has_record.load(Ordering::Acquire)
    }

    /// The result: built here if no one has started it, else awaited.
    fn take(&self, build: impl FnOnce() -> Built) -> Built {
        self.try_build(build);
        self.wait().take().expect("heavy snarl taken twice")
    }
}

/// What the writer receives, in file order.
enum Item {
    Line { pos: usize, end: usize, min_allele: usize, max_walk: usize, line: Vec<u8> },
    Giant { record: Record, lv: usize, ps: Option<String>, len: usize },
    Block(pgsites::Block),
}

impl Item {
    fn pos(&self) -> usize {
        match self {
            Item::Line { pos, .. } => *pos,
            Item::Giant { record, .. } => record.pos,
            Item::Block(b) => b.pos,
        }
    }
}

pub struct CatalogOptions {
    pub threads: usize,
    /// zstd level of the blocks.
    pub level: i32,
}

pub fn vcf_header(chrom: &str, length: usize) -> String {
    format!(
        "##fileformat=VCFv4.2\n\
         ##INFO=<ID=AT,Number=R,Type=String,Description=\"Allele Traversal as path in graph\">\n\
         ##INFO=<ID=LV,Number=1,Type=Integer,Description=\"Level in the snarl tree (0=top level)\">\n\
         ##INFO=<ID=PS,Number=1,Type=String,Description=\"ID of variant corresponding to parent snarl\">\n\
         ##contig=<ID={chrom},length={length}>\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )
}

pub struct Summary {
    pub snarls: usize,
    pub on_reference: usize,
    pub records: usize,
    pub blocks: usize,
    pub seconds: f64,
}

/// LV / PS, as vg's update_nesting_info_tags: the ancestors that have a record
/// count, the nearest names the parent. None if `emitted` cannot tell for an
/// ancestor.
fn nesting(i: usize, parent: &[u32], snarls: &[crate::snarls::Snarl], emitted: impl Fn(usize) -> Option<bool>) -> Option<(usize, Option<String>)> {
    let (mut lv, mut ps) = (0usize, None);
    let mut cur = parent[i];
    while cur != u32::MAX {
        let pi = cur as usize;
        if emitted(pi)? {
            lv += 1;
            if ps.is_none() {
                ps = Some(snarl_name(snarls[pi].start, snarls[pi].end));
            }
        }
        cur = parent[pi];
    }
    Some((lv, ps))
}

/// `region`: 1-based inclusive; only records whose [POS, POS + len(REF) - 1]
/// overlaps it are written, as `tabix contig:start-end` would return them.
pub fn build_catalog(
    g: &GBZ,
    reference: &ReferencePath,
    snarls_path: &str,
    region: Option<(usize, usize)>,
    out_prefix: &str,
    opts: &CatalogOptions,
) -> Result<Summary, String> {
    let t0 = Instant::now();
    let snarls = read_snarls(snarls_path).map_err(|e| format!("{snarls_path}: {e}"))?;
    let parent: Vec<u32> = {
        let idx: HashMap<_, usize> =
            snarls.iter().enumerate().map(|(i, s)| (canonical((s.start, s.end)), i)).collect();
        snarls
            .iter()
            .map(|s| s.parent.and_then(|p| idx.get(&canonical(p))).map(|&i| i as u32).unwrap_or(u32::MAX))
            .collect()
    };

    // Pass 1: reference traversals, kept when there is exactly one (a snarl the
    // reference crosses twice would need one record per crossing; none of the
    // HPRC v2.1 CHM13 contigs has one) and, for a region, when it can overlap.
    let mut refs: Vec<RefTrav> = snarls
        .par_iter()
        .enumerate()
        .filter_map(|(i, sn)| {
            let mut t = reference.traversals(g, sn.start, sn.end);
            if t.len() != 1 {
                return None;
            }
            let (trav, s0, e0) = t.pop().unwrap();
            let (left, right) = if reference.pos[s0] <= reference.pos[e0] { (s0, e0) } else { (e0, s0) };
            let key = reference.pos[left];
            if let Some((beg, end)) = region {
                // The record lies within the boundary nodes' span, or starts one
                // base before it; exact overlap is checked on the record itself.
                let span_end = reference.pos[right] + g.sequence_len(reference.steps[right].0).unwrap();
                if key > end || span_end + 1 < beg {
                    return None;
                }
            }
            Some(RefTrav { snarl: i, trav, start_step: s0, end_step: e0, key })
        })
        .collect();
    let on_reference = refs.len();
    let depth = |mut i: usize| {
        let mut d = 0;
        while parent[i] != u32::MAX {
            i = parent[i] as usize;
            d += 1;
        }
        d
    };
    refs.par_sort_by_key(|r| (r.key, depth(r.snarl), r.snarl));

    // Heavy snarls, in position order (ancestors first), by reference
    // traversal and by snarl.
    let mut heavy_of = vec![u32::MAX; refs.len()];
    let mut heavy_of_snarl: HashMap<usize, usize> = HashMap::new();
    let mut heavy_refs = Vec::new();
    let mut on_ref = vec![false; snarls.len()];
    for (k, rt) in refs.iter().enumerate() {
        on_ref[rt.snarl] = true;
        if reference.pos[rt.start_step].abs_diff(reference.pos[rt.end_step]) >= HEAVY_SPAN {
            heavy_of[k] = heavy_refs.len() as u32;
            heavy_of_snarl.insert(rt.snarl, heavy_refs.len());
            heavy_refs.push(k);
        }
    }
    let heavy: Vec<Heavy> = heavy_refs.iter().map(|_| Heavy::new()).collect();
    let next_heavy = AtomicUsize::new(0);

    let header = vcf_header(&reference.chrom, reference.length);
    let path = format!("{out_prefix}.pgsites");
    let mut pgs = pgsites::Writer::create(&path, &header, opts.level, opts.threads).map_err(|e| format!("{path}: {e}"))?;
    let chrom = reference.chrom.as_str();
    let (level, threads) = (opts.level, opts.threads);
    // A giant record's block, compressed from its walks.
    let giant_block = |r: &Record, lv: usize, ps: Option<&str>, len: usize| -> std::io::Result<pgsites::Block> {
        let bytes = pgsites::compress_record(len, level, threads, |out| r.write_line(g, chrom, lv, ps, out))?;
        Ok(pgsites::Block { pos: r.pos, max_end: r.end(), min_allele: r.min_allele(), max_walk: r.max_walk(), raw_len: len + 1, bytes })
    };
    // A heavy snarl as a helper builds it: a giant record whose ancestors are
    // known is compressed here, and its walks freed.
    let build_heavy = |rt: &RefTrav| -> Built {
        let Some(r) = build_record(g, reference, &snarls[rt.snarl], rt, true) else { return Built::Record(None) };
        if region.is_some_and(|(beg, end)| r.pos > end || r.end() < beg) {
            return Built::Record(Some(r));
        }
        let ancestor = |pi: usize| match heavy_of_snarl.get(&pi) {
            Some(&h) => Some(heavy[h].emitted()),
            None if on_ref[pi] => None,
            None => Some(false),
        };
        match nesting(rt.snarl, &parent, &snarls, ancestor) {
            Some((lv, ps)) => {
                let len = r.line_len(chrom, lv, ps.as_deref());
                if len <= STREAM_BYTES {
                    return Built::Record(Some(r));
                }
                match giant_block(&r, lv, ps.as_deref(), len) {
                    Ok(b) => Built::Block(b),
                    Err(_) => Built::Record(Some(r)), // the writer reports it
                }
            }
            None => Built::Record(Some(r)),
        }
    };

    let mut emitted = vec![false; snarls.len()];
    let (mut records, mut lo) = (0usize, 0usize);
    let blocks = std::thread::scope(|scope| -> Result<usize, String> {
        for _ in 0..(threads / 2).clamp(1, heavy_refs.len().max(1)) {
            scope.spawn(|| loop {
                let h = next_heavy.fetch_add(1, Ordering::Relaxed);
                let Some(&k) = heavy_refs.get(h) else { break };
                heavy[h].try_build(|| build_heavy(&refs[k]));
            });
        }
        let (tx, rx) = std::sync::mpsc::sync_channel::<Vec<Item>>(2);
        // The writer: lines and blocks in order.
        let writer = scope.spawn(move || -> std::io::Result<usize> {
            for items in rx {
                for it in items {
                    match it {
                        Item::Line { pos, end, min_allele, max_walk, line } => pgs.add(pos, end, min_allele, max_walk, &line),
                        Item::Giant { record, lv, ps, len } => pgs.add_block(giant_block(&record, lv, ps.as_deref(), len)?)?,
                        Item::Block(b) => pgs.add_block(b)?,
                    }
                }
                pgs.flush()?;
            }
            pgs.finish()
        });
        while lo < refs.len() {
            // A batch ends where the next key is at least two bases on, so records
            // (POS is the key or one base later) never interleave across batches.
            let mut hi = (lo + BATCH).min(refs.len());
            while hi < refs.len() && refs[hi].key < refs[hi - 1].key + 2 {
                hi += 1;
            }
            // The batch's snarls in parallel, but its heavy ones collected here,
            // off the pool: a pool thread waiting for a heavy snarl could leave
            // no thread for that snarl's own parallel walk search.
            let mut built: Vec<(usize, Built)> = (lo..hi)
                .into_par_iter()
                .filter(|&k| heavy_of[k] == u32::MAX)
                .map(|k| {
                    let rt = &refs[k];
                    (rt.snarl, Built::Record(build_record(g, reference, &snarls[rt.snarl], rt, false)))
                })
                .collect();
            for k in (lo..hi).filter(|&k| heavy_of[k] != u32::MAX) {
                let rt = &refs[k];
                built.push((rt.snarl, heavy[heavy_of[k] as usize].take(|| build_heavy(rt))));
            }
            for (i, b) in &built {
                emitted[*i] = !matches!(b, Built::Record(None));
            }
            let emitted = &emitted;
            let mut items: Vec<(usize, String, Item)> = built
                .into_par_iter()
                .filter_map(|(i, b)| {
                    let r = match b {
                        Built::Block(blk) => {
                            let id = snarl_name(snarls[i].start, snarls[i].end);
                            return Some((blk.pos, format!("{id}\t"), Item::Block(blk)));
                        }
                        Built::Record(r) => r?,
                    };
                    if let Some((beg, end)) = region {
                        if r.pos > end || r.end() < beg {
                            return None;
                        }
                    }
                    let (lv, ps) = nesting(i, &parent, &snarls, |pi| Some(emitted[pi])).unwrap();
                    // Order within a position is the lines' byte order, which the
                    // ID (unique, then a tab) decides.
                    let key = format!("{}\t", r.id);
                    let len = r.line_len(chrom, lv, ps.as_deref());
                    let item = if len > STREAM_BYTES {
                        Item::Giant { record: r, lv, ps, len }
                    } else {
                        let mut line = Vec::with_capacity(len);
                        r.write_line(g, chrom, lv, ps.as_deref(), &mut line).unwrap();
                        assert_eq!(line.len(), len, "record {} length", r.id);
                        Item::Line { pos: r.pos, end: r.end(), min_allele: r.min_allele(), max_walk: r.max_walk(), line }
                    };
                    Some((item.pos(), key, item))
                })
                .collect();
            items.par_sort_unstable_by(|a, b| (a.0, &a.1).cmp(&(b.0, &b.1)));
            records += items.len();
            if tx.send(items.into_iter().map(|x| x.2).collect()).is_err() {
                break; // the writer failed; its error is reported below
            }
            lo = hi;
        }
        drop(tx);
        writer.join().map_err(|_| "writer thread panicked".to_string())?.map_err(|e| format!("{path}: {e}"))
    })?;
    Ok(Summary {
        snarls: snarls.len(),
        on_reference,
        records,
        blocks,
        seconds: t0.elapsed().as_secs_f64(),
    })
}
