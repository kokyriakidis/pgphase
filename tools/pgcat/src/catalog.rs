//! Building one contig's catalog, or one region of it, from its snarls.
//!
//! Pass 1 finds each snarl's reference traversal; snarls the reference does not
//! pass through have no position and no record (vg deconstruct skips them
//! too). The rest are sorted by position, ancestors first at a tie, and built
//! in batches: in parallel, then LV/PS from the ancestors already built, then
//! written in POS order by a writer thread while the next batch builds. Only one
//! or two batches are ever in memory.

use crate::pgsites;
use crate::record::{build_record, RefTrav, Record};
use crate::reference::ReferencePath;
use crate::snarls::{canonical, read_snarls, snarl_name};
use gbz::GBZ;
use rayon::prelude::*;
use std::collections::HashMap;
use std::time::Instant;

/// Snarls built per batch: enough to keep every thread busy, few enough that
/// one or two batches of records fit in memory.
const BATCH: usize = 50_000;

/// Records with longer lines are not formatted in memory: the writer streams
/// them from their walks, one at a time, into a block of their own.
const STREAM_BYTES: usize = 64 << 20;

/// What the writer receives, in file order.
enum Item {
    Line { pos: usize, end: usize, min_allele: usize, max_walk: usize, line: Vec<u8> },
    Giant { record: Record, lv: usize, ps: Option<String>, len: usize },
}

impl Item {
    fn pos(&self) -> usize {
        match self {
            Item::Line { pos, .. } => *pos,
            Item::Giant { record, .. } => record.pos,
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

    let header = vcf_header(&reference.chrom, reference.length);
    let path = format!("{out_prefix}.pgsites");
    let mut pgs = pgsites::Writer::create(&path, &header, opts.level, opts.threads).map_err(|e| format!("{path}: {e}"))?;
    let chrom = reference.chrom.as_str();
    let mut emitted = vec![false; snarls.len()];
    let (mut records, mut lo) = (0usize, 0usize);
    let blocks = std::thread::scope(|scope| -> Result<usize, String> {
        let (tx, rx) = std::sync::mpsc::sync_channel::<Vec<Item>>(2);
        // The writer: lines in order, a giant record streamed from its walks.
        let writer = scope.spawn(move || -> std::io::Result<usize> {
            for items in rx {
                for it in items {
                    match it {
                        Item::Line { pos, end, min_allele, max_walk, line } => pgs.add(pos, end, min_allele, max_walk, &line),
                        Item::Giant { record, lv, ps, len } => {
                            let ps = ps.as_deref();
                            pgs.add_streamed(record.pos, record.end(), record.min_allele(), record.max_walk(), len, |out| {
                                record.write_line(g, chrom, lv, ps, out)
                            })?
                        }
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
            let built: Vec<(usize, Option<Record>)> = refs[lo..hi]
                .par_iter()
                .map(|rt| (rt.snarl, build_record(g, reference, &snarls[rt.snarl], rt)))
                .collect();
            for (i, r) in &built {
                emitted[*i] = r.is_some();
            }
            let emitted = &emitted;
            let mut items: Vec<(usize, String, Item)> = built
                .into_par_iter()
                .filter_map(|(i, r)| {
                    let r = r?;
                    if let Some((beg, end)) = region {
                        if r.pos > end || r.end() < beg {
                            return None;
                        }
                    }
                    // LV / PS: ancestors that have a record, as vg's update_nesting_info_tags.
                    let (mut lv, mut ps) = (0usize, None);
                    let mut cur = parent[i];
                    while cur != u32::MAX {
                        let pi = cur as usize;
                        if emitted[pi] {
                            lv += 1;
                            if ps.is_none() {
                                ps = Some(snarl_name(snarls[pi].start, snarls[pi].end));
                            }
                        }
                        cur = parent[pi];
                    }
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
