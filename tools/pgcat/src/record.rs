//! One snarl's record: its haplotype walks, alleles and VCF line, written the
//! way `vg deconstruct -a` writes them, without what needs haplotype names
//! (genotypes, AC/AF/AN/NS). ALTs are ordered by how many haplotypes take them.
//!
//! An allele is its walk: its bases are the walk's inner node sequences. A
//! record keeps the walks and writes the line from them, so an allele's bases
//! exist only while it is written (or hashed, to find equal alleles). A snarl
//! spanning megabases with hundreds of haplotype alleles (chr2 of HPRC v2.1:
//! 11.4 Mb, 435 alleles, a 6.5 GB line) costs its walks, not its alleles.

use crate::reference::ReferencePath;
use crate::snarls::{snarl_name, Handle, Snarl};
use gbz::support::{self, Orientation};
use gbz::GBZ;
use std::collections::HashMap;
use std::hash::Hasher;
use std::io::Write;

/// A snarl's single traversal by the reference path.
pub struct RefTrav {
    pub snarl: usize,
    pub trav: Vec<Handle>,
    pub start_step: usize,
    pub end_step: usize,
    /// 0-based position of the snarl's leftmost boundary node: the sort key.
    pub key: usize,
}

/// A node visit as one integer: node id * 2, plus 1 if reverse.
type Visit = u64;

fn visit(h: Handle) -> Visit {
    ((h.0 as u64) << 1) | (h.1 == Orientation::Reverse) as u64
}

fn node(v: Visit) -> usize {
    (v >> 1) as usize
}

fn reverse(v: Visit) -> bool {
    v & 1 == 1
}

fn digits(mut n: usize) -> usize {
    let mut d = 1;
    while n >= 10 {
        n /= 10;
        d += 1;
    }
    d
}

/// A node's bases in the visit's orientation, upper case.
fn oriented_seq(g: &GBZ, v: Visit) -> Vec<u8> {
    let s = g.sequence(node(v)).unwrap();
    let s = if reverse(v) { support::reverse_complement(s) } else { s.to_vec() };
    s.to_ascii_uppercase()
}

/// A walk's nodes between its two boundary nodes.
fn inner(w: &[Visit]) -> &[Visit] {
    if w.len() < 2 { &[] } else { &w[1..w.len() - 1] }
}

/// The bases between a walk's boundary nodes, in the snarl's orientation.
fn walk_string(g: &GBZ, w: &[Visit]) -> Vec<u8> {
    let inner = inner(w);
    let mut s = Vec::with_capacity(inner.iter().map(|&v| g.sequence_len(node(v)).unwrap()).sum());
    for &v in inner {
        s.extend(oriented_seq(g, v));
    }
    s
}

/// A walk's bases' length and hash (that of walk_string, without building it:
/// the hasher takes the bases node by node).
fn walk_key(g: &GBZ, w: &[Visit]) -> (usize, u64) {
    let mut h = std::collections::hash_map::DefaultHasher::new();
    let mut len = 0;
    for &v in inner(w) {
        let s = oriented_seq(g, v);
        len += s.len();
        h.write(&s);
    }
    (len, h.finish())
}

/// Snarls with more walk steps than this hash their alleles in parallel.
const PARALLEL_STEPS: usize = 1 << 20;

pub struct Record {
    pub pos: usize,
    pub id: String,
    /// The reference's walk, then one walk per distinct ALT sequence.
    walks: Vec<Vec<Visit>>,
    /// Each allele's length as written, padding included.
    lens: Vec<usize>,
    /// The base written before every allele when they are not all one length.
    pad: Option<u8>,
    /// The reference reaches the snarl's end node first: alleles are written
    /// reverse-complemented and traversals backwards.
    reversed: bool,
}

/// vg's list_haplotypes: the distinct haplotype walks from `start` to `end`,
/// with how many haplotypes take each, by forward GBWT search. Walks share
/// prefixes in an arena and are materialised only on reaching `end`. A walk
/// that steps onto a boundary node backwards has left the snarl and is dropped:
/// it would otherwise follow its haplotype to the end of the chromosome.
fn haplotype_walks(g: &GBZ, start: Handle, end: Handle) -> Vec<(Vec<Visit>, usize)> {
    let mut out = Vec::new();
    let Some(first) = g.search_state(start.0, start.1) else { return out };
    let (start, end) = (visit(start), visit(end));
    let (start_back, end_back) = (start ^ 1, end ^ 1);
    let mut arena: Vec<(Visit, u32)> = vec![(start, u32::MAX)];
    let mut stack = vec![(0u32, 1usize, first)];
    while let Some((at, depth, state)) = stack.pop() {
        let Some(it) = g.follow_forward(&state) else { continue };
        for next in it {
            let v = visit(support::decode_node(next.forward.node));
            if v == end_back || v == start_back {
                continue;
            }
            arena.push((v, at));
            let idx = (arena.len() - 1) as u32;
            if v == end {
                let mut w = Vec::with_capacity(depth + 1);
                let mut k = idx;
                while k != u32::MAX {
                    w.push(arena[k as usize].0);
                    k = arena[k as usize].1;
                }
                w.reverse();
                out.push((w, next.len()));
            } else {
                stack.push((idx, depth + 1, next));
            }
        }
    }
    out
}

/// The visits a parallel walk search had made when it split off a task: the
/// parent's prefix, then these.
struct Prefix {
    parent: Option<std::sync::Arc<Prefix>>,
    steps: Vec<Visit>,
}

impl Prefix {
    fn extend_into(&self, out: &mut Vec<Visit>) {
        if let Some(p) = &self.parent {
            p.extend_into(out);
        }
        out.extend_from_slice(&self.steps);
    }
}

/// Tasks a parallel walk search may split into.
const SEARCH_TASKS: usize = 512;

/// haplotype_walks, its search tree split across rayon tasks at branches (the
/// same walks and counts, in another order).
fn haplotype_walks_parallel(g: &GBZ, start: Handle, end: Handle) -> Vec<(Vec<Visit>, usize)> {
    use std::sync::atomic::{AtomicUsize, Ordering};
    use std::sync::{Arc, Mutex};
    let Some(first) = g.search_state(start.0, start.1) else { return Vec::new() };
    let (start, end) = (visit(start), visit(end));
    let out = Mutex::new(Vec::new());
    let tasks = AtomicUsize::new(1);

    // One task: a serial search below `prefix`, whose last visit `state` is at.
    fn task<'s>(
        scope: &rayon::Scope<'s>,
        g: &'s GBZ,
        prefix: Arc<Prefix>,
        state: gbz::BidirectionalState,
        (start, end): (Visit, Visit),
        out: &'s Mutex<Vec<(Vec<Visit>, usize)>>,
        tasks: &'s AtomicUsize,
    ) {
        let (start_back, end_back) = (start ^ 1, end ^ 1);
        // Local arena: (visit, parent index), u32::MAX for the prefix.
        let mut arena: Vec<(Visit, u32)> = Vec::new();
        let mut stack = vec![(u32::MAX, state)];
        let local_path = |arena: &[(Visit, u32)], mut k: u32| {
            let mut w = Vec::new();
            while k != u32::MAX {
                w.push(arena[k as usize].0);
                k = arena[k as usize].1;
            }
            w.reverse();
            w
        };
        while let Some((at, state)) = stack.pop() {
            let Some(it) = g.follow_forward(&state) else { continue };
            let nexts: Vec<_> = it
                .filter(|n| {
                    let v = visit(support::decode_node(n.forward.node));
                    v != end_back && v != start_back
                })
                .collect();
            let branching = nexts.len() > 1;
            for (j, next) in nexts.into_iter().enumerate() {
                let v = visit(support::decode_node(next.forward.node));
                if v == end {
                    let mut w = Vec::new();
                    prefix.extend_into(&mut w);
                    w.extend(local_path(&arena, at));
                    w.push(v);
                    out.lock().unwrap().push((w, next.len()));
                } else if branching && j > 0 && tasks.fetch_add(1, Ordering::Relaxed) < SEARCH_TASKS {
                    let mut steps = local_path(&arena, at);
                    steps.push(v);
                    let p = Arc::new(Prefix { parent: Some(prefix.clone()), steps });
                    scope.spawn(move |s| task(s, g, p, next, (start, end), out, tasks));
                } else {
                    arena.push((v, at));
                    stack.push(((arena.len() - 1) as u32, next));
                }
            }
        }
    }

    let root = Arc::new(Prefix { parent: None, steps: vec![start] });
    rayon::scope(|s| task(s, g, root, first, (start, end), &out, &tasks));
    out.into_inner().unwrap()
}

/// The record of a snarl the reference passes through, or None when every
/// haplotype spells the reference. `parallel`: search its walks in parallel
/// (for snarls spanning tens of kilobases or more).
pub fn build_record(g: &GBZ, reference: &ReferencePath, sn: &Snarl, rt: &RefTrav, parallel: bool) -> Option<Record> {
    let ref_walk: Vec<Visit> = rt.trav.iter().map(|&h| visit(h)).collect();
    let mut walks = if parallel { haplotype_walks_parallel(g, sn.start, sn.end) } else { haplotype_walks(g, sn.start, sn.end) };
    let start_pos = reference.pos[rt.start_step];
    let end_pos = reference.pos[rt.end_step];
    // As vg: VCF position and padding base come from the boundary node the
    // reference reaches first.
    let use_start = start_pos < end_pos;
    let first_path = visit(reference.steps[if use_start { rt.start_step } else { rt.end_step }]);
    let mut first_path_pos = if use_start { start_pos } else { end_pos };
    let first_trav = if use_start { ref_walk[0] } else { *ref_walk.last().unwrap() };
    let same = (reverse(first_trav) == reverse(first_path)) == use_start;
    let fp_seq = oriented_seq(g, first_path);
    let prev_char = if same {
        first_path_pos += fp_seq.len();
        fp_seq[fp_seq.len() - 1]
    } else {
        support::reverse_complement(&fp_seq[..1])[0]
    };
    let mut pos = first_path_pos + 1;

    // Alleles: the reference, then distinct sequences, most haplotypes first.
    // Equal sequences are found by hash and confirmed byte by byte.
    walks.retain(|(w, _)| *w != ref_walk);
    walks.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
    let (ref_len, ref_hash) = walk_key(g, &ref_walk);
    let keys: Vec<(usize, u64)> = if walks.iter().map(|(w, _)| w.len()).sum::<usize>() > PARALLEL_STEPS {
        use rayon::prelude::*;
        walks.par_iter().map(|(w, _)| walk_key(g, w)).collect()
    } else {
        walks.iter().map(|(w, _)| walk_key(g, w)).collect()
    };
    let mut by_hash: HashMap<u64, Vec<usize>> = HashMap::new();
    by_hash.entry(ref_hash).or_default().push(0);
    let mut kept: Vec<Vec<Visit>> = vec![ref_walk];
    let mut lens: Vec<usize> = vec![ref_len];
    let mut substitution = true;
    for ((w, _), (len, h)) in walks.into_iter().zip(keys) {
        // Equal length and hash are confirmed byte by byte.
        let seen = by_hash.get(&h).is_some_and(|ks| {
            ks.iter().any(|&k| lens[k] == len && walk_string(g, &kept[k]) == walk_string(g, &w))
        });
        if seen {
            continue;
        }
        substitution = substitution && len == ref_len;
        by_hash.entry(h).or_default().push(kept.len());
        lens.push(len);
        kept.push(w);
    }
    if kept.len() < 2 {
        return None;
    }
    let pad = if substitution { None } else { Some(prev_char) };
    if pad.is_some() {
        pos -= 1;
        for l in &mut lens {
            *l += 1;
        }
    }
    Some(Record { pos, id: snarl_name(sn.start, sn.end), walks: kept, lens, pad, reversed: !use_start })
}

impl Record {
    /// Last reference base the record covers, as tabix computes it for VCF.
    pub fn end(&self) -> usize {
        self.pos + self.lens[0].max(1) - 1
    }

    /// The shortest allele's bases, the padding base not counted: no read
    /// aligned along fewer graph bases can traverse any allele whole.
    pub fn min_allele(&self) -> usize {
        self.lens.iter().min().unwrap() - self.pad.is_some() as usize
    }

    /// The most nodes any allele's walk visits, boundary nodes included.
    pub fn max_walk(&self) -> usize {
        self.walks.iter().map(|w| w.len()).max().unwrap()
    }

    /// The length of the line `write_line` writes, newline excluded.
    pub fn line_len(&self, chrom: &str, lv: usize, ps: Option<&str>) -> usize {
        let at: usize = self.walks.iter().map(|w| w.iter().map(|&v| 1 + digits(node(v))).sum::<usize>()).sum();
        chrom.len() + 1 + digits(self.pos) + 1 + self.id.len() + 1
            + self.lens.iter().sum::<usize>() + (self.lens.len() - 2) + "\t\t60\t.\tAT=".len()
            + at + (self.walks.len() - 1)
            + ";LV=".len() + digits(lv) + ps.map(|p| ";PS=".len() + p.len()).unwrap_or(0)
    }

    fn write_allele<W: Write + ?Sized>(&self, g: &GBZ, k: usize, out: &mut W) -> std::io::Result<()> {
        if let Some(p) = self.pad {
            out.write_all(&[p])?;
        }
        let inner = inner(&self.walks[k]);
        if self.reversed {
            // The bases reverse-complemented: the inner nodes backwards, each
            // written in the opposite orientation.
            for &v in inner.iter().rev() {
                out.write_all(&oriented_seq(g, v ^ 1))?;
            }
        } else {
            for &v in inner {
                out.write_all(&oriented_seq(g, v))?;
            }
        }
        Ok(())
    }

    fn write_at<W: Write + ?Sized>(&self, k: usize, out: &mut W) -> std::io::Result<()> {
        let w = &self.walks[k];
        let n = w.len();
        for i in 0..n {
            let v = if self.reversed { w[n - 1 - i] } else { w[i] };
            out.write_all(if reverse(v) != self.reversed { b"<" } else { b">" })?;
            write!(out, "{}", node(v))?;
        }
        Ok(())
    }

    /// The VCF line, without its newline.
    pub fn write_line<W: Write + ?Sized>(&self, g: &GBZ, chrom: &str, lv: usize, ps: Option<&str>, out: &mut W) -> std::io::Result<()> {
        write!(out, "{chrom}\t{}\t{}\t", self.pos, self.id)?;
        self.write_allele(g, 0, out)?;
        out.write_all(b"\t")?;
        for k in 1..self.walks.len() {
            if k > 1 {
                out.write_all(b",")?;
            }
            self.write_allele(g, k, out)?;
        }
        out.write_all(b"\t60\t.\tAT=")?;
        for k in 0..self.walks.len() {
            if k > 0 {
                out.write_all(b",")?;
            }
            self.write_at(k, out)?;
        }
        write!(out, ";LV={lv}")?;
        if let Some(p) = ps {
            write!(out, ";PS={p}")?;
        }
        Ok(())
    }
}
