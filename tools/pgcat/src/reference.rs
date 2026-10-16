//! A contig's reference path: its steps, their positions, and the steps of each
//! node, built from every fragment of `<sample>#0#<contig>` in the GBZ.
//!
//! A subgraph GBZ stores a reference contig as one or more fragments, each with
//! the offset of its first base (`CHM13#0#chr6[31350872]`). The fragments are
//! laid end to end in offset order; positions include the offsets, and a
//! traversal never runs from one fragment into the next.

use crate::snarls::{flip, Handle};
use gbz::support::Orientation;
use gbz::GBZ;
use std::collections::{HashMap, HashSet};

pub struct ReferencePath {
    /// `<sample>#0#<contig>`: the CHROM written in the catalog.
    pub chrom: String,
    pub steps: Vec<Handle>,
    /// 0-based position of each step's first base on the contig.
    pub pos: Vec<usize>,
    /// For each step, one past the last step of its fragment.
    frag_end: Vec<u32>,
    /// For each step, the first step of its fragment.
    frag_beg: Vec<u32>,
    steps_of: HashMap<usize, Vec<u32>>,
    /// One past the contig's last reference base.
    pub length: usize,
}

impl ReferencePath {
    pub fn load(g: &GBZ, sample: &str, contig: &str) -> Result<Self, String> {
        let md = g.metadata().ok_or("the GBZ has no path metadata")?;
        let missing = || format!("no reference path {sample}#0#{contig} in the graph");
        let sample_id = md.sample_id(sample).ok_or_else(missing)?;
        let contig_id = md.contig_id(contig).ok_or_else(missing)?;
        let mut fragments: Vec<(usize, usize)> = Vec::new(); // (offset, path id)
        for id in 0..md.paths() {
            let Some(name) = md.path(id) else { continue };
            if name.sample() == sample_id && name.contig() == contig_id && name.phase() == 0 {
                fragments.push((name.fragment as usize, id));
            }
        }
        if fragments.is_empty() {
            return Err(missing());
        }
        fragments.sort();
        let mut r = ReferencePath {
            chrom: format!("{sample}#0#{contig}"),
            steps: Vec::new(),
            pos: Vec::new(),
            frag_end: Vec::new(),
            frag_beg: Vec::new(),
            steps_of: HashMap::new(),
            length: 0,
        };
        for (offset, id) in fragments {
            let beg = r.steps.len() as u32;
            let mut at = offset;
            for h in g.path(id, Orientation::Forward).ok_or("unreadable reference path")? {
                r.pos.push(at);
                at += g.sequence_len(h.0).unwrap();
                r.steps.push(h);
                r.frag_beg.push(beg);
            }
            let end = r.steps.len() as u32;
            r.frag_end.resize(r.steps.len(), end);
            r.length = r.length.max(at);
        }
        for (k, h) in r.steps.iter().enumerate() {
            r.steps_of.entry(h.0).or_default().push(k as u32);
        }
        Ok(r)
    }

    /// vg's PathTraversalFinder on this path: every traversal from `start` to
    /// `end`, as (traversal in snarl orientation, start step, end step). From
    /// each visit of the start node the path is walked forward, and if that does
    /// not reach `end`, backward; it may not step onto a node just outside the
    /// snarl, or past the end of its fragment.
    pub fn traversals(&self, g: &GBZ, start: Handle, end: Handle) -> Vec<(Vec<Handle>, usize, usize)> {
        let mut out = Vec::new();
        let Some(starts) = self.steps_of.get(&start.0) else { return out };
        if !self.steps_of.contains_key(&end.0) {
            return out;
        }
        let mut block: HashSet<usize> = HashSet::new();
        if let Some(it) = g.predecessors(start.0, start.1) {
            block.extend(it.map(|(n, _)| n));
        }
        if let Some(it) = g.successors(end.0, end.1) {
            block.extend(it.map(|(n, _)| n));
        }
        for &s0 in starts {
            let s0 = s0 as usize;
            let (fbeg, fend) = (self.frag_beg[s0] as usize, self.frag_end[s0] as usize);
            let mut trav = Vec::new();
            let mut k = s0;
            loop {
                let h = self.steps[k];
                trav.push(h);
                if h == end || k + 1 >= fend {
                    break;
                }
                let next = self.steps[k + 1];
                if next == end || !block.contains(&next.0) {
                    k += 1;
                } else {
                    break;
                }
            }
            if self.steps[k] == end {
                out.push((trav, s0, k));
                continue;
            }
            let end_check = flip(end);
            trav.clear();
            k = s0;
            loop {
                let h = flip(self.steps[k]);
                trav.push(h);
                if h == end || k == fbeg {
                    break;
                }
                let prev = self.steps[k - 1];
                if prev == end_check || !block.contains(&prev.0) {
                    k -= 1;
                } else {
                    break;
                }
            }
            if self.steps[k] == end_check {
                out.push((trav, s0, k));
            }
        }
        out
    }
}
