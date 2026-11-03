//! `pgcat index-gaf`: a GAF made searchable by reference position, as pgphase
//! reads it. Each alignment is placed on the reference by the nodes of its walk
//! that lie on the reference sample's paths (pggaf's rule, below), written as
//! `contig  start+1  end  <the GAF line>`, sorted, bgzipped and tabix-indexed.
//!
//! Alignments are projected in parallel batches as the input is read; sorted
//! runs that exceed the memory budget are spilled, one zstd frame per BUCKET of
//! the reference, and merged bucket by bucket: each bucket's frames from every
//! run are sorted on a thread of their own, many buckets at once. Header
//! lines `#max_span <contig> <bases>` give the most graph bases any alignment
//! on each contig spans.

use crate::bgzf::IndexedWriter;
use gbz::support::Orientation;
use gbz::GBZ;
use rayon::prelude::*;
use std::collections::{BinaryHeap, HashMap};
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

pub struct IndexGafOptions {
    pub graph: String,
    pub gaf: String,
    pub out: String,
    pub reference_sample: Option<String>,
    /// Bytes of alignments sorted in memory before a run is spilled.
    pub memory: usize,
    pub tmp_dir: Option<String>,
    /// BGZF compression level.
    pub level: u32,
    /// Where the alignments with no reference interval are written (gzip GAF, each line tagged `pu:Z:REASON`),
    /// rather than dropped. One without anchors is also tagged with the reference beside it (RefMap::nearest):
    /// `pb:Z:CONTIG:OFFSET:BASES` before its walk and `pa:Z:...` after, OFFSET the 0-based reference position where
    /// the graph search joins the reference and BASES the graph bases between.
    pub unplaced: Option<String>,
    /// Align each alignment without anchors to the reference near the nodes beside it (rescue::rescue) and place it
    /// there, tagged rs:Z:CONTIG:POS:STRAND:CIGAR and rq:i:MAPQ, rather than leave it unplaced.
    pub rescue: bool,
    /// Write each alignment's 12 GAF columns and SLIM_TAGS alone: what pgphase reads (cs:Z:, bq:Z: and the like, most
    /// of a long read's line, it never does).
    pub slim: bool,
    /// Place alignments by this ref map (`pgcat ref-map`) instead of the graph's paths, the graph not loaded: for an
    /// index without --unplaced or --rescue, which search the graph.
    pub ref_map: Option<String>,
}

/// `pgcat ref-map`: the reference map of a graph's reference sample, saved for index-gaf --ref-map.
pub fn save_ref_map(graph: &str, sample: Option<&str>, out: &str) -> Result<(), String> {
    let g: GBZ = simple_sds::serialize::load_from(graph).map_err(|e| format!("{graph}: {e}"))?;
    let sample = crate::graph::resolve_reference_sample(graph, sample)?;
    RefMap::build(&g, &sample)?.save(&sample, out)
}

/// The tags a slim index keeps: rq:i: (the rescue's MAPQ, which pgphase reads), pu:Z: (why a line had no place) and
/// rs:Z: (where the rescue put it).
const SLIM_TAGS: [&[u8]; 3] = [b"rq:i:", b"pu:Z:", b"rs:Z:"];

/// Appends `line` (no newline) to `out`: whole, or with `slim` its 12 columns and SLIM_TAGS alone.
fn append_line(out: &mut Vec<u8>, line: &[u8], slim: bool) {
    if !slim {
        out.extend_from_slice(line);
        return;
    }
    let mut fields = line.split(|&c| c == b'\t');
    for (k, f) in fields.by_ref().take(12).enumerate() {
        if k > 0 {
            out.push(b'\t');
        }
        out.extend_from_slice(f);
    }
    for f in fields {
        if SLIM_TAGS.iter().any(|t| f.starts_with(t)) {
            out.push(b'\t');
            out.extend_from_slice(f);
        }
    }
}

/// Why an alignment with a walk has no reference interval.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Unplaced {
    /// No node of the aligned part of the walk lies on the reference.
    NoAnchor,
    /// A node of the walk is not in the graph.
    MissingNode,
    /// The anchors fit the reference at more than one place.
    Ambiguous,
    /// The anchors fit no single path in one orientation and order (an inversion, a rearrangement).
    Inconsistent,
}

impl Unplaced {
    const ALL: [Unplaced; 4] = [Unplaced::NoAnchor, Unplaced::MissingNode, Unplaced::Ambiguous, Unplaced::Inconsistent];

    pub fn name(self) -> &'static str {
        match self {
            Unplaced::NoAnchor => "no_anchor",
            Unplaced::MissingNode => "missing_node",
            Unplaced::Ambiguous => "ambiguous",
            Unplaced::Inconsistent => "inconsistent",
        }
    }
}

/// A node's place on a reference path.
#[derive(Clone, Copy)]
struct Placement {
    path: u32,
    step: u32,
    start: u64,
    reverse: bool,
}

const NONE: u32 = u32::MAX;

/// Node lengths and the reference placements of every node, by node id.
/// A cursor over a byte slice, handing out its next bytes.
struct Cursor<'a> {
    data: &'a [u8],
    at: usize,
}

impl<'a> Cursor<'a> {
    fn take(&mut self, n: usize) -> Option<&'a [u8]> {
        let s = self.data.get(self.at..self.at.checked_add(n)?)?;
        self.at += n;
        Some(s)
    }
}

struct RefMap {
    node_len: Vec<u32>,
    first: Vec<u32>,
    placements: Vec<Placement>,
    /// Placements after a node's first, for nodes the reference visits twice.
    more: HashMap<usize, Vec<u32>>,
    path_contig: Vec<u32>,
    contigs: Vec<String>,
}

impl RefMap {
    fn build(g: &GBZ, sample: &str) -> Result<Self, String> {
        let md = g.metadata().ok_or("the GBZ has no path metadata")?;
        let sample_id = md.sample_id(sample).ok_or_else(|| format!("no paths of sample {sample} in the graph"))?;
        let size = g.max_node() + 1;
        let mut node_len = vec![NONE; size];
        for id in g.node_iter() {
            node_len[id] = g.sequence_len(id).unwrap() as u32;
        }
        let mut m = RefMap {
            node_len,
            first: vec![NONE; size],
            placements: Vec::new(),
            more: HashMap::new(),
            path_contig: Vec::new(),
            contigs: Vec::new(),
        };
        let mut contig_ids: HashMap<String, u32> = HashMap::new();
        for id in 0..md.paths() {
            let Some(p) = md.path(id) else { continue };
            if p.sample() != sample_id {
                continue;
            }
            let name = md.contig_name(p.contig());
            let next = m.contigs.len() as u32;
            let contig = *contig_ids.entry(name.clone()).or_insert_with(|| {
                m.contigs.push(name);
                next
            });
            let path = m.path_contig.len() as u32;
            m.path_contig.push(contig);
            let mut pos = p.fragment as u64;
            for (step, (node, o)) in g.path(id, Orientation::Forward).ok_or("unreadable reference path")?.enumerate() {
                let k = m.placements.len() as u32;
                m.placements.push(Placement { path, step: step as u32, start: pos, reverse: o == Orientation::Reverse });
                if m.first[node] == NONE {
                    m.first[node] = k;
                } else {
                    m.more.entry(node).or_default().push(k);
                }
                pos += m.node_len[node] as u64;
            }
        }
        if m.path_contig.is_empty() {
            return Err(format!("no paths of sample {sample} in the graph"));
        }
        Ok(m)
    }

    /// Writes the map to `path` (a ref map: `pgcat ref-map`), for index-gaf to load instead of the graph. Little-endian:
    /// "PGCATRM1"; the reference sample; node_len and first (u32 each, per node id); the placements (path, step,
    /// start, reverse); the further placements of nodes visited twice (node, then its list); path_contig; contigs.
    fn save(&self, sample: &str, path: &str) -> Result<(), String> {
        use std::io::Write;
        let io = |e: std::io::Error| format!("{path}: {e}");
        let mut w = BufWriter::with_capacity(1 << 22, std::fs::File::create(path).map_err(io)?);
        let mut put = |b: &[u8]| w.write_all(b).map_err(io);
        let text = |s: &str| [(s.len() as u64).to_le_bytes().to_vec(), s.as_bytes().to_vec()].concat();
        put(b"PGCATRM1")?;
        put(&text(sample))?;
        put(&(self.node_len.len() as u64).to_le_bytes())?;
        for v in [&self.node_len, &self.first] {
            let mut buf = Vec::with_capacity(4 * v.len());
            for x in v.iter() {
                buf.extend_from_slice(&x.to_le_bytes());
            }
            put(&buf)?;
        }
        put(&(self.placements.len() as u64).to_le_bytes())?;
        let mut buf = Vec::with_capacity(17 * self.placements.len());
        for p in &self.placements {
            buf.extend_from_slice(&p.path.to_le_bytes());
            buf.extend_from_slice(&p.step.to_le_bytes());
            buf.extend_from_slice(&p.start.to_le_bytes());
            buf.push(p.reverse as u8);
        }
        put(&buf)?;
        let mut more: Vec<(&usize, &Vec<u32>)> = self.more.iter().collect();
        more.sort_unstable();
        put(&(more.len() as u64).to_le_bytes())?;
        for (node, ks) in more {
            put(&(*node as u64).to_le_bytes())?;
            put(&(ks.len() as u64).to_le_bytes())?;
            for k in ks {
                put(&k.to_le_bytes())?;
            }
        }
        put(&(self.path_contig.len() as u64).to_le_bytes())?;
        for c in &self.path_contig {
            put(&c.to_le_bytes())?;
        }
        put(&(self.contigs.len() as u64).to_le_bytes())?;
        for c in &self.contigs {
            put(&text(c))?;
        }
        drop(put);
        w.flush().map_err(io)
    }

    /// The map `save` wrote, and the reference sample it was made for.
    fn load(path: &str) -> Result<(Self, String), String> {
        let data = std::fs::read(path).map_err(|e| format!("{path}: {e}"))?;
        let bad = || format!("{path}: not a pgcat ref map (pgcat ref-map writes one)");
        let mut c = Cursor { data: &data, at: 0 };
        let u64_ = |b: &[u8]| u64::from_le_bytes(b.try_into().unwrap());
        let u32_ = |b: &[u8]| u32::from_le_bytes(b.try_into().unwrap());
        if c.take(8).ok_or_else(bad)? != b"PGCATRM1" {
            return Err(bad());
        }
        let n = u64_(c.take(8).ok_or_else(bad)?) as usize;
        let sample = String::from_utf8(c.take(n).ok_or_else(bad)?.to_vec()).map_err(|_| bad())?;
        let size = u64_(c.take(8).ok_or_else(bad)?) as usize;
        let node_len: Vec<u32> = c.take(4 * size).ok_or_else(bad)?.chunks_exact(4).map(u32_).collect();
        let first: Vec<u32> = c.take(4 * size).ok_or_else(bad)?.chunks_exact(4).map(u32_).collect();
        let count = u64_(c.take(8).ok_or_else(bad)?) as usize;
        let placements: Vec<Placement> = c
            .take(17 * count)
            .ok_or_else(bad)?
            .chunks_exact(17)
            .map(|b| Placement { path: u32_(&b[0..4]), step: u32_(&b[4..8]), start: u64_(&b[8..16]), reverse: b[16] != 0 })
            .collect();
        let mut more: HashMap<usize, Vec<u32>> = HashMap::new();
        for _ in 0..u64_(c.take(8).ok_or_else(bad)?) {
            let node = u64_(c.take(8).ok_or_else(bad)?) as usize;
            let k = u64_(c.take(8).ok_or_else(bad)?) as usize;
            more.insert(node, c.take(4 * k).ok_or_else(bad)?.chunks_exact(4).map(u32_).collect());
        }
        let k = u64_(c.take(8).ok_or_else(bad)?) as usize;
        let path_contig: Vec<u32> = c.take(4 * k).ok_or_else(bad)?.chunks_exact(4).map(u32_).collect();
        let mut contigs = Vec::new();
        for _ in 0..u64_(c.take(8).ok_or_else(bad)?) {
            let n = u64_(c.take(8).ok_or_else(bad)?) as usize;
            contigs.push(String::from_utf8(c.take(n).ok_or_else(bad)?.to_vec()).map_err(|_| bad())?);
        }
        Ok((RefMap { node_len, first, placements, more, path_contig, contigs }, sample))
    }

    /// Each contig's letters, from its reference paths (N where no path covers it).
    fn letters(&self, g: &GBZ) -> Vec<Vec<u8>> {
        let mut out: Vec<Vec<u8>> = vec![Vec::new(); self.contigs.len()];
        for node in 0..self.first.len() {
            for p in self.placements_of(node) {
                let (Some(s), seq) = (g.sequence(node), &mut out[self.path_contig[p.path as usize] as usize]) else { continue };
                let end = p.start as usize + s.len();
                if seq.len() < end {
                    seq.resize(end, b'N');
                }
                let at = &mut seq[p.start as usize..end];
                if p.reverse {
                    for (x, &c) in at.iter_mut().zip(s.iter().rev()) {
                        *x = crate::surject::comp(c);
                    }
                } else {
                    for (x, &c) in at.iter_mut().zip(s) {
                        *x = c.to_ascii_uppercase();
                    }
                }
            }
        }
        out
    }

    fn placements_of(&self, node: usize) -> impl Iterator<Item = &Placement> {
        let first = self.first.get(node).copied().unwrap_or(NONE);
        let more = if first == NONE { None } else { self.more.get(&node) };
        (first != NONE)
            .then(|| &self.placements[first as usize])
            .into_iter()
            .chain(more.into_iter().flatten().map(|&k| &self.placements[k as usize]))
    }

    /// The walk's nodes that overlap the aligned part [ts, te) of the walk and lie on the reference, in walk order.
    fn anchors(&self, walk: &[(usize, bool)], ts: u64, te: u64) -> Result<Vec<Anchor>, Unplaced> {
        if walk.is_empty() || te <= ts {
            return Err(Unplaced::NoAnchor);
        }
        let mut anchors: Vec<Anchor> = Vec::new();
        let mut off = 0u64;
        for &(node, reverse) in walk {
            let len = *self.node_len.get(node).ok_or(Unplaced::MissingNode)?;
            if len == NONE {
                return Err(Unplaced::MissingNode);
            }
            let len = len as u64;
            let next = off + len;
            if len != 0 && next > ts && off < te {
                let local_start = ts.saturating_sub(off);
                let local_end = if te < next { te - off } else { len };
                if local_start < local_end && self.first[node] != NONE {
                    anchors.push(Anchor { node, reverse, local_start, local_end });
                }
            }
            off = next;
            if off >= te {
                break;
            }
        }
        if anchors.is_empty() { Err(Unplaced::NoAnchor) } else { Ok(anchors) }
    }

    /// The reference bases an anchor's aligned part covers at placement p.
    fn interval(&self, a: &Anchor, p: &Placement) -> (u64, u64) {
        let len = self.node_len[a.node] as u64;
        if a.reverse == p.reverse {
            (p.start + a.local_start, p.start + a.local_end)
        } else {
            (p.start + (len - a.local_end), p.start + (len - a.local_start))
        }
    }

    /// pggaf's walk_ref_interval. Starting from each placement of the first
    /// anchor, every later anchor must have exactly one placement on the same
    /// path, in the same relative orientation and strictly further along it;
    /// the interval spans the anchors' aligned bases. Any node missing from the
    /// graph, or two starts giving different intervals, leaves the alignment
    /// unplaced, and why.
    fn project(&self, walk: &[(usize, bool)], ts: u64, te: u64) -> Result<(u32, u64, u64), Unplaced> {
        let anchors = self.anchors(walk, ts, te)?;
        let a0 = &anchors[0];
        let try_from = |p0: &Placement| -> Result<(u32, u64, u64), Unplaced> {
            let forward = a0.reverse == p0.reverse;
            let mut prev = p0.step;
            let (mut lo, mut hi) = self.interval(a0, p0);
            for a in &anchors[1..] {
                let mut found: Option<&Placement> = None;
                for p in self.placements_of(a.node) {
                    if p.path != p0.path || (a.reverse == p.reverse) != forward {
                        continue;
                    }
                    if (forward && p.step <= prev) || (!forward && p.step >= prev) {
                        continue;
                    }
                    if found.is_some() {
                        return Err(Unplaced::Ambiguous);
                    }
                    found = Some(p);
                }
                let p = found.ok_or(Unplaced::Inconsistent)?;
                prev = p.step;
                let (s, e) = self.interval(a, p);
                lo = lo.min(s);
                hi = hi.max(e);
            }
            if lo < hi { Ok((self.path_contig[p0.path as usize], lo, hi)) } else { Err(Unplaced::Inconsistent) }
        };
        let mut resolved: Option<(u32, u64, u64)> = None;
        let mut why = Unplaced::Inconsistent;
        for p0 in self.placements_of(a0.node) {
            match try_from(p0) {
                Ok(c) => match resolved {
                    None => resolved = Some(c),
                    Some(r) if r != c => return Err(Unplaced::Ambiguous),
                    _ => {}
                },
                Err(Unplaced::Ambiguous) => why = Unplaced::Ambiguous,
                Err(_) => {}
            }
        }
        resolved.ok_or(why)
    }
}

/// How far, in graph bases, the search for a reference node next to an alignment without anchors goes from each
/// end of its walk (RefMap::nearest), and how many nodes it visits at most.
const NEAREST_REACH: u64 = 200_000;
const NEAREST_NODES: usize = 200_000;

impl RefMap {
    /// The reference beside a walk with no anchor (Unplaced::NoAnchor), before and after its aligned part [ts, te). A
    /// reference node of the walk itself outside that part is the side, nearest the part first; else the graph is
    /// searched outward from that end of the walk -- forward from its last node, backward from its first -- for the
    /// reference node fewest graph bases away, within NEAREST_REACH and NEAREST_NODES. Only a node the reference visits
    /// once is a side: one it visits more often does not say where the read is, and leaves that side unknown. For each
    /// side, the contig, the reference position where the walk or the search joins the reference, and the bases
    /// between; `next(node, reverse)` is the handles that follow a handle.
    fn nearest(&self, walk: &[(usize, bool)], ts: u64, te: u64, next: &dyn Fn(usize, bool) -> Vec<(usize, bool)>) -> [Option<(u32, u64, u64)>; 2] {
        // The node's place on the reference if the reference visits it once; None if never or more often.
        let once = |node: usize| -> Option<&Placement> {
            let mut all = self.placements_of(node);
            let p = all.next()?;
            all.next().is_none().then_some(p)
        };
        let search = |from: (usize, bool)| -> Option<(u32, u64, u64)> {
            let mut heap = BinaryHeap::new();
            let mut seen: crate::mixhash::MixMap<(usize, bool), u64> = Default::default();
            for h in next(from.0, from.1) {
                heap.push(std::cmp::Reverse((0u64, h)));
            }
            while let Some(std::cmp::Reverse((d, (node, reverse)))) = heap.pop() {
                // Popped nearest first: past the reach or the node budget, nothing nearer is left.
                if d > NEAREST_REACH || seen.len() >= NEAREST_NODES {
                    break;
                }
                if seen.contains_key(&(node, reverse)) {
                    continue;
                }
                seen.insert((node, reverse), d);
                if self.first.get(node).is_some_and(|&f| f != NONE) {
                    let p = once(node)?;
                    let len = self.node_len[node] as u64;
                    let at = if reverse == p.reverse { p.start } else { p.start + len };
                    return Some((self.path_contig[p.path as usize], at, d));
                }
                let len = self.node_len.get(node).copied().filter(|&l| l != NONE)? as u64;
                for h in next(node, reverse) {
                    heap.push(std::cmp::Reverse((d + len, h)));
                }
            }
            None
        };
        // The walk's own reference nodes before and after [ts, te): the last before it, the first after it.
        let mut own: [Option<(u32, u64, u64)>; 2] = [None, None];
        let mut off = 0u64;
        for &(node, reverse) in walk {
            let Some(&len) = self.node_len.get(node) else { break };
            if len == NONE {
                break;
            }
            let end = off + len as u64;
            if let Some(p) = once(node) {
                let contig = self.path_contig[p.path as usize];
                if end <= ts {
                    let at = if reverse == p.reverse { p.start + len as u64 } else { p.start };
                    own[0] = Some((contig, at, ts - end));
                } else if off >= te && own[1].is_none() {
                    let at = if reverse == p.reverse { p.start } else { p.start + len as u64 };
                    own[1] = Some((contig, at, off - te));
                }
            }
            off = end;
        }
        let (Some(&first), Some(&last)) = (walk.first(), walk.last()) else { return [None, None] };
        [own[0].or_else(|| search((first.0, !first.1))), own[1].or_else(|| search(last))]
    }
}

/// A walk node on the reference, and the part of it the alignment covers (in walk orientation).
struct Anchor {
    node: usize,
    reverse: bool,
    local_start: u64,
    local_end: u64,
}

/// Why an alignment has no place.
#[derive(Default, Clone, Copy)]
struct Counts {
    lines: u64,
    placed: u64,
    unmapped: u64,
    unplaced: u64,
    /// Of the unplaced, by reason (Unplaced::ALL's order).
    why: [u64; 4],
    /// Of those without anchors, written to --unplaced with the reference found beside them.
    beside: u64,
    /// Of the placed, those without anchors placed by --rescue.
    rescued: u64,
    /// Of the unplaced without anchors, those --rescue aligned where the reference is at full depth (rescue::Depth).
    held: u64,
    /// And those it aligned where the reads it places would bring the reference to three copies' depth.
    copies: u64,
    malformed: u64,
}

impl Counts {
    fn add(mut self, o: Counts) -> Counts {
        self.lines += o.lines;
        self.placed += o.placed;
        self.unmapped += o.unmapped;
        self.unplaced += o.unplaced;
        self.beside += o.beside;
        self.rescued += o.rescued;
        self.held += o.held;
        self.copies += o.copies;
        for (a, b) in self.why.iter_mut().zip(o.why) {
            *a += b;
        }
        self.malformed += o.malformed;
        self
    }
}

fn parse_u64(b: &[u8]) -> Option<u64> {
    if b.is_empty() || b.len() > 19 {
        return None;
    }
    b.iter().try_fold(0u64, |v, &c| c.is_ascii_digit().then(|| v * 10 + (c - b'0') as u64))
}

fn parse_walk(b: &[u8], out: &mut Vec<(usize, bool)>) -> bool {
    out.clear();
    let mut i = 0;
    while i < b.len() {
        let reverse = match b[i] {
            b'>' => false,
            b'<' => true,
            _ => return false,
        };
        i += 1;
        let s = i;
        while i < b.len() && b[i] != b'>' && b[i] != b'<' {
            i += 1;
        }
        match parse_u64(&b[s..i]) {
            Some(n) => out.push((n as usize, reverse)),
            None => return false,
        }
    }
    true
}


/// A placed alignment in a run: its key and where its line is.
#[derive(Clone, Copy)]
struct Rec {
    rank: u32,
    start: u64,
    end: u64,
    off: u64,
    len: u32,
}

/// One input line: placed (key), or the reason it is not.
enum Line {
    /// Contig rank, reference start and end, and the graph bases aligned.
    Placed(u32, u64, u64, u64),
    Skip,
    Unmapped,
    Unplaced(Unplaced),
    Malformed,
}

/// The aligned part of a GAF line's walk: path start and path end (columns 8 and 9).
fn path_range(line: &[u8]) -> Option<(u64, u64)> {
    let mut f = line.split(|&c| c == b'\t').skip(7);
    Some((parse_u64(f.next()?)?, parse_u64(f.next()?)?))
}

/// The graph bases a GAF line's alignment covers (path end less path start).
fn walk_aligned(line: &[u8]) -> Option<u64> {
    let (ts, te) = path_range(line)?;
    te.checked_sub(ts)
}

fn classify(m: &RefMap, rank: &[u32], line: &[u8], walk: &mut Vec<(usize, bool)>) -> Line {
    if line.is_empty() || line[0] == b'#' || line[0] == b'@' {
        return Line::Skip;
    }
    let mut f: [&[u8]; 12] = [&[]; 12];
    let mut it = line.split(|&c| c == b'\t');
    for x in f.iter_mut() {
        match it.next() {
            Some(v) => *x = v,
            None => return Line::Malformed,
        }
    }
    if f[4] == b"*" {
        return Line::Unmapped;
    }
    if f[4] != b"+" && f[4] != b"-" {
        return Line::Malformed;
    }
    let (Some(ts), Some(te)) = (parse_u64(f[7]), parse_u64(f[8])) else { return Line::Malformed };
    if !parse_walk(f[5], walk) {
        return Line::Malformed;
    }
    match m.project(walk, ts, te) {
        Ok((contig, s, e)) => Line::Placed(rank[contig as usize], s, e, te - ts),
        Err(why) => Line::Unplaced(why),
    }
}

/// An alignment without anchors aligned to the reference by rescue::rescue, in the stretch rescue::stretch gives for
/// the reference beside its walk (`sides`): its contig and the alignment.
fn rescue_line(line: &[u8], walk: &[(usize, bool)], sides: &[Option<(u32, u64, u64)>; 2], letters: &[Vec<u8>], index: &crate::rescue::Index, seq: &dyn Fn(usize) -> Option<Vec<u8>>) -> Option<(u32, crate::rescue::Rescue)> {
    let text = std::str::from_utf8(line).ok()?;
    let pstart: usize = text.split('\t').nth(7)?.parse().ok()?;
    let cs = text.split('\t').find_map(|t| t.strip_prefix("cs:Z:"))?;
    let (contig, lo, hi) = crate::rescue::stretch(sides)?;
    let read = crate::rescue::spell(walk, seq, pstart, cs)?;
    let r = crate::rescue::rescue(&read, index, contig, &letters[contig as usize], lo, hi)?;
    Some((contig, r))
}

/// The pb:Z: and pa:Z: tags of an unplaced alignment without anchors (RefMap::nearest's sides).
fn beside_tags(m: &RefMap, sides: &[Option<(u32, u64, u64)>; 2]) -> String {
    let mut tags = String::new();
    for (tag, x) in ["pb", "pa"].iter().zip(sides) {
        if let Some((contig, at, d)) = x {
            tags.push_str(&format!("\t{tag}:Z:{}:{at}:{d}", m.contigs[*contig as usize]));
        }
    }
    tags
}

/// An alignment without anchors set aside for --rescue: its line, walk and the reference beside it.
type Pending = (Vec<u8>, Vec<(usize, bool)>, [Option<(u32, u64, u64)>; 2]);

/// The input, decompressed if gzip or BGZF.
fn open_input(path: &str) -> Result<Box<dyn Read + Send>, String> {
    let mut f = BufReader::with_capacity(1 << 20, std::fs::File::open(path).map_err(|e| format!("{path}: {e}"))?);
    let head = f.fill_buf().map_err(|e| format!("{path}: {e}"))?;
    let gz = head.starts_with(&[0x1f, 0x8b]);
    // BGZF: gzip members with the extra field BC (the SAM spec's), each block a deflate stream of its own.
    let bgzf = gz && head.len() >= 16 && head[3] & 4 != 0 && head[12] == b'B' && head[13] == b'C';
    Ok(if bgzf {
        Box::new(BgzfReader { f, buf: Vec::new(), pos: 0 })
    } else if gz {
        Box::new(flate2::bufread::MultiGzDecoder::new(f))
    } else {
        Box::new(f)
    })
}

/// Blocks of a BGZF file decompressed at a time, on the thread pool.
const BGZF_BATCH: usize = 512;

/// A BGZF file as its letters, BGZF_BATCH blocks decompressed in parallel at a time, each block's CRC checked.
struct BgzfReader {
    f: BufReader<std::fs::File>,
    buf: Vec<u8>,
    pos: usize,
}

impl BgzfReader {
    /// The next block, whole (header to ISIZE), or None at the end of the file.
    fn block(&mut self) -> std::io::Result<Option<Vec<u8>>> {
        let mut head = [0u8; 12];
        match self.f.read_exact(&mut head) {
            Ok(()) => {}
            Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => return Ok(None),
            Err(e) => return Err(e),
        }
        let bad = |what: &str| std::io::Error::new(std::io::ErrorKind::InvalidData, format!("BGZF: {what}"));
        if head[0] != 0x1f || head[1] != 0x8b || head[3] & 4 == 0 {
            return Err(bad("not a BGZF block"));
        }
        let xlen = u16::from_le_bytes([head[10], head[11]]) as usize;
        let mut extra = vec![0u8; xlen];
        self.f.read_exact(&mut extra)?;
        let mut bsize = None;
        let mut i = 0;
        while i + 4 <= xlen {
            let slen = u16::from_le_bytes([extra[i + 2], extra[i + 3]]) as usize;
            if extra[i] == b'B' && extra[i + 1] == b'C' && slen == 2 && i + 6 <= xlen {
                bsize = Some(u16::from_le_bytes([extra[i + 4], extra[i + 5]]) as usize + 1);
            }
            i += 4 + slen;
        }
        let total = bsize.ok_or_else(|| bad("no BC field"))?;
        let rest = total.checked_sub(12 + xlen).filter(|&r| r >= 8).ok_or_else(|| bad("block size"))?;
        let mut block = Vec::with_capacity(total);
        block.extend_from_slice(&head);
        block.extend_from_slice(&extra);
        block.resize(total, 0);
        self.f.read_exact(&mut block[12 + xlen..])?;
        debug_assert_eq!(block.len(), 12 + xlen + rest);
        Ok(Some(block))
    }

    /// Decompresses the next batch of blocks into `buf`; false at the end of the file.
    fn fill(&mut self) -> std::io::Result<bool> {
        let mut blocks = Vec::with_capacity(BGZF_BATCH);
        while blocks.len() < BGZF_BATCH {
            match self.block()? {
                Some(b) => blocks.push(b),
                None => break,
            }
        }
        if blocks.is_empty() {
            return Ok(false);
        }
        let letters: Vec<std::io::Result<Vec<u8>>> = blocks.par_iter().map(|b| inflate_block(b)).collect();
        self.buf.clear();
        self.pos = 0;
        for l in letters {
            self.buf.extend_from_slice(&l?);
        }
        Ok(true)
    }
}

/// One BGZF block's letters: its deflate stream inflated and checked against its CRC32 and ISIZE.
fn inflate_block(block: &[u8]) -> std::io::Result<Vec<u8>> {
    let bad = |what: &str| std::io::Error::new(std::io::ErrorKind::InvalidData, format!("BGZF: {what}"));
    let xlen = u16::from_le_bytes([block[10], block[11]]) as usize;
    let n = block.len();
    let size = u32::from_le_bytes(block[n - 4..].try_into().unwrap()) as usize;
    let crc = u32::from_le_bytes(block[n - 8..n - 4].try_into().unwrap());
    let mut out = Vec::with_capacity(size);
    let mut d = flate2::Decompress::new(false);
    d.decompress_vec(&block[12 + xlen..n - 8], &mut out, flate2::FlushDecompress::Finish).map_err(|e| bad(&e.to_string()))?;
    if out.len() != size {
        return Err(bad("block length"));
    }
    let mut c = flate2::Crc::new();
    c.update(&out);
    if c.sum() != crc {
        return Err(bad("CRC"));
    }
    Ok(out)
}

impl Read for BgzfReader {
    fn read(&mut self, out: &mut [u8]) -> std::io::Result<usize> {
        while self.pos == self.buf.len() {
            if !self.fill()? {
                return Ok(0);
            }
        }
        let k = out.len().min(self.buf.len() - self.pos);
        out[..k].copy_from_slice(&self.buf[self.pos..self.pos + k]);
        self.pos += k;
        Ok(k)
    }
}

/// Reads the input in chunks of whole lines on a thread of its own, so
/// decompression overlaps the projection. Every chunk ends with a newline.
fn read_chunks(path: &str) -> Result<std::sync::mpsc::Receiver<std::io::Result<Vec<u8>>>, String> {
    const CHUNK: usize = 32 << 20;
    let mut input = open_input(path)?;
    let (tx, rx) = std::sync::mpsc::sync_channel(3);
    std::thread::spawn(move || {
        let mut buf: Vec<u8> = Vec::new();
        loop {
            // Fill up to CHUNK more bytes after what is carried over.
            let have = buf.len();
            buf.resize(have + CHUNK, 0);
            let mut n = have;
            let mut eof = false;
            while n < buf.len() {
                match input.read(&mut buf[n..]) {
                    Ok(0) => {
                        eof = true;
                        break;
                    }
                    Ok(k) => n += k,
                    Err(e) if e.kind() == std::io::ErrorKind::Interrupted => {}
                    Err(e) => {
                        let _ = tx.send(Err(e));
                        return;
                    }
                }
            }
            buf.truncate(n);
            if eof {
                if !buf.is_empty() {
                    if buf.last() != Some(&b'\n') {
                        buf.push(b'\n');
                    }
                    let _ = tx.send(Ok(buf));
                }
                return;
            }
            // Send the whole lines; carry the partial last one (a line longer
            // than a chunk is carried until it is whole).
            let Some(nl) = buf.iter().rposition(|&c| c == b'\n') else { continue };
            let carry = buf.split_off(nl + 1);
            if tx.send(Ok(std::mem::replace(&mut buf, carry))).is_err() {
                return;
            }
        }
    });
    Ok(rx)
}

/// Reference bases per bucket of a spilled run: each bucket is a zstd frame of its own, and the merge sorts a bucket
/// of every run at a time, on a thread of its own.
const BUCKET: u64 = 1 << 20;

/// Sorted runs: kept in memory, or spilled as zstd files.
struct Runs {
    dir: PathBuf,
    spilled: Vec<Spilled>,
}

/// A spilled run: its file, and where each bucket's frame lies in it, by (contig rank, start / BUCKET), in order.
struct Spilled {
    path: PathBuf,
    frames: Vec<((u32, u64), u64, u64)>,
}

impl Runs {
    /// Writes `recs` (sorted, not empty) as the next run.
    fn spill(&mut self, data: &[u8], recs: &[Rec]) -> std::io::Result<()> {
        let path = self.dir.join(format!("run{}.zst", self.spilled.len() + 1));
        let run = spill_run(&self.dir, path, data, recs)?;
        self.spilled.push(run);
        Ok(())
    }
}

/// Writes `recs` (sorted, not empty) as a run at `path` in `dir`: per record its rank, start, end and length, then its
/// line, one zstd frame per bucket.
fn spill_run(dir: &Path, path: PathBuf, data: &[u8], recs: &[Rec]) -> std::io::Result<Spilled> {
    std::fs::create_dir_all(dir)?;
    let bucket = |r: &Rec| (r.rank, r.start / BUCKET);
    let mut cuts = vec![0];
    cuts.extend((1..recs.len()).filter(|&i| bucket(&recs[i]) != bucket(&recs[i - 1])));
    cuts.push(recs.len());
    let frames: Vec<std::io::Result<Vec<u8>>> = cuts
        .par_windows(2)
        .map(|w| {
            let mut raw = Vec::new();
            for r in &recs[w[0]..w[1]] {
                raw.extend_from_slice(&r.rank.to_le_bytes());
                raw.extend_from_slice(&r.start.to_le_bytes());
                raw.extend_from_slice(&r.end.to_le_bytes());
                raw.extend_from_slice(&r.len.to_le_bytes());
                raw.extend_from_slice(&data[r.off as usize..r.off as usize + r.len as usize]);
            }
            zstd::bulk::compress(&raw, 1)
        })
        .collect();
    let mut out = BufWriter::with_capacity(1 << 22, std::fs::File::create(&path)?);
    let mut run = Spilled { path, frames: Vec::with_capacity(frames.len()) };
    let mut at = 0u64;
    for (w, frame) in cuts.windows(2).zip(frames) {
        let frame = frame?;
        out.write_all(&frame)?;
        run.frames.push((bucket(&recs[w[0]]), at, frame.len() as u64));
        at += frame.len() as u64;
    }
    out.flush()?;
    Ok(run)
}

impl Drop for Runs {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.dir);
    }
}

fn sort_run(data: &[u8], recs: &mut [Rec]) {
    let line = |r: &Rec| &data[r.off as usize..r.off as usize + r.len as usize];
    recs.par_sort_unstable_by(|a, b| {
        (a.rank, a.start, a.end).cmp(&(b.rank, b.start, b.end)).then_with(|| line(a).cmp(line(b)))
    });
}

/// One bucket of the merge: its lines in order, laid out as `emit` writes them, and each line's contig, start, end
/// and length.
struct Bucket {
    text: Vec<u8>,
    recs: Vec<(u32, u64, u64, usize)>,
}

/// Bucket `key` of every spilled run (`files` open on them), its records sorted as `sort_run` sorts them.
fn merge_bucket(runs: &[Spilled], files: &[std::fs::File], key: (u32, u64), names: &[String], by_rank: &[u32]) -> Result<Bucket, String> {
    use std::os::unix::fs::FileExt;
    let mut raw: Vec<Vec<u8>> = Vec::new();
    for (run, f) in runs.iter().zip(files) {
        let Ok(i) = run.frames.binary_search_by(|x| x.0.cmp(&key)) else { continue };
        let (_, off, len) = run.frames[i];
        let mut frame = vec![0u8; len as usize];
        f.read_exact_at(&mut frame, off).map_err(|e| format!("{}: {e}", run.path.display()))?;
        raw.push(zstd::stream::decode_all(&frame[..]).map_err(|e| format!("{}: {e}", run.path.display()))?);
    }
    let mut recs: Vec<((u32, u64, u64), &[u8])> = Vec::new();
    for r in &raw {
        let mut p = 0;
        while p < r.len() {
            let rank = u32::from_le_bytes(r[p..p + 4].try_into().unwrap());
            let start = u64::from_le_bytes(r[p + 4..p + 12].try_into().unwrap());
            let end = u64::from_le_bytes(r[p + 12..p + 20].try_into().unwrap());
            let len = u32::from_le_bytes(r[p + 20..p + 24].try_into().unwrap()) as usize;
            recs.push(((rank, start, end), &r[p + 24..p + 24 + len]));
            p += 24 + len;
        }
    }
    recs.sort_unstable();
    let mut b = Bucket { text: Vec::new(), recs: Vec::with_capacity(recs.len()) };
    for (key, line) in recs {
        let at = b.text.len();
        let contig = lay_out(&mut b.text, names, by_rank, key, line);
        b.recs.push((contig, key.1, key.2, b.text.len() - at));
    }
    Ok(b)
}

/// Appends one placed alignment's output line to `buf`: `contig  start+1  end  line`. Returns its contig.
fn lay_out(buf: &mut Vec<u8>, names: &[String], by_rank: &[u32], key: (u32, u64, u64), line: &[u8]) -> u32 {
    let contig = by_rank[key.0 as usize];
    buf.extend_from_slice(names[contig as usize].as_bytes());
    buf.extend_from_slice(format!("\t{}\t{}\t", key.1 + 1, key.2).as_bytes());
    buf.extend_from_slice(line);
    buf.push(b'\n');
    contig
}

/// Writes one placed alignment (`lay_out`).
fn emit(w: &mut IndexedWriter, names: &[String], by_rank: &[u32], key: (u32, u64, u64), line: &[u8], buf: &mut Vec<u8>) -> std::io::Result<()> {
    buf.clear();
    let contig = lay_out(buf, names, by_rank, key, line);
    w.add(buf, contig, key.1, key.2)
}

pub fn run(o: &IndexGafOptions, sample: &str) -> Result<(), String> {
    let t0 = Instant::now();
    let log = |msg: &str| eprintln!("[pgcat {:7.1}s] {msg}", t0.elapsed().as_secs_f64());
    // The reference map from the graph, or from a saved ref map without the graph (which --unplaced and --rescue
    // search, so they take the graph).
    let (g, m, sample) = match &o.ref_map {
        Some(path) => {
            if o.unplaced.is_some() || o.rescue {
                return Err("--ref-map places alignments without the graph; --unplaced and --rescue search it (give --graph)".into());
            }
            let (m, made_for) = RefMap::load(path)?;
            if !sample.is_empty() && sample != made_for {
                return Err(format!("{path}: the ref map of {made_for}, not of {sample}"));
            }
            (None, m, made_for)
        }
        None => {
            let g: GBZ = simple_sds::serialize::load_from(&o.graph).map_err(|e| format!("{}: {e}", o.graph))?;
            let m = RefMap::build(&g, sample)?;
            (Some(g), m, sample.to_string())
        }
    };
    let sample = sample.as_str();
    // The graph stays loaded only to find the reference beside the alignments without anchors, and the reference's
    // letters are read only to align them there.
    let letters = if o.rescue { m.letters(g.as_ref().unwrap()) } else { Vec::new() };
    let g = if o.unplaced.is_some() || o.rescue { g } else { None };
    let seq = |node: usize| g.as_ref().and_then(|g| g.sequence(node)).map(|s| s.to_vec());
    let next = |node: usize, reverse: bool| -> Vec<(usize, bool)> {
        let o = if reverse { Orientation::Reverse } else { Orientation::Forward };
        let edges = g.as_ref().and_then(|g| g.successors(node, o));
        edges.map(|it| it.map(|(n, o)| (n, o == Orientation::Reverse)).collect()).unwrap_or_default()
    };
    log(&format!("{} loaded: {} reference paths of {sample} on {} contigs", if o.ref_map.is_some() { "reference map" } else { "graph" }, m.path_contig.len(), m.contigs.len()));
    // Output order is by contig name, then position, then the line itself.
    let mut by_rank: Vec<u32> = (0..m.contigs.len() as u32).collect();
    by_rank.sort_by(|&a, &b| m.contigs[a as usize].as_bytes().cmp(m.contigs[b as usize].as_bytes()));
    let mut rank = vec![0u32; m.contigs.len()];
    for (r, &c) in by_rank.iter().enumerate() {
        rank[c as usize] = r as u32;
    }

    let tmp_parent = o.tmp_dir.clone().unwrap_or_else(|| {
        Path::new(&o.out).parent().filter(|p| !p.as_os_str().is_empty()).map(|p| p.to_string_lossy().into_owned()).unwrap_or_else(|| ".".into())
    });
    crate::build::remove_abandoned_work(&tmp_parent, ".pgcat-sort-");
    let mut runs = Runs { dir: PathBuf::from(&tmp_parent).join(format!(".pgcat-sort-{}", std::process::id())), spilled: Vec::new() };
    let mut data: Vec<u8> = Vec::new();
    let mut recs: Vec<Rec> = Vec::new();
    let mut counts = Counts::default();
    // The most graph bases any alignment on each contig (by rank) spans.
    let mut max_span = vec![0u64; m.contigs.len()];
    let mut unplaced_out = match &o.unplaced {
        Some(path) => Some(flate2::write::GzEncoder::new(
            BufWriter::with_capacity(1 << 20, std::fs::File::create(path).map_err(|e| format!("{path}: {e}"))?),
            flate2::Compression::new(o.level),
        )),
        None => None,
    };
    let mut pending: Vec<Pending> = Vec::new();
    // The sorted run being spilled on its own thread, if any.
    let mut spilling: Option<std::thread::JoinHandle<Result<Spilled, String>>> = None;
    let mut depth = crate::rescue::Depth::default();
    let threads = rayon::current_num_threads();
    for chunk in read_chunks(&o.gaf)? {
        let chunk = chunk.map_err(|e| format!("{}: {e}", o.gaf))?;
        // Split at line ends into pieces for the threads.
        let mut pieces: Vec<&[u8]> = Vec::new();
        let target = chunk.len() / (threads * 4) + 1;
        let mut s = 0;
        while s < chunk.len() {
            let mut e = (s + target).min(chunk.len());
            while e < chunk.len() && chunk[e - 1] != b'\n' {
                e += 1;
            }
            pieces.push(&chunk[s..e]);
            s = e;
        }
        let keep_unplaced = unplaced_out.is_some();
        let results: Vec<(Vec<u8>, Vec<Rec>, Counts, Vec<u64>, Vec<(usize, usize, Unplaced, String)>, Vec<Pending>)> = pieces
            .par_iter()
            .map(|piece| {
                let base = piece.as_ptr() as usize - chunk.as_ptr() as usize;
                let mut walk = Vec::new();
                // The piece's placed lines as the sorted runs hold them, offsets from the piece's own bytes.
                let mut bytes: Vec<u8> = Vec::new();
                let mut recs_here: Vec<Rec> = Vec::new();
                let mut dropped = Vec::new();
                let mut aside: Vec<Pending> = Vec::new();
                let mut c = Counts::default();
                let mut span = vec![0u64; rank.len()];
                let mut p = 0;
                for line in piece.split_inclusive(|&b| b == b'\n') {
                    let body = &line[..line.len() - 1];
                    let body = body.strip_suffix(b"\r").unwrap_or(body);
                    match classify(&m, &rank, body, &mut walk) {
                        Line::Placed(r, s, e, aligned) => {
                            c.lines += 1;
                            c.placed += 1;
                            span[r as usize] = span[r as usize].max(aligned);
                            let off = bytes.len();
                            append_line(&mut bytes, body, o.slim);
                            recs_here.push(Rec { rank: r, start: s, end: e, off: off as u64, len: (bytes.len() - off) as u32 });
                        }
                        Line::Skip => {}
                        Line::Unmapped => {
                            c.lines += 1;
                            c.unmapped += 1;
                        }
                        Line::Unplaced(why) => {
                            c.lines += 1;
                            let sides = if why == Unplaced::NoAnchor && (keep_unplaced || o.rescue) {
                                let (ts, te) = path_range(body).unwrap_or((0, 0));
                                m.nearest(&walk, ts, te, &next)
                            } else {
                                [None, None]
                            };
                            if o.rescue && why == Unplaced::NoAnchor {
                                // Rescued after the pass, from one index over every such alignment's stretch.
                                aside.push((body.to_vec(), walk.clone(), sides));
                                p += line.len();
                                continue;
                            }
                            c.unplaced += 1;
                            c.why[Unplaced::ALL.iter().position(|&w| w == why).unwrap()] += 1;
                            if keep_unplaced {
                                let mut tags = format!("\tpu:Z:{}", why.name());
                                if why == Unplaced::NoAnchor {
                                    c.beside += sides.iter().any(|x| x.is_some()) as u64;
                                    tags.push_str(&beside_tags(&m, &sides));
                                }
                                dropped.push((base + p, body.len(), why, tags));
                            }
                        }
                        Line::Malformed => {
                            c.lines += 1;
                            c.malformed += 1;
                        }
                    }
                    p += line.len();
                }
                (bytes, recs_here, c, span, dropped, aside)
            })
            .collect();
        for (bytes, recs_here, c, span, dropped, aside) in results {
            pending.extend(aside);
            if let Some(w) = unplaced_out.as_mut() {
                for (at, len, _, tags) in dropped {
                    w.write_all(&chunk[at..at + len]).map_err(|e| format!("unplaced: {e}"))?;
                    w.write_all(tags.as_bytes()).and_then(|_| w.write_all(b"\n")).map_err(|e| format!("unplaced: {e}"))?;
                }
            }
            counts = counts.add(c);
            for (m, s) in max_span.iter_mut().zip(span) {
                *m = (*m).max(s);
            }
            let base = data.len() as u64;
            data.extend_from_slice(&bytes);
            for r in recs_here {
                if o.rescue {
                    depth.add(r.rank as usize, r.start, r.end);
                }
                recs.push(Rec { off: r.off + base, ..r });
            }
        }
        // A full buffer is sorted and spilled on a thread of its own while reading goes on into the next: each holds
        // half the memory budget.
        if data.len() + recs.len() * std::mem::size_of::<Rec>() >= o.memory / 2 {
            if let Some(h) = spilling.take() {
                let run = h.join().map_err(|_| "a sorted run's thread failed".to_string())??;
                runs.spilled.push(run);
            }
            let (mut d, mut r) = (std::mem::take(&mut data), std::mem::take(&mut recs));
            let path = runs.dir.join(format!("run{}.zst", runs.spilled.len() + 1));
            let dir = runs.dir.clone();
            spilling = Some(std::thread::spawn(move || -> Result<Spilled, String> {
                sort_run(&d, &mut r);
                let run = spill_run(&dir, path, &d, &r).map_err(|e| format!("{}: {e}", dir.display()));
                d.clear();
                r.clear();
                run
            }));
            log(&format!("{} alignments read; sorted run {} spilling", counts.lines, runs.spilled.len() + 1));
        }
    }
    if let Some(h) = spilling.take() {
        let run = h.join().map_err(|_| "a sorted run's thread failed".to_string())??;
        runs.spilled.push(run);
    }
    if !pending.is_empty() {
        depth.finish();
        let stretches: Vec<(u32, u64, u64)> = pending.iter().filter_map(|p| crate::rescue::stretch(&p.2)).collect();
        let index = crate::rescue::Index::build(&letters, &stretches);
        log(&format!("rescue: {} alignments without anchors, an index of {} reference minimizers", pending.len(), index.len()));
        let done: Vec<Option<(u32, crate::rescue::Rescue)>> =
            pending.par_iter().map(|(line, walk, sides)| rescue_line(line, walk, sides, &letters, &index, &seq)).collect();
        let at: Vec<Option<(usize, u64, u64)>> =
            done.iter().map(|r| r.as_ref().map(|(contig, r)| (rank[*contig as usize] as usize, r.pos, r.pos + r.span()))).collect();
        // Held where the reference is at full depth already; then the depth with the rest placed, against the median
        // before them.
        let full: Vec<bool> = at.iter().map(|a| a.is_some_and(|(c, s, e)| depth.full(c, s, e))).collect();
        let mut after = depth.clone();
        for (a, &f) in at.iter().zip(&full) {
            if let (Some((c, s, e)), false) = (*a, f) {
                after.add(c, s, e);
            }
        }
        for ((((line, _, sides), r), a), f) in pending.into_iter().zip(done).zip(at).zip(full) {
            let r = if f {
                counts.held += 1;
                None
            } else if a.is_some_and(|(c, s, e)| after.share(c, s, e) >= crate::rescue::MAX_RESCUED_DEPTH) {
                counts.copies += 1;
                None
            } else {
                r
            };
            match r {
                Some((contig, r)) => {
                    let rk = rank[contig as usize];
                    counts.placed += 1;
                    counts.rescued += 1;
                    max_span[rk as usize] = max_span[rk as usize].max(walk_aligned(&line).unwrap_or(0));
                    let off = data.len();
                    append_line(&mut data, &line, o.slim);
                    data.extend_from_slice(format!("\tpu:Z:no_anchor\trs:Z:{}\trq:i:{}", r.tag(&m.contigs[contig as usize]), r.mapq).as_bytes());
                    recs.push(Rec { rank: rk, start: r.pos, end: r.pos + r.span(), off: off as u64, len: (data.len() - off) as u32 });
                }
                None => {
                    counts.unplaced += 1;
                    counts.why[0] += 1;
                    if let Some(w) = unplaced_out.as_mut() {
                        counts.beside += sides.iter().any(|x| x.is_some()) as u64;
                        let tags = format!("\tpu:Z:no_anchor{}\n", beside_tags(&m, &sides));
                        w.write_all(&line).and_then(|_| w.write_all(tags.as_bytes())).map_err(|e| format!("unplaced: {e}"))?;
                    }
                }
            }
        }
    }
    sort_run(&data, &mut recs);
    if let Some(w) = unplaced_out.take() {
        w.finish().and_then(|mut b| b.flush()).map_err(|e| format!("unplaced: {e}"))?;
        log(&format!("{} of {} alignments without anchors have the reference beside them (pb/pa tags)", counts.beside, counts.why[0]));
    }
    let why: Vec<String> = Unplaced::ALL.iter().zip(counts.why).map(|(w, n)| format!("{} {n}", w.name())).collect();
    drop(g);
    drop(letters);
    log(&format!(
        "{} alignments: {} placed ({} without anchors, by --rescue; {} more aligned where the reference is at full depth, {} where the reads placed would bring it to three copies' depth), {} unmapped, {} without a unique reference interval ({}), {} malformed",
        counts.lines, counts.placed, counts.rescued, counts.held, counts.copies, counts.unmapped, counts.unplaced, why.join(", "), counts.malformed
    ));

    let tmp_out = format!("{}.tmp", o.out);
    let mut w = IndexedWriter::create(&tmp_out, m.contigs.clone(), o.level).map_err(|e| format!("{tmp_out}: {e}"))?;
    let mut buf = Vec::new();
    let io = |e: std::io::Error| format!("{tmp_out}: {e}");
    // The longest alignment on each contig, for pgphase: a site whose every
    // allele spans more graph bases cannot be traversed whole by any read.
    let mut header = String::new();
    // Every contig is listed, 0 for one without alignments.
    for (r, &c) in by_rank.iter().enumerate() {
        header.push_str(&format!("#max_span\t{}\t{}\n", m.contigs[c as usize], max_span[r]));
    }
    w.header(header.as_bytes()).map_err(io)?;
    if runs.spilled.is_empty() {
        for r in &recs {
            let line = &data[r.off as usize..r.off as usize + r.len as usize];
            emit(&mut w, &m.contigs, &by_rank, (r.rank, r.start, r.end), line, &mut buf).map_err(io)?;
        }
    } else {
        if !recs.is_empty() {
            runs.spill(&data, &recs).map_err(|e| format!("{}: {e}", runs.dir.display()))?;
        }
        drop(std::mem::take(&mut data));
        drop(std::mem::take(&mut recs));
        log(&format!("merging {} sorted runs", runs.spilled.len()));
        let files = runs
            .spilled
            .iter()
            .map(|r| std::fs::File::open(&r.path).map_err(|e| format!("{}: {e}", r.path.display())))
            .collect::<Result<Vec<_>, String>>()?;
        let mut keys: Vec<(u32, u64)> = runs.spilled.iter().flat_map(|r| r.frames.iter().map(|f| f.0)).collect();
        keys.sort_unstable();
        keys.dedup();
        // A window of buckets at a time: sorted in parallel, then written in order (the writer compresses in parallel).
        for window in keys.chunks(rayon::current_num_threads() * 2) {
            let buckets: Vec<Result<Bucket, String>> =
                window.par_iter().map(|&k| merge_bucket(&runs.spilled, &files, k, &m.contigs, &by_rank)).collect();
            for b in buckets {
                let b = b?;
                let mut at = 0;
                for &(contig, beg, end, len) in &b.recs {
                    w.add(&b.text[at..at + len], contig, beg, end).map_err(io)?;
                    at += len;
                }
            }
        }
    }
    let tbi = format!("{}.tbi", o.out);
    w.finish(&format!("{tbi}.part")).map_err(io)?;
    std::fs::rename(&tmp_out, &o.out).map_err(io)?;
    std::fs::rename(format!("{tbi}.part"), &tbi).map_err(|e| format!("{tbi}: {e}"))?;
    log(&format!("{} and {tbi}: {} alignments", o.out, counts.placed));
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bgzf_input_read_in_batches() {
        // More blocks than a batch, written as pgcat's own BGZF writer writes them, read back whole.
        let text: Vec<u8> = (0..3 * BGZF_BATCH * 700).map(|i| b"ACGT\t\n"[i % 6]).collect();
        let path = std::env::temp_dir().join(format!("pgcat-bgzf-test-{}.gz", std::process::id()));
        let mut file: Vec<u8> = Vec::new();
        for chunk in text.chunks(900) {
            file.extend_from_slice(&crate::bgzf::compress_block(chunk, flate2::Compression::new(6)));
        }
        file.extend_from_slice(&crate::bgzf::compress_block(&[], flate2::Compression::new(6)));
        std::fs::write(&path, &file).unwrap();
        let mut back = Vec::new();
        open_input(path.to_str().unwrap()).unwrap().read_to_end(&mut back).unwrap();
        std::fs::remove_file(&path).unwrap();
        assert_eq!(back, text);
        // A corrupted letter fails its block's CRC.
        let mut bad = crate::bgzf::compress_block(b"ACGT", flate2::Compression::new(0));
        let n = bad.len();
        bad[n - 9] ^= 1;
        assert!(inflate_block(&bad).is_err());
    }

    #[test]
    fn slim_lines_keep_what_pgphase_reads() {
        let line = b"r1\t100\t0\t100\t+\t>1>2\t200\t5\t105\t99\t100\t60\tcs:Z::100\tbq:Z:IIII\trq:i:7\tpu:Z:no_anchor";
        let mut out = Vec::new();
        append_line(&mut out, line, true);
        assert_eq!(out, b"r1\t100\t0\t100\t+\t>1>2\t200\t5\t105\t99\t100\t60\trq:i:7\tpu:Z:no_anchor".to_vec());
        out.clear();
        append_line(&mut out, line, false);
        assert_eq!(out, line.to_vec());
    }

    #[test]
    fn bucketed_merge_is_one_sort() {
        // Three runs of records on two contigs over five buckets, merged bucket by bucket, against all of them sorted.
        let mut runs = Runs { dir: std::env::temp_dir().join(format!("pgcat-merge-test-{}", std::process::id())), spilled: Vec::new() };
        let (names, by_rank) = (vec!["chrA".to_string(), "chrB".to_string()], vec![1u32, 0]);
        let mut all: Vec<((u32, u64, u64), Vec<u8>)> = Vec::new();
        for run in 0..3u64 {
            let (mut data, mut recs) = (Vec::new(), Vec::new());
            for i in 0..300u64 {
                let (rank, start) = (((i + run) % 2) as u32, (i * 17_389 + run * 104_729) % (5 * BUCKET));
                // Some lines repeat a key, with other letters.
                let line = format!("read{}_{run}", i % 7).into_bytes();
                recs.push(Rec { rank, start, end: start + 15_000, off: data.len() as u64, len: line.len() as u32 });
                data.extend_from_slice(&line);
                all.push(((rank, start, start + 15_000), line));
            }
            sort_run(&data, &mut recs);
            runs.spill(&data, &recs).unwrap();
        }
        all.sort();
        let files: Vec<std::fs::File> = runs.spilled.iter().map(|r| std::fs::File::open(&r.path).unwrap()).collect();
        let mut keys: Vec<(u32, u64)> = runs.spilled.iter().flat_map(|r| r.frames.iter().map(|f| f.0)).collect();
        keys.sort_unstable();
        keys.dedup();
        assert!(keys.len() >= 8);
        let mut got = Vec::new();
        for k in keys {
            got.extend_from_slice(&merge_bucket(&runs.spilled, &files, k, &names, &by_rank).unwrap().text);
        }
        let mut want = Vec::new();
        for (key, line) in &all {
            lay_out(&mut want, &names, &by_rank, *key, line);
        }
        assert_eq!(got, want);
    }

    fn map() -> RefMap {
        // Reference path: >1 >2 >3 (lengths 10, 5, 10) at 100; node 4 (5 bp)
        // is off the reference; node 2 is also visited by a second path.
        let mut m = RefMap {
            node_len: vec![NONE, 10, 5, 10, 5],
            first: vec![NONE; 5],
            placements: Vec::new(),
            more: HashMap::new(),
            path_contig: vec![0, 1],
            contigs: vec!["chrA".into(), "chrB".into()],
        };
        let add = |node: usize, p: Placement, m: &mut RefMap| {
            let k = m.placements.len() as u32;
            m.placements.push(p);
            if m.first[node] == NONE { m.first[node] = k } else { m.more.entry(node).or_default().push(k) }
        };
        add(1, Placement { path: 0, step: 0, start: 100, reverse: false }, &mut m);
        add(2, Placement { path: 0, step: 1, start: 110, reverse: false }, &mut m);
        add(3, Placement { path: 0, step: 2, start: 115, reverse: false }, &mut m);
        m
    }

    #[test]
    fn projection() {
        let m = map();
        // Forward, all three nodes, aligned from base 3 to 22 of the walk.
        assert_eq!(m.project(&[(1, false), (2, false), (3, false)], 3, 22), Ok((0, 103, 122)));
        // Reverse walk <3 <2 <1: the same bases, read backwards.
        assert_eq!(m.project(&[(3, true), (2, true), (1, true)], 3, 22), Ok((0, 103, 122)));
        // An off-reference node between anchors is skipped.
        assert_eq!(m.project(&[(1, false), (4, false), (3, false)], 0, 25), Ok((0, 100, 125)));
        // Steps going backwards along a forward walk: no place.
        assert_eq!(m.project(&[(3, false), (1, false)], 0, 20), Err(Unplaced::Inconsistent));
        // The orientation flipping part way (an inversion): no place.
        assert_eq!(m.project(&[(1, false), (3, true)], 0, 20), Err(Unplaced::Inconsistent));
        // A node missing from the graph: no place.
        assert_eq!(m.project(&[(1, false), (9, false)], 0, 12), Err(Unplaced::MissingNode));
        // Only off-reference nodes: no place.
        assert_eq!(m.project(&[(4, false)], 0, 5), Err(Unplaced::NoAnchor));
    }

    #[test]
    fn reference_beside_an_unanchored_walk() {
        let mut m = map();
        // Node 7 (10 bp) off the reference, on a path 1 > 4 > 7 > 3 that leaves out node 2.
        m.node_len.extend([NONE, NONE, 10]);
        m.first.extend([NONE, NONE, NONE]);
        let fwd = [(1, 4), (4, 7), (7, 3)];
        let next = |n: usize, r: bool| -> Vec<(usize, bool)> {
            fwd.iter()
                .filter_map(|&(a, b)| if !r && a == n { Some((b, false)) } else if r && b == n { Some((a, true)) } else { None })
                .collect()
        };
        // Before the walk, node 4's 5 bp back to node 1's end; after it, straight onto node 3.
        assert_eq!(m.nearest(&[(7, false)], 0, 10, &next), [Some((0, 110, 5)), Some((0, 115, 0))]);
        // The same walk reversed: the sides swap.
        assert_eq!(m.nearest(&[(7, true)], 0, 10, &next), [Some((0, 115, 0)), Some((0, 110, 5))]);
        // No edges: nothing beside it.
        assert_eq!(m.nearest(&[(7, false)], 0, 10, &|_, _| Vec::new()), [None, None]);
        // The walk's own reference nodes outside its aligned part are the sides, with no search: walk >1 >4 >7 >3
        // aligned on node 4's last base and node 7 (bases 14..25), so node 1 ends 4 bases before and node 3 starts
        // right after.
        assert_eq!(m.nearest(&[(1, false), (4, false), (7, false), (3, false)], 14, 25, &|_, _| Vec::new()),
                   [Some((0, 110, 4)), Some((0, 115, 0))]);
        // Reversed (<3 <7 <4 <1, the same bases): the sides swap.
        assert_eq!(m.nearest(&[(3, true), (7, true), (4, true), (1, true)], 10, 21, &|_, _| Vec::new()),
                   [Some((0, 115, 0)), Some((0, 110, 4))]);
        // A side the reference visits twice says nothing: node 1 on a second path too, the left side is unknown.
        let mut twice = map();
        twice.node_len.extend([NONE, NONE, 10]);
        twice.first.extend([NONE, NONE, NONE]);
        twice.placements.push(Placement { path: 1, step: 0, start: 500, reverse: false });
        twice.more.entry(1).or_default().push(twice.placements.len() as u32 - 1);
        assert_eq!(twice.nearest(&[(7, false)], 0, 10, &next), [None, Some((0, 115, 0))]);
    }

    #[test]
    fn ambiguous_placements() {
        let mut m = map();
        // Node 2 also on path 1 (chrB): a walk of node 2 alone has two places.
        let k = m.placements.len() as u32;
        m.placements.push(Placement { path: 1, step: 0, start: 500, reverse: false });
        m.more.entry(2).or_default().push(k);
        assert_eq!(m.project(&[(2, false)], 0, 5), Err(Unplaced::Ambiguous));
        // With node 1 before it, only path 0 fits.
        assert_eq!(m.project(&[(1, false), (2, false)], 0, 15), Ok((0, 100, 115)));
    }

    #[test]
    fn lines() {
        let m = map();
        let rank = vec![0, 1];
        let mut w = Vec::new();
        let ok = b"r1\t20\t0\t20\t+\t>1>2>3\t25\t2\t22\t20\t20\t60";
        assert!(matches!(classify(&m, &rank, ok, &mut w), Line::Placed(0, 102, 122, 20)));
        assert!(matches!(classify(&m, &rank, b"r2\t20\t*\t*\t*\t*\t*\t*\t*\t0\t0\t255", &mut w), Line::Unmapped));
        assert!(matches!(classify(&m, &rank, b"r3\t20\t0\t20\t+\t>1>2", &mut w), Line::Malformed));
        assert!(matches!(classify(&m, &rank, b"r4\t20\t0\t20\t+\tchr1:5-9\t25\t2\t22\t20\t20\t60", &mut w), Line::Malformed));
        assert!(matches!(classify(&m, &rank, b"@HD\tVN:1", &mut w), Line::Skip));
    }
}
