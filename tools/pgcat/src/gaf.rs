//! `pgcat index-gaf`: a GAF made searchable by reference position, as pgphase
//! reads it. Each alignment is placed on the reference by the nodes of its walk
//! that lie on the reference sample's paths (pggaf's rule, below), written as
//! `contig  start+1  end  <the GAF line>`, sorted, bgzipped and tabix-indexed.
//!
//! Alignments are projected in parallel batches as the input is read; sorted
//! runs that exceed the memory budget are spilled (zstd) and merged. Header
//! lines `#max_span <contig> <bases>` give the most graph bases any alignment
//! on each contig spans.

use crate::bgzf::IndexedWriter;
use gbz::support::Orientation;
use gbz::GBZ;
use rayon::prelude::*;
use std::cmp::Ordering;
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
    /// The reference beside a walk with no anchor (Unplaced::NoAnchor): the graph is searched outward from each end of
    /// the walk -- forward from its last node, and backward from its first -- for the reference node fewest graph
    /// bases away, within NEAREST_REACH. For each side, the contig, the reference position where the search joins the
    /// reference, and the bases between; `next(node, reverse)` is the handles that follow a handle.
    fn nearest(&self, walk: &[(usize, bool)], next: &dyn Fn(usize, bool) -> Vec<(usize, bool)>) -> [Option<(u32, u64, u64)>; 2] {
        let search = |from: (usize, bool)| -> Option<(u32, u64, u64)> {
            let mut heap = BinaryHeap::new();
            let mut seen: HashMap<(usize, bool), u64> = HashMap::new();
            for h in next(from.0, from.1) {
                heap.push(std::cmp::Reverse((0u64, h)));
            }
            while let Some(std::cmp::Reverse((d, (node, reverse)))) = heap.pop() {
                if seen.contains_key(&(node, reverse)) || d > NEAREST_REACH || seen.len() >= NEAREST_NODES {
                    continue;
                }
                seen.insert((node, reverse), d);
                if let Some(p) = self.placements_of(node).next() {
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
        let (Some(&first), Some(&last)) = (walk.first(), walk.last()) else { return [None, None] };
        [search((first.0, !first.1)), search(last)]
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
    malformed: u64,
}

impl Counts {
    fn add(mut self, o: Counts) -> Counts {
        self.lines += o.lines;
        self.placed += o.placed;
        self.unmapped += o.unmapped;
        self.unplaced += o.unplaced;
        self.beside += o.beside;
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

/// A placed line found in a chunk: (offset, length, contig rank, start, end).
type Found = (usize, usize, u32, u64, u64);

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

/// The input, decompressed if gzip or BGZF.
fn open_input(path: &str) -> Result<Box<dyn Read + Send>, String> {
    let mut f = BufReader::with_capacity(1 << 20, std::fs::File::open(path).map_err(|e| format!("{path}: {e}"))?);
    let gz = f.fill_buf().map(|b| b.starts_with(&[0x1f, 0x8b])).map_err(|e| format!("{path}: {e}"))?;
    Ok(if gz { Box::new(flate2::bufread::MultiGzDecoder::new(f)) } else { Box::new(f) })
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

/// Sorted runs: kept in memory, or spilled as zstd files.
struct Runs {
    dir: PathBuf,
    spilled: Vec<PathBuf>,
}

impl Runs {
    fn spill(&mut self, data: &[u8], recs: &[Rec]) -> std::io::Result<()> {
        std::fs::create_dir_all(&self.dir)?;
        let path = self.dir.join(format!("run{}.zst", self.spilled.len()));
        let mut w = zstd::stream::write::Encoder::new(BufWriter::with_capacity(1 << 22, std::fs::File::create(&path)?), 1)?;
        w.multithread(rayon::current_num_threads() as u32)?;
        for r in recs {
            w.write_all(&r.rank.to_le_bytes())?;
            w.write_all(&r.start.to_le_bytes())?;
            w.write_all(&r.end.to_le_bytes())?;
            w.write_all(&r.len.to_le_bytes())?;
            w.write_all(&data[r.off as usize..r.off as usize + r.len as usize])?;
        }
        w.finish()?.flush()?;
        self.spilled.push(path);
        Ok(())
    }
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

/// The next record of a spilled run.
struct RunReader {
    r: BufReader<zstd::stream::read::Decoder<'static, BufReader<std::fs::File>>>,
}

impl RunReader {
    fn next(&mut self, line: &mut Vec<u8>) -> std::io::Result<Option<(u32, u64, u64)>> {
        let mut h = [0u8; 24];
        match self.r.read_exact(&mut h) {
            Ok(()) => {}
            Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => return Ok(None),
            Err(e) => return Err(e),
        }
        let rank = u32::from_le_bytes(h[0..4].try_into().unwrap());
        let start = u64::from_le_bytes(h[4..12].try_into().unwrap());
        let end = u64::from_le_bytes(h[12..20].try_into().unwrap());
        let len = u32::from_le_bytes(h[20..24].try_into().unwrap()) as usize;
        line.resize(len, 0);
        self.r.read_exact(line)?;
        Ok(Some((rank, start, end)))
    }
}

struct Head {
    key: (u32, u64, u64),
    line: Vec<u8>,
    run: usize,
}

impl PartialEq for Head {
    fn eq(&self, o: &Self) -> bool {
        self.cmp(o) == Ordering::Equal
    }
}
impl Eq for Head {}
impl PartialOrd for Head {
    fn partial_cmp(&self, o: &Self) -> Option<Ordering> {
        Some(self.cmp(o))
    }
}
impl Ord for Head {
    // Reversed: BinaryHeap pops the smallest.
    fn cmp(&self, o: &Self) -> Ordering {
        (o.key, &o.line, o.run).cmp(&(self.key, &self.line, self.run))
    }
}

/// Writes one placed alignment: `contig  start+1  end  line`.
fn emit(w: &mut IndexedWriter, names: &[String], by_rank: &[u32], key: (u32, u64, u64), line: &[u8], buf: &mut Vec<u8>) -> std::io::Result<()> {
    let contig = by_rank[key.0 as usize];
    buf.clear();
    buf.extend_from_slice(names[contig as usize].as_bytes());
    buf.extend_from_slice(format!("\t{}\t{}\t", key.1 + 1, key.2).as_bytes());
    buf.extend_from_slice(line);
    buf.push(b'\n');
    w.add(buf, contig, key.1, key.2)
}

pub fn run(o: &IndexGafOptions, sample: &str) -> Result<(), String> {
    let t0 = Instant::now();
    let log = |msg: &str| eprintln!("[pgcat {:7.1}s] {msg}", t0.elapsed().as_secs_f64());
    let g: GBZ = simple_sds::serialize::load_from(&o.graph).map_err(|e| format!("{}: {e}", o.graph))?;
    let m = RefMap::build(&g, sample)?;
    // The graph stays loaded only to find the reference beside the unplaced alignments without anchors.
    let g = o.unplaced.is_some().then_some(g);
    let next = |node: usize, reverse: bool| -> Vec<(usize, bool)> {
        let o = if reverse { Orientation::Reverse } else { Orientation::Forward };
        let edges = g.as_ref().and_then(|g| g.successors(node, o));
        edges.map(|it| it.map(|(n, o)| (n, o == Orientation::Reverse)).collect()).unwrap_or_default()
    };
    log(&format!("graph loaded: {} reference paths of {sample} on {} contigs", m.path_contig.len(), m.contigs.len()));
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
        let results: Vec<(Vec<Found>, Counts, Vec<u64>, Vec<(usize, usize, Unplaced, String)>)> = pieces
            .par_iter()
            .map(|piece| {
                let base = piece.as_ptr() as usize - chunk.as_ptr() as usize;
                let mut walk = Vec::new();
                let mut placed = Vec::new();
                let mut dropped = Vec::new();
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
                            placed.push((base + p, body.len(), r, s, e));
                        }
                        Line::Skip => {}
                        Line::Unmapped => {
                            c.lines += 1;
                            c.unmapped += 1;
                        }
                        Line::Unplaced(why) => {
                            c.lines += 1;
                            c.unplaced += 1;
                            c.why[Unplaced::ALL.iter().position(|&w| w == why).unwrap()] += 1;
                            if keep_unplaced {
                                let mut tags = format!("\tpu:Z:{}", why.name());
                                if why == Unplaced::NoAnchor {
                                    let sides = m.nearest(&walk, &next);
                                    c.beside += sides.iter().any(|x| x.is_some()) as u64;
                                    for (tag, x) in ["pb", "pa"].iter().zip(sides) {
                                        if let Some((contig, at, d)) = x {
                                            tags.push_str(&format!("\t{tag}:Z:{}:{at}:{d}", m.contigs[contig as usize]));
                                        }
                                    }
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
                (placed, c, span, dropped)
            })
            .collect();
        for (placed, c, span, dropped) in results {
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
            for (at, len, r, s, e) in placed {
                recs.push(Rec { rank: r, start: s, end: e, off: data.len() as u64, len: len as u32 });
                data.extend_from_slice(&chunk[at..at + len]);
            }
        }
        if data.len() + recs.len() * std::mem::size_of::<Rec>() >= o.memory {
            sort_run(&data, &mut recs);
            runs.spill(&data, &recs).map_err(|e| format!("{}: {e}", runs.dir.display()))?;
            log(&format!("{} alignments read; sorted run {} spilled", counts.lines, runs.spilled.len()));
            data.clear();
            recs.clear();
        }
    }
    sort_run(&data, &mut recs);
    if let Some(w) = unplaced_out.take() {
        w.finish().and_then(|mut b| b.flush()).map_err(|e| format!("unplaced: {e}"))?;
        log(&format!("{} of {} alignments without anchors have the reference beside them (pb/pa tags)", counts.beside, counts.why[0]));
    }
    let why: Vec<String> = Unplaced::ALL.iter().zip(counts.why).map(|(w, n)| format!("{} {n}", w.name())).collect();
    drop(g);
    log(&format!(
        "{} alignments: {} placed, {} unmapped, {} without a unique reference interval ({}), {} malformed",
        counts.lines, counts.placed, counts.unmapped, counts.unplaced, why.join(", "), counts.malformed
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
        let mut readers = Vec::new();
        for p in &runs.spilled {
            let f = BufReader::with_capacity(1 << 20, std::fs::File::open(p).map_err(|e| format!("{}: {e}", p.display()))?);
            let d = zstd::stream::read::Decoder::with_buffer(f).map_err(|e| e.to_string())?;
            readers.push(RunReader { r: BufReader::with_capacity(1 << 20, d) });
        }
        let mut heap = BinaryHeap::new();
        for (run, rd) in readers.iter_mut().enumerate() {
            let mut line = Vec::new();
            if let Some(key) = rd.next(&mut line).map_err(|e| e.to_string())? {
                heap.push(Head { key, line, run });
            }
        }
        while let Some(mut h) = heap.pop() {
            emit(&mut w, &m.contigs, &by_rank, h.key, &h.line, &mut buf).map_err(io)?;
            if let Some(key) = readers[h.run].next(&mut h.line).map_err(|e| e.to_string())? {
                h.key = key;
                heap.push(h);
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
        assert_eq!(m.nearest(&[(7, false)], &next), [Some((0, 110, 5)), Some((0, 115, 0))]);
        // The same walk reversed: the sides swap.
        assert_eq!(m.nearest(&[(7, true)], &next), [Some((0, 115, 0)), Some((0, 110, 5))]);
        // No edges: nothing beside it.
        assert_eq!(m.nearest(&[(7, false)], &|_, _| Vec::new()), [None, None]);
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
