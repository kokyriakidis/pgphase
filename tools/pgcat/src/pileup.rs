//! `pgcat pileup` (experimental): candidate small variants from the reads' graph alignments alone -- a pileup of
//! each read's own edits (its `cs` mismatches, insertions and deletions) where the graph puts them, no alignment to
//! the reference.
//!
//! An edit on a node the reference path visits once is placed on the reference: a mismatch at its position, an
//! insertion between two reference bases next to each other, a deletion of reference bases in a row; indels are
//! left-aligned on the reference, so the same event spelled at different places in a repeat is one candidate. An
//! edit on any other node (another haplotype's or copy's sequence) is counted by node, not placed (a nested
//! variant: later). Every edit at least MIN_ALT reads carry is a candidate; each read covering its reference span
//! (one base either side for an indel) without it is REF, with another edit there OTHER.
//!
//! With --catalog (the contig's .pgsites), the reads' alleles at the catalog's snarls are their own sequences: a
//! window is a snarl's reference bases strictly between its two boundary nodes (the outermost snarls of at most
//! MAX_SITE reference bases, nested ones in larger snarls), and a read's allele there is its letters between its
//! aligned last base of one boundary node and first base of the other, whatever its walk and edits -- a walk through
//! another haplotype's allele with edits of the read's own on top is one sequence. Each distinct allele sequence is
//! aligned to the reference's and its SNVs and indels (left-aligned) are candidates, ALT the reads whose allele has
//! the variant and REF the others spanning the window; edits inside windows are not counted again as edits.
//!
//! Output: candidates.tsv's columns as pgphase writes them where they apply (CHROM POS TYPE REF ALT DP REF_COUNT
//! ALT_COUNT OTHER_COUNT FORWARD_ALT REVERSE_ALT MEAN_ALT_MAPQ) and SOURCE (edit, or snarl), sorted by position.

use crate::pgsites;
use crate::surject::{align_global_with, comp, load_contig, parse_walk, Contig, ASM_COSTS, I, M};
use gbz::GBZ;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::io::{BufRead, BufReader, Read, Write};

pub struct PileupOptions {
    pub graph: String,
    pub gaf: String,
    pub contig: String,
    pub reference_sample: String,
    pub out_contig: Option<String>,
    pub min_alt: u32,
    /// Keep only reads covering this 0-based half-open reference range, and candidates in it.
    pub region: Option<(u64, u64)>,
    /// Write each read's evidence here (name, reverse, MAPQ, reference ranges, placed edits, off-reference edits).
    pub reads_out: Option<String>,
    /// The contig's catalog (.pgsites): reads' alleles at its snarls from their sequences.
    pub catalog: Option<String>,
}

/// The most reference bases a snarl may hold to be a window of its own (larger ones give way to their children).
const MAX_SITE: u64 = 500;
/// Reference bases a window takes either side of its snarl; windows that then touch are merged (a repeat a chain of
/// snarls spans is one window: a walk through one snarl's deletion and the next one's insertion, with the bases
/// between, is the substitutions the truth spells), up to MAX_WINDOW bases.
const WINDOW_PAD: u64 = 10;
const MAX_WINDOW: u64 = 1000;

/// A placed edit: 0-based VCF position (the base before an indel), REF and ALT as VCF writes them.
#[derive(Clone, PartialEq, Eq, Hash, PartialOrd, Ord, Debug)]
pub struct Edit {
    pub pos: u64,
    pub reference: Vec<u8>,
    pub alt: Vec<u8>,
}

impl Edit {
    fn kind(&self) -> &'static str {
        match (self.reference.len(), self.alt.len()) {
            (1, 1) => "SNV",
            (r, a) if a > r => "INS",
            _ => "DEL",
        }
    }
    /// The reference span a read must cover to say REF or OTHER: the edit's bases and one either side.
    fn span(&self) -> (u64, u64) {
        (self.pos.saturating_sub(1), self.pos + self.reference.len() as u64 + 1)
    }
}

/// A deletion of `seq[s..s + len]` left-aligned, as a VCF record.
pub fn deletion(seq: &[u8], mut s: usize, len: usize) -> Option<Edit> {
    if s == 0 || s + len > seq.len() || len == 0 {
        return None;
    }
    while s > 1 && seq[s - 1].eq_ignore_ascii_case(&seq[s + len - 1]) {
        s -= 1;
    }
    Some(Edit { pos: (s - 1) as u64, reference: seq[s - 1..s + len].to_ascii_uppercase(), alt: vec![seq[s - 1].to_ascii_uppercase()] })
}

/// An insertion of `ins` after `seq[p]` left-aligned, as a VCF record.
pub fn insertion(seq: &[u8], mut p: usize, ins: &[u8]) -> Option<Edit> {
    if ins.is_empty() || p >= seq.len() {
        return None;
    }
    let mut ins: Vec<u8> = ins.to_ascii_uppercase();
    while p > 0 && seq[p].eq_ignore_ascii_case(ins.last().unwrap()) {
        ins.rotate_right(1);
        p -= 1;
    }
    let mut alt = vec![seq[p].to_ascii_uppercase()];
    alt.extend_from_slice(&ins);
    Some(Edit { pos: p as u64, reference: vec![seq[p].to_ascii_uppercase()], alt })
}

/// One read's evidence: the reference ranges its bases on reference nodes cover (merged), its placed edits, how
/// many edits it has on other nodes, its MAPQ and whether it runs against the reference.
pub struct ReadEvidence {
    pub name: String,
    /// The read's allele (its letters, in the reference's orientation) at each window it spans: (window, letters).
    pub alleles: Vec<(usize, Vec<u8>)>,
    pub ranges: Vec<(u64, u64)>,
    pub edits: Vec<Edit>,
    pub node_edits: u32,
    pub mapq: u32,
    pub reverse: bool,
}

pub fn read_evidence(g: &GBZ, c: &Contig, line: &str, windows: &[(u64, u64)]) -> Option<ReadEvidence> {
    let f: Vec<&str> = line.split('\t').collect();
    let k = if f.len() > 4 && (f[4] == "+" || f[4] == "-") { 0 } else { 3 };
    if f.len() < k + 12 || f[k + 5] == "*" {
        return None;
    }
    let pstart: usize = f[k + 7].parse().ok()?;
    let mapq: u32 = f[k + 11].parse().unwrap_or(255);
    let cs = f.iter().find_map(|t| t.strip_prefix("cs:Z:"))?;
    let walk = parse_walk(f[k + 5])?;
    let lens: Vec<usize> = walk.iter().map(|&(n, _)| g.sequence_len(n).unwrap_or(0)).collect();
    let (mut w, mut off) = (0usize, pstart);
    while w < walk.len() && off >= lens[w] {
        off -= lens[w];
        w += 1;
    }
    // The reference place of a walk base: (position, direction the walk runs along the reference).
    let place = |w: usize, off: usize| -> Option<(u64, i64)> {
        let (node, rev) = walk[w];
        let &(start, node_rev) = c.place.get(&node)?;
        let len = lens[w];
        let fo = if rev { len - 1 - off } else { off };
        let r = if node_rev { start + (len - 1 - fo) as u64 } else { start + fo as u64 };
        Some((r, if rev == node_rev { 1 } else { -1 }))
    };
    let advance = |w: &mut usize, off: &mut usize| {
        *off += 1;
        if *off == lens[*w] {
            *w += 1;
            *off = 0;
        }
    };
    let mut covered: Vec<(u64, i64)> = Vec::new();
    let mut at_q: Vec<usize> = Vec::new(); // the read offset of each covered base
    let mut read: Vec<u8> = Vec::new();
    let mut edits: Vec<Edit> = Vec::new();
    let mut node_edits = 0u32;
    let mut last: Option<(u64, i64)> = None; // the last walk base consumed, if placed
    let mut last_placed_walk = true;
    let b = cs.as_bytes();
    let mut i = 0;
    while i < b.len() {
        let op = b[i];
        i += 1;
        let st = i;
        match op {
            b':' => {
                while i < b.len() && b[i].is_ascii_digit() {
                    i += 1;
                }
                let n: usize = cs[st..i].parse().ok()?;
                for _ in 0..n {
                    if w >= walk.len() {
                        return None;
                    }
                    let p = place(w, off);
                    if let Some(pd) = p {
                        covered.push(pd);
                        at_q.push(read.len());
                    }
                    let (node, rev) = walk[w];
                    let s = g.sequence(node)?;
                    read.push(if rev { comp(s[lens[w] - 1 - off]) } else { s[off].to_ascii_uppercase() });
                    last = p;
                    last_placed_walk = p.is_some();
                    advance(&mut w, &mut off);
                }
            }
            b'*' => {
                if i + 2 > b.len() || w >= walk.len() {
                    return None;
                }
                let alt = b[i + 1].to_ascii_uppercase();
                i += 2;
                let p = place(w, off);
                match p {
                    Some((r, d)) => {
                        covered.push((r, d));
                        at_q.push(read.len());
                        let a = if d == 1 { alt } else { comp(alt) };
                        if (r as usize) < c.seq.len() && !c.seq[r as usize].eq_ignore_ascii_case(&a) {
                            edits.push(Edit { pos: r, reference: vec![c.seq[r as usize].to_ascii_uppercase()], alt: vec![a] });
                        }
                    }
                    None => node_edits += 1,
                }
                read.push(alt);
                last = p;
                last_placed_walk = p.is_some();
                advance(&mut w, &mut off);
            }
            b'+' => {
                while i < b.len() && b[i].is_ascii_alphabetic() {
                    i += 1;
                }
                let ins = &b[st..i];
                read.extend(ins.iter().map(|x| x.to_ascii_uppercase()));
                // Placed only between two reference bases next to each other along the walk.
                let next = if w < walk.len() { place(w, off) } else { None };
                match (last, next) {
                    (Some((r0, d0)), Some((r1, d1))) if d0 == d1 && r1 as i64 == r0 as i64 + d0 => {
                        let (p, seq): (u64, Vec<u8>) = if d0 == 1 { (r0, ins.to_vec()) } else { (r1, ins.iter().rev().map(|&x| comp(x)).collect()) };
                        if let Some(e) = insertion(&c.seq, p as usize, &seq) {
                            edits.push(e);
                        }
                    }
                    _ => node_edits += 1,
                }
            }
            b'-' => {
                while i < b.len() && b[i].is_ascii_alphabetic() {
                    i += 1;
                }
                let mut places = Vec::with_capacity(i - st);
                for _ in st..i {
                    if w >= walk.len() {
                        return None;
                    }
                    places.push(place(w, off));
                    advance(&mut w, &mut off);
                }
                let all: Option<Vec<(u64, i64)>> = places.iter().copied().collect();
                match all {
                    Some(v) if last_placed_walk && v.windows(2).all(|x| x[1].1 == x[0].1 && x[1].0 as i64 == x[0].0 as i64 + x[0].1) => {
                        let lo = v.iter().map(|x| x.0).min().unwrap();
                        if let Some(e) = deletion(&c.seq, lo as usize, v.len()) {
                            edits.push(e);
                        }
                        last = v.last().copied();
                    }
                    _ => {
                        node_edits += 1;
                        last = None;
                    }
                }
            }
            _ => return None,
        }
    }
    let fwd = covered.iter().filter(|x| x.1 == 1).count();
    let reverse = fwd * 2 < covered.len();
    let mut pos: Vec<u64> = covered.iter().map(|x| x.0).collect();
    pos.sort_unstable();
    let mut ranges: Vec<(u64, u64)> = Vec::new();
    for p in pos {
        match ranges.last_mut() {
            Some(r) if p <= r.1 => r.1 = r.1.max(p + 1),
            _ => ranges.push((p, p + 1)),
        }
    }
    edits.sort();
    edits.dedup();
    // Alleles at the windows the read spans: its letters between its two anchors, in the reference's orientation.
    let mut alleles = Vec::new();
    if !windows.is_empty() && !covered.is_empty() {
        let dir: i64 = if reverse { -1 } else { 1 };
        let qmap: HashMap<u64, usize> = covered.iter().zip(&at_q).filter(|(c, _)| c.1 == dir).map(|(c, &q)| (c.0, q)).collect();
        let (lo, hi) = (ranges.first().unwrap().0, ranges.last().unwrap().1);
        let first = windows.partition_point(|w| w.1 < lo);
        for (k, &(a, b)) in windows.iter().enumerate().skip(first) {
            if a > hi {
                break;
            }
            if a == 0 {
                continue;
            }
            let (Some(&qa), Some(&qb)) = (qmap.get(&(a - 1)), qmap.get(&b)) else { continue };
            let seq: Vec<u8> = if !reverse && qa < qb {
                read[qa + 1..qb].to_vec()
            } else if reverse && qb < qa {
                read[qb + 1..qa].iter().rev().map(|&x| comp(x)).collect()
            } else {
                continue;
            };
            alleles.push((k, seq));
        }
    }
    Some(ReadEvidence { name: f[k].to_string(), alleles, ranges, edits, node_edits, mapq, reverse })
}

/// The windows: the catalog's snarls' reference bases strictly between their boundary nodes, (start, end) 0-based
/// half-open, the outermost of those holding at most MAX_SITE bases, sorted.
pub(crate) fn load_windows(g: &GBZ, c: &Contig, path: &str, region: Option<(u64, u64)>) -> Result<Vec<(u64, u64)>, String> {
    let r = pgsites::Reader::open(path).map_err(|e| format!("{path}: {e}"))?;
    let (beg, end) = region.map_or((0usize, usize::MAX / 2), |(a, b)| (a.saturating_sub(MAX_SITE) as usize, (b + MAX_SITE) as usize));
    let mut all: Vec<(u64, u64)> = Vec::new();
    r.for_each_record(beg, end, |line| {
        let line = String::from_utf8_lossy(line);
        let Some(id) = line.split('\t').nth(2) else { return };
        let ends: Vec<usize> = id.split(['>', '<']).filter(|x| !x.is_empty()).filter_map(|x| x.parse().ok()).collect();
        if ends.len() != 2 {
            return;
        }
        let (Some(&(sa, _)), Some(&(sb, _))) = (c.place.get(&ends[0]), c.place.get(&ends[1])) else { return };
        let (la, lb) = (g.sequence_len(ends[0]).unwrap_or(0) as u64, g.sequence_len(ends[1]).unwrap_or(0) as u64);
        let (lo, hi) = if sa < sb { (sa + la, sb) } else { (sb + lb, sa) };
        if hi >= lo && hi - lo <= MAX_SITE {
            all.push((lo, hi));
        }
    })
    .map_err(|e| e.to_string())?;
    all.sort_by(|x, y| x.0.cmp(&y.0).then(y.1.cmp(&x.1)));
    all.dedup();
    let mut outer: Vec<(u64, u64)> = Vec::new();
    for w in all {
        if let Some(&last) = outer.last() {
            if w.0 >= last.0 && w.1 <= last.1 && !(w.0 == w.1 && w.0 == last.1) {
                continue;
            }
        }
        outer.push(w);
    }
    let len = c.seq.len() as u64;
    let mut out: Vec<(u64, u64)> = Vec::new();
    for (lo, hi) in outer {
        let (lo, hi) = (lo.saturating_sub(WINDOW_PAD).max(1), (hi + WINDOW_PAD).min(len.saturating_sub(1)));
        match out.last_mut() {
            Some(last) if lo <= last.1 && hi.max(last.1) - last.0 <= MAX_WINDOW => last.1 = last.1.max(hi),
            _ => out.push((lo, hi)),
        }
    }
    Ok(out)
}

/// The SNVs and indels (left-aligned) an allele makes of the reference's letters `seq[lo..hi]`.
fn decompose(seq: &[u8], lo: u64, hi: u64, allele: &[u8]) -> Vec<Edit> {
    let reference = &seq[lo as usize..hi as usize];
    let mut out = Vec::new();
    let (mut i, mut j) = (0usize, 0usize); // allele, reference
    for (op, n) in align_global_with(allele, reference, ASM_COSTS) {
        let n = n as usize;
        match op {
            M => {
                for t in 0..n {
                    if !allele[i + t].eq_ignore_ascii_case(&reference[j + t]) {
                        let p = lo as usize + j + t;
                        out.push(Edit { pos: p as u64, reference: vec![seq[p].to_ascii_uppercase()], alt: vec![allele[i + t].to_ascii_uppercase()] });
                    }
                }
                i += n;
                j += n;
            }
            I => {
                if let Some(e) = insertion(seq, lo as usize + j - 1, &allele[i..i + n]) {
                    out.push(e);
                }
                i += n;
            }
            _ => {
                if let Some(e) = deletion(seq, lo as usize + j, n) {
                    out.push(e);
                }
                j += n;
            }
        }
    }
    out
}

#[derive(Default, Clone)]
struct Tally {
    alt: u32,
    fwd_alt: u32,
    rev_alt: u32,
    mapq_alt: u64,
    reference: u32,
    other: u32,
    snarl: bool,
}

pub fn run(o: &PileupOptions) -> Result<(), String> {
    let g: GBZ = simple_sds::serialize::load_from(&o.graph).map_err(|e| format!("{}: {e}", o.graph))?;
    let (contig, chrom) = load_contig(&g, &o.reference_sample, &o.contig)?;
    let out_contig = o.out_contig.clone().unwrap_or(chrom);
    let file = std::fs::File::open(&o.gaf).map_err(|e| format!("{}: {e}", o.gaf))?;
    let input: Box<dyn Read> = if o.gaf.ends_with(".gz") { Box::new(flate2::read::MultiGzDecoder::new(file)) } else { Box::new(file) };
    let mut lines = BufReader::with_capacity(1 << 22, input).lines();
    let windows = match &o.catalog {
        Some(p) => load_windows(&g, &contig, p, o.region)?,
        None => Vec::new(),
    };
    eprintln!("pgcat pileup: {} snarl windows", windows.len());
    let mut reads: Vec<ReadEvidence> = Vec::new();
    loop {
        let batch: Vec<String> = lines.by_ref().take(20_000).filter_map(|l| l.ok()).filter(|l| !l.starts_with('#') && !l.is_empty()).collect();
        if batch.is_empty() {
            break;
        }
        let ev: Vec<Option<ReadEvidence>> = batch.par_iter().map(|l| read_evidence(&g, &contig, l, &windows)).collect();
        reads.extend(ev.into_iter().flatten().filter(|r| match o.region {
            Some((a, b)) => r.ranges.iter().any(|&(s, e)| s < b && e > a),
            None => true,
        }));
    }
    if let Some(path) = &o.reads_out {
        let mut w = std::io::BufWriter::new(std::fs::File::create(path).map_err(|e| format!("{path}: {e}"))?);
        for r in &reads {
            let ranges: Vec<String> = r.ranges.iter().map(|x| format!("{}-{}", x.0, x.1)).collect();
            let edits: Vec<String> = r.edits.iter().map(|e| format!("{}:{}>{}", e.pos + 1, String::from_utf8_lossy(&e.reference), String::from_utf8_lossy(&e.alt))).collect();
            writeln!(w, "{}\t{}\t{}\t{}\t{}\t{}", r.name, r.reverse as u8, r.mapq, ranges.join(","), edits.join(","), r.node_edits).map_err(|e| e.to_string())?;
        }
    }
    // Candidates: edits enough reads carry, outside the windows; and the variants of the reads' window alleles.
    let in_window = |e: &Edit| {
        let i = windows.partition_point(|w| w.1 + 1 < e.pos);
        windows[i..].iter().take_while(|w| w.0 <= e.pos + 1).any(|w| e.pos + 1 >= w.0 && e.pos <= w.1)
    };
    let mut tally: HashMap<Edit, Tally> = HashMap::new();
    let mut by_window: Vec<HashMap<Vec<u8>, (u32, u32, u32, u64)>> = vec![HashMap::new(); windows.len()];
    for r in &reads {
        for (k, a) in &r.alleles {
            let t = by_window[*k].entry(a.clone()).or_default();
            t.0 += 1;
            if r.reverse { t.2 += 1 } else { t.1 += 1 }
            t.3 += r.mapq as u64;
        }
    }
    for (k, alleles) in by_window.iter().enumerate() {
        let (lo, hi) = windows[k];
        let total: u32 = alleles.values().map(|v| v.0).sum();
        let mut found: HashMap<Edit, (u32, u32, u32, u64)> = HashMap::new();
        for (a, v) in alleles {
            if a.as_slice() == &contig.seq[lo as usize..hi as usize] {
                continue;
            }
            let mut vs = decompose(&contig.seq, lo, hi, a);
            vs.sort();
            vs.dedup();
            for e in vs {
                let f = found.entry(e).or_default();
                f.0 += v.0;
                f.1 += v.1;
                f.2 += v.2;
                f.3 += v.3;
            }
        }
        for (e, (alt, fwd, rev, mq)) in found {
            let t = tally.entry(e).or_default();
            t.alt += alt;
            t.fwd_alt += fwd;
            t.rev_alt += rev;
            t.mapq_alt += mq;
            t.reference += total - alt;
            t.snarl = true;
        }
    }
    for r in &reads {
        for e in r.edits.iter().filter(|e| !in_window(e)) {
            let t = tally.entry(e.clone()).or_default();
            t.alt += 1;
            if r.reverse { t.rev_alt += 1 } else { t.fwd_alt += 1 }
            t.mapq_alt += r.mapq as u64;
        }
    }
    tally.retain(|e, t| t.alt >= o.min_alt && o.region.map_or(true, |(a, b)| e.pos >= a && e.pos < b));
    let mut cands: Vec<Edit> = tally.keys().cloned().collect();
    cands.sort();
    let edit_cands: Vec<Edit> = cands.iter().filter(|e| !tally[*e].snarl).cloned().collect();
    let cands_all = cands;
    let cands = edit_cands;
    let starts: Vec<u64> = cands.iter().map(|e| e.span().0).collect();
    let widest = cands.iter().map(|e| e.span().1 - e.span().0).max().unwrap_or(0);
    // REF and OTHER: the reads covering each candidate's span without it.
    let counts: Vec<(usize, bool)> = reads
        .par_iter()
        .flat_map_iter(|r| {
            let has: HashSet<&Edit> = r.edits.iter().collect();
            let at: HashSet<u64> = r.edits.iter().map(|e| e.pos).collect();
            let mut out = Vec::new();
            for &(a, b) in &r.ranges {
                let lo = starts.partition_point(|&s| s + widest < a);
                for (k, e) in cands.iter().enumerate().skip(lo) {
                    let (s0, s1) = e.span();
                    if s0 >= b {
                        break;
                    }
                    if s0 < a || s1 > b || has.contains(e) {
                        continue;
                    }
                    out.push((k, at.contains(&e.pos)));
                }
            }
            out.into_iter()
        })
        .collect();
    for (k, other) in counts {
        let t = tally.get_mut(&cands[k]).unwrap();
        if other { t.other += 1 } else { t.reference += 1 }
    }
    let stdout = std::io::stdout();
    let mut w = std::io::BufWriter::with_capacity(1 << 22, stdout.lock());
    let io = |e: std::io::Error| e.to_string();
    writeln!(w, "CHROM\tPOS\tTYPE\tREF\tALT\tDP\tREF_COUNT\tALT_COUNT\tOTHER_COUNT\tFORWARD_ALT\tREVERSE_ALT\tMEAN_ALT_MAPQ\tSOURCE").map_err(io)?;
    for e in &cands_all {
        let t = &tally[e];
        writeln!(
            w,
            "{out_contig}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.1}\t{}",
            e.pos + 1,
            e.kind(),
            String::from_utf8_lossy(&e.reference),
            String::from_utf8_lossy(&e.alt),
            t.alt + t.reference + t.other,
            t.reference,
            t.alt,
            t.other,
            t.fwd_alt,
            t.rev_alt,
            t.mapq_alt as f64 / t.alt.max(1) as f64,
            if t.snarl { "snarl" } else { "edit" }
        )
        .map_err(io)?;
    }
    w.flush().map_err(io)?;
    let node_edits: u64 = reads.iter().map(|r| r.node_edits as u64).sum();
    let snarl = cands_all.len() - cands.len();
    eprintln!("pgcat pileup: {} reads, {} candidates (>= {} reads; {} from snarl alleles, {} edits), {} read edits off reference nodes (not placed)", reads.len(), cands_all.len(), o.min_alt, snarl, cands.len(), node_edits);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn window_alleles_decomposed() {
        let seq = b"ACGTACGTAC";
        // The window GTAC (2..6): a substitution, a deletion (left-aligned), an insertion.
        let e = decompose(seq, 2, 6, b"GTTC");
        assert_eq!(e, vec![Edit { pos: 4, reference: b"A".to_vec(), alt: b"T".to_vec() }]);
        let e = decompose(seq, 2, 6, b"GTC");
        assert_eq!(e, vec![Edit { pos: 3, reference: b"TA".to_vec(), alt: b"T".to_vec() }]);
        let e = decompose(seq, 2, 6, b"GTAAC");
        assert_eq!(e, vec![Edit { pos: 3, reference: b"T".to_vec(), alt: b"TA".to_vec() }]);
        // Three substitutions in a repeat stay substitutions (not a deletion and an insertion).
        let seq = b"CAGAGAGAAAGAGAGAGC";
        let mut allele = seq[1..17].to_vec();
        allele[2] = b'G';
        allele[6] = b'G';
        allele[10] = b'G';
        let e = decompose(seq, 1, 17, &allele);
        assert_eq!(e.len(), 3, "{e:?}");
        assert!(e.iter().all(|x| x.reference.len() == 1 && x.alt.len() == 1));
    }

    #[test]
    fn indels_left_aligned() {
        let seq = b"GCAAAAT";
        // Deleting any one A is the same event: after the C.
        let d = deletion(seq, 4, 1).unwrap();
        assert_eq!((d.pos, d.reference.as_slice(), d.alt.as_slice()), (1, &b"CA"[..], &b"C"[..]));
        // Inserting an A after any A: after the C.
        let i = insertion(seq, 4, b"A").unwrap();
        assert_eq!((i.pos, i.reference.as_slice(), i.alt.as_slice()), (1, &b"C"[..], &b"CA"[..]));
        // A dinucleotide repeat insertion moves left by whole units' rotation.
        let seq = b"GTACACACG";
        let i = insertion(seq, 6, b"CA").unwrap();
        assert_eq!((i.pos, i.alt.as_slice()), (1, &b"TAC"[..]));
    }
}
