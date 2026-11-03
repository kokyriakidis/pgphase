//! Rescue (`index-gaf --rescue`): a reference position for an alignment whose walk cannot give one, a read whose
//! aligned nodes all lie off the reference (no_anchor). The graph holds such sequence apart from the reference's --
//! another assembly's copy of an inversion or a segmental duplication, never aligned to it -- so the reference nodes
//! beside the walk (RefMap::nearest) say only roughly where the reference holds the same letters: of HG002's reads
//! that minimap2 places confidently, 13% lie between those two nodes, 86% within 300 kb of them and 95% within 1 Mb.
//!
//! So the read's letters (its walk spelled with its cs edits, `spell`) are aligned to the reference within REACH of
//! the nodes beside the walk. One `Index` holds the reference's canonical minimizers (K, W) over every such stretch of
//! a run's reads, less those occurring more than MAX_OCCURRENCES times (repeats, as minimap2 drops them); the read's
//! minimizers are looked up in it, and the hits in its own stretch kept. A read whose best chain's blocks hold less
//! than MIN_CHAINED of its letters is not aligned at all: a read of the same sequence (99% identity or better) is
//! mostly blocks, one of a paralog or a satellite (its minimizers repeats) gaps, and aligning those gaps was nearly
//! all the cost. Hits on one diagonal are
//! joined into blocks; the blocks of each strand are chained as surject chains a read's
//! placed blocks (surject::chain), and the best chain is aligned as surject aligns those (surject::linearize). The
//! best chain again, without the blocks within a read's length of the first one, scores the alignment's uniqueness
//! (its mapping quality). An alignment is kept if it holds MIN_COVERAGE of the read at MIN_IDENTITY or better.
//!
//! Its uniqueness is within the stretch searched, not the genome: a copy further away is not seen. A sample's copy the
//! reference lacks (an extra copy of a copy-number variable duplication: HG002's at chr17:38.8 Mb, 17q12) is such a
//! read: it aligns at 60 to the reference's one copy, whose own reads are all there already. So a read is placed only
//! where the placed alignments cover the reference below its contig's median depth (`Depth`): a locus missing reads
//! (one haplotype's diverged from the reference's, so its reads had no anchor) takes them back; one at full depth
//! takes none. Placed at full depth, the 17q12 reads (1.15-1.42x of the median before) wrote 36 false hets there and
//! made 29 reads discordant; the reads that brought calls back (chr16:28.6-29.0, chr18:47.0-47.3, chr19:57.7-58.0 Mb)
//! lay at 0.0-0.7x. Nor is a read placed where the reads placed together would bring the reference to three copies'
//! depth (`MAX_RESCUED_DEPTH`): a locus one haplotype's reads had missed, filled back past two copies with the reads
//! of an extra copy beside them -- 16p13.11 (chr16:14.9-15.1 Mb), 0.15x of the median before and 1.28x after, wrote
//! 56 false hets there while every 100 kb the reads brought calls back in stayed at 0.4-1.21x.

use crate::surject::{chain, comp, join_cost, D, I, M, S};
use rayon::prelude::*;

/// How far either side of the reference nodes beside a walk the read is looked for.
pub const REACH: u64 = 1_000_000;
/// Minimizer length and window (k-mers per window).
const K: usize = 19;
const W: usize = 10;
/// A minimizer found more often than this in the index is a repeat's and is left out of it.
const MAX_OCCURRENCES: usize = 32;
/// The least share of the read's letters the best chain's blocks must hold for the read to be aligned.
const MIN_CHAINED: f64 = 0.5;
/// The least share of the read an alignment must hold, and the least identity (matching columns over all aligned
/// columns, gaps included), to be kept.
const MIN_COVERAGE: f64 = 0.8;
const MIN_IDENTITY: f64 = 0.99;

/// Reference bases per bin of `Depth`.
const DEPTH_BIN: u64 = 1_000;
/// The share of its contig's median depth a placement's reference may be covered at, at most, to take a read.
const MAX_DEPTH: f64 = 1.0;
/// The share of its contig's median depth at which a locus holds three copies' reads rather than two (half way from
/// two copies, 1.0, to three, 1.5): the reads --rescue places may bring a placement's reference below it, not to it.
pub const MAX_RESCUED_DEPTH: f64 = 1.25;

/// How deep the placed alignments cover each contig, in DEPTH_BIN bins: the bases covered per bin.
#[derive(Default, Clone)]
pub struct Depth {
    bins: Vec<Vec<u64>>,
    median: Vec<f64>,
}

impl Depth {
    /// An alignment placed on `contig` over [start, end).
    pub fn add(&mut self, contig: usize, start: u64, end: u64) {
        if contig >= self.bins.len() {
            self.bins.resize_with(contig + 1, Vec::new);
        }
        let bins = &mut self.bins[contig];
        let last = (end.max(1) - 1) / DEPTH_BIN;
        if bins.len() <= last as usize {
            bins.resize(last as usize + 1, 0);
        }
        for b in start / DEPTH_BIN..=last {
            let (lo, hi) = (start.max(b * DEPTH_BIN), end.min((b + 1) * DEPTH_BIN));
            bins[b as usize] += hi.saturating_sub(lo);
        }
    }

    /// Each contig's median depth over its covered bins, once every alignment is added.
    pub fn finish(&mut self) {
        self.median = self
            .bins
            .iter()
            .map(|bins| {
                let mut v: Vec<u64> = bins.iter().copied().filter(|&x| x > 0).collect();
                if v.is_empty() {
                    return 0.0;
                }
                let mid = v.len() / 2;
                v.select_nth_unstable(mid);
                v[mid] as f64 / DEPTH_BIN as f64
            })
            .collect();
    }

    /// Whether [start, end) of `contig` is covered at MAX_DEPTH of its median or deeper.
    pub fn full(&self, contig: usize, start: u64, end: u64) -> bool {
        self.share(contig, start, end) >= MAX_DEPTH
    }

    /// How deep [start, end) of `contig` is covered, as a share of the contig's median (0 where it has none). The
    /// median stays the one `finish` found: alignments added after it count in the depth, not in the median.
    pub fn share(&self, contig: usize, start: u64, end: u64) -> f64 {
        let (Some(bins), Some(&median)) = (self.bins.get(contig), self.median.get(contig)) else { return 0.0 };
        if median == 0.0 || end <= start {
            return 0.0;
        }
        let covered: u64 = (start / DEPTH_BIN..=(end - 1) / DEPTH_BIN)
            .map(|b| {
                let (lo, hi) = (start.max(b * DEPTH_BIN), end.min((b + 1) * DEPTH_BIN));
                bins.get(b as usize).map_or(0, |&x| x * (hi - lo) / DEPTH_BIN)
            })
            .sum();
        covered as f64 / (end - start) as f64 / median
    }
}

/// A read aligned to the reference by `rescue`.
#[derive(Debug, PartialEq)]
pub struct Rescue {
    /// 0-based reference position of the first aligned base.
    pub pos: u64,
    /// The read's letters (`spell`'s, in the walk's orientation) reverse-complemented to align.
    pub reverse: bool,
    pub cigar: Vec<(u8, u32)>,
    pub mapq: u8,
}

impl Rescue {
    /// The reference bases the alignment covers.
    pub fn span(&self) -> u64 {
        self.cigar.iter().filter(|o| o.0 == M || o.0 == D).map(|o| o.1 as u64).sum()
    }

    /// `CONTIG:POS:STRAND:CIGAR`, POS 0-based and STRAND + or - (the rs:Z: tag).
    pub fn tag(&self, contig: &str) -> String {
        let cigar: String = self.cigar.iter().map(|&(op, n)| format!("{n}{}", op as char)).collect();
        format!("{contig}:{}:{}:{cigar}", self.pos, if self.reverse { '-' } else { '+' })
    }
}

/// The rs:Z: tag's fields: contig, 0-based position, reverse, CIGAR.
pub fn parse_tag(tag: &str) -> Option<(&str, u64, bool, Vec<(u8, u32)>)> {
    let mut f = tag.rsplitn(4, ':');
    let (cigar, strand, pos, contig) = (f.next()?, f.next()?, f.next()?, f.next()?);
    let mut ops = Vec::new();
    let mut n = 0u32;
    for c in cigar.bytes() {
        if c.is_ascii_digit() {
            n = n.checked_mul(10)?.checked_add((c - b'0') as u32)?;
        } else {
            if n == 0 || ![M, I, D, S].contains(&c) {
                return None;
            }
            ops.push((c, n));
            n = 0;
        }
    }
    let reverse = match strand {
        "+" => false,
        "-" => true,
        _ => return None,
    };
    (n == 0 && !ops.is_empty()).then_some((contig, pos.parse().ok()?, reverse, ops))
}

/// The read's letters as its graph alignment spells them, in the walk's orientation: the walk's letters from
/// `pstart` with the cs edits applied (a substitution's and an insertion's letters the read's; a deletion's skipped).
/// `seq(node)` is a node's letters. None if the cs is malformed or runs past the walk.
pub fn spell(walk: &[(usize, bool)], seq: &dyn Fn(usize) -> Option<Vec<u8>>, pstart: usize, cs: &str) -> Option<Vec<u8>> {
    let mut path: Vec<u8> = Vec::new();
    for &(node, reverse) in walk {
        let s = seq(node)?;
        if reverse {
            path.extend(s.iter().rev().map(|&c| comp(c)));
        } else {
            path.extend(s.iter().map(|c| c.to_ascii_uppercase()));
        }
    }
    let b = cs.as_bytes();
    let (mut i, mut at) = (0usize, pstart);
    let mut read = Vec::new();
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
                read.extend_from_slice(path.get(at..at + n)?);
                at += n;
            }
            b'*' => {
                if i + 2 > b.len() || at >= path.len() {
                    return None;
                }
                read.push(b[i + 1].to_ascii_uppercase());
                i += 2;
                at += 1;
            }
            b'+' => {
                while i < b.len() && b[i].is_ascii_alphabetic() {
                    i += 1;
                }
                read.extend(b[st..i].iter().map(|x| x.to_ascii_uppercase()));
            }
            b'-' => {
                while i < b.len() && b[i].is_ascii_alphabetic() {
                    i += 1;
                }
                at += i - st;
                if at > path.len() {
                    return None;
                }
            }
            _ => return None,
        }
    }
    Some(read)
}

/// minimap2's invertible integer hash.
fn hash64(key: u64, mask: u64) -> u64 {
    let mut k = (!key).wrapping_add(key << 21) & mask;
    k ^= k >> 24;
    k = (k.wrapping_add(k << 3)).wrapping_add(k << 8) & mask;
    k ^= k >> 14;
    k = (k.wrapping_add(k << 2)).wrapping_add(k << 4) & mask;
    k ^= k >> 28;
    k.wrapping_add(k << 31) & mask
}

fn code(c: u8) -> Option<u64> {
    match c {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

/// Every k-mer of `s` without an N, canonical: (hash, position of its first letter, whether the reverse complement
/// is the one hashed), calling `f` on each.
fn canonical_kmers(s: &[u8], mut f: impl FnMut(u64, usize, bool)) {
    let mask = (1u64 << (2 * K)) - 1;
    let shift = 2 * (K as u64 - 1);
    let (mut fwd, mut rev, mut len) = (0u64, 0u64, 0usize);
    for (i, &c) in s.iter().enumerate() {
        let Some(x) = code(c) else {
            len = 0;
            continue;
        };
        fwd = ((fwd << 2) | x) & mask;
        rev = (rev >> 2) | ((3 - x) << shift);
        len += 1;
        if len >= K && fwd != rev {
            let reverse = rev < fwd;
            f(hash64(if reverse { rev } else { fwd }, mask), i + 1 - K, reverse);
        }
    }
}

/// The canonical minimizers of `s`, the smallest hash of every W consecutive k-mers (one window if there are
/// fewer): (hash, position, whether its reverse complement was hashed), in position order.
fn minimizers(s: &[u8]) -> Vec<(u64, u32, bool)> {
    let mut kmers: Vec<(u64, usize, bool)> = Vec::new();
    canonical_kmers(s, |h, p, r| kmers.push((h, p, r)));
    let mut out = Vec::new();
    let mut last = usize::MAX;
    let windows = if kmers.len() >= W { kmers.len() - W + 1 } else { kmers.len().min(1) };
    for w in 0..windows {
        let end = (w + W).min(kmers.len());
        let k = (w..end).min_by_key(|&k| (kmers[k].0, kmers[k].1)).unwrap();
        if k != last {
            out.push((kmers[k].0, kmers[k].1 as u32, kmers[k].2));
            last = k;
        }
    }
    out
}

/// The reference's minimizers over the stretches a run's reads are looked for in, repeats left out: (hash, contig,
/// position with the reverse flag in its top bit), sorted.
pub struct Index {
    entries: Vec<(u64, u32, u32)>,
    /// Each minimizer's run of entries: its first and how many.
    runs: crate::mixhash::MixMap<u64, (u32, u32)>,
}

const REVERSE_BIT: u32 = 1 << 31;

impl Index {
    /// The index of `letters` (each contig's) over `stretches` (contig, start, end; overlapping ones merged).
    pub fn build(letters: &[Vec<u8>], stretches: &[(u32, u64, u64)]) -> Index {
        let mut merged: Vec<(u32, u64, u64)> = Vec::new();
        let mut sorted = stretches.to_vec();
        sorted.sort_unstable();
        for (c, lo, hi) in sorted {
            let hi = hi.min(letters[c as usize].len() as u64);
            match merged.last_mut() {
                Some(last) if last.0 == c && lo <= last.2 => last.2 = last.2.max(hi),
                _ if lo < hi => merged.push((c, lo, hi)),
                _ => {}
            }
        }
        let mut entries: Vec<(u64, u32, u32)> = merged
            .par_iter()
            .flat_map_iter(|&(c, lo, hi)| {
                minimizers(&letters[c as usize][lo as usize..hi as usize])
                    .into_iter()
                    .map(move |(h, p, r)| (h, c, (lo as u32 + p) | if r { REVERSE_BIT } else { 0 }))
            })
            .collect();
        entries.par_sort_unstable();
        // Repeats out.
        let mut kept = Vec::with_capacity(entries.len());
        let distinct = entries.par_windows(2).filter(|w| w[0].0 != w[1].0).count() + 1;
        let mut runs = crate::mixhash::MixMap::with_capacity_and_hasher(distinct, Default::default());
        let mut i = 0;
        while i < entries.len() {
            let j = i + entries[i..].partition_point(|e| e.0 == entries[i].0);
            if j - i <= MAX_OCCURRENCES {
                runs.insert(entries[i].0, (kept.len() as u32, (j - i) as u32));
                kept.extend_from_slice(&entries[i..j]);
            }
            i = j;
        }
        Index { entries: kept, runs }
    }

    /// Minimizers in the index.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    fn find(&self, h: u64) -> &[(u64, u32, u32)] {
        match self.runs.get(&h) {
            Some(&(at, n)) => &self.entries[at as usize..(at + n) as usize],
            None => &[],
        }
    }
}

/// Where a read whose walk has no anchor is looked for (contig, start, end): within REACH of the reference beside its
/// walk (`sides`, RefMap::nearest's) -- both sides when they are on one contig, else the nearer.
pub fn stretch(sides: &[Option<(u32, u64, u64)>; 2]) -> Option<(u32, u64, u64)> {
    let (contig, lo, hi) = match sides {
        [Some(a), Some(b)] if a.0 == b.0 => (a.0, a.1.min(b.1), a.1.max(b.1)),
        [Some(a), Some(b)] => if a.2 <= b.2 { (a.0, a.1, a.1) } else { (b.0, b.1, b.1) },
        [Some(a), None] | [None, Some(a)] => (a.0, a.1, a.1),
        [None, None] => return None,
    };
    Some((contig, lo.saturating_sub(REACH), hi + REACH))
}

/// A chain's score as surject::chain scores it: its bases less the join costs between consecutive blocks.
fn chain_score(blocks: &[(usize, u64, usize)]) -> f64 {
    let bases: usize = blocks.iter().map(|b| b.2).sum();
    let joins: f64 = blocks
        .windows(2)
        .map(|w| {
            let (dq, dr) = ((w[1].0 - (w[0].0 + w[0].2)) as u64, w[1].1 - (w[0].1 + w[0].2 as u64));
            join_cost(dq.abs_diff(dr))
        })
        .sum();
    bases as f64 - joins
}

/// Hits on one diagonal (read offset, reference position, sorted by diagonal then read offset) joined where they
/// overlap or touch, then sorted by read offset: surject::chain's blocks.
fn blocks(mut hits: Vec<(usize, u64)>) -> Vec<(usize, u64, usize)> {
    hits.sort_unstable_by_key(|&(q, r)| (r as i64 - q as i64, q));
    let mut out: Vec<(usize, u64, usize)> = Vec::new();
    for (q, r) in hits {
        if let Some(last) = out.last_mut() {
            if last.1 as i64 - last.0 as i64 == r as i64 - q as i64 && q <= last.0 + last.2 {
                last.2 = last.2.max(q + K - last.0);
                continue;
            }
        }
        out.push((q, r, K));
    }
    out.sort_unstable();
    out
}

/// `letters` (the read in the chain `c`'s orientation) aligned whole by karan, at least cost, inside the reference
/// the chain spans, widened by a tenth of the read and 100 letters either side, on the chain's diagonals and as many
/// more either side as a fiftieth of the read and 100; then each end soft-clipped at the column where it scores lowest,
/// if below -CLIP (BWA-MEM's rule, as surject clips). The reference position of the first aligned letter, and the
/// CIGAR.
/// BWA-MEM's clipping penalty: a read end is clipped only where clipping it gains more.
#[cfg(feature = "karan")]
const CLIP: i64 = 5;

#[cfg(feature = "karan")]
fn place_read(letters: &[u8], c: &[(usize, u64, usize)], reference: &[u8]) -> Option<(u64, Vec<(u8, u32)>)> {
    let n = letters.len() as u64;
    let pad = n / 10 + 100;
    let (first, last) = (c[0], c[c.len() - 1]);
    let lo = first.1.saturating_sub(first.0 as u64 + pad);
    let hi = (last.1 + (n - last.0 as u64) + pad).min(reference.len() as u64);
    if hi <= lo {
        return None;
    }
    let diagonal = |b: &(usize, u64, usize)| b.1 as i64 - lo as i64 - b.0 as i64;
    let slack = (n / 50 + 100) as i64;
    let band = (c.iter().map(diagonal).min()? - slack, c.iter().map(diagonal).max()? + slack);
    let p = crate::karan::place(&reference[lo as usize..hi as usize], letters, band)?;
    // Each read end clipped at the column where the end scores lowest, if below -CLIP (BWA-MEM's local scoring:
    // match 1, mismatch 4, a gap of k letters 6 + 2k or 26 + k, the lesser).
    let mut cols: Vec<(u8, u32)> = Vec::new();
    for &(op, len) in &p.cigar {
        if op == b'=' || op == b'X' {
            cols.extend(std::iter::repeat_n((op, 1), len as usize));
        } else {
            cols.push((op, len));
        }
    }
    let gain = |&(op, len): &(u8, u32)| -> i64 {
        match op {
            b'=' => 1,
            b'X' => -4,
            _ => -(6 + 2 * len as i64).min(26 + len as i64),
        }
    };
    let cut = |cols: &mut dyn Iterator<Item = &(u8, u32)>| -> usize {
        let (mut sum, mut low, mut at) = (0i64, -CLIP, 0usize);
        for (k, c) in cols.enumerate() {
            sum += gain(c);
            if sum < low {
                (low, at) = (sum, k + 1);
            }
        }
        at
    };
    let head = cut(&mut cols.iter());
    let tail = cut(&mut cols.iter().rev().take(cols.len() - head));
    // The read's and the reference's letters a run of columns holds.
    let read_letters = |cs: &[(u8, u32)]| cs.iter().filter(|c| c.0 != b'D').map(|c| c.1 as u64).sum::<u64>();
    let reference_letters = |cs: &[(u8, u32)]| cs.iter().filter(|c| c.0 != b'I').map(|c| c.1 as u64).sum::<u64>();
    let mut cigar: Vec<(u8, u32)> = Vec::new();
    let mut push = |op: u8, len: u64| match cigar.last_mut() {
        Some(last) if last.0 == op => last.1 += len as u32,
        _ if len > 0 => cigar.push((op, len as u32)),
        _ => {}
    };
    push(S, p.query_start as u64 + read_letters(&cols[..head]));
    for &(op, len) in &cols[head..cols.len() - tail] {
        push(if op == b'=' || op == b'X' { M } else { op }, len as u64);
    }
    push(S, (letters.len() - p.query_end) as u64 + read_letters(&cols[cols.len() - tail..]));
    Some((lo + p.reference_start as u64 + reference_letters(&cols[..head]), cigar))
}

/// Without karan the rescue aligns nothing (--rescue is refused there).
#[cfg(not(feature = "karan"))]
fn place_read(_letters: &[u8], _c: &[(usize, u64, usize)], _reference: &[u8]) -> Option<(u64, Vec<(u8, u32)>)> {
    None
}

/// `read` (`spell`'s letters) aligned to contig `contig` (`reference` its letters) between `lo` and `hi`, by the
/// hits of its minimizers in `index`; None if no alignment holds MIN_COVERAGE of it at MIN_IDENTITY.
pub fn rescue(read: &[u8], index: &Index, contig: u32, reference: &[u8], lo: u64, hi: u64) -> Option<Rescue> {
    let n = read.len();
    if n < K + W || lo >= hi {
        return None;
    }
    // Hits per strand: (read offset in that strand's letters, reference position).
    let mut hits: [Vec<(usize, u64)>; 2] = [Vec::new(), Vec::new()];
    for (h, q, q_rev) in minimizers(read) {
        for &(_, c, packed) in index.find(h) {
            let (pos, r_rev) = ((packed & !REVERSE_BIT) as u64, packed & REVERSE_BIT != 0);
            if c != contig || pos < lo || pos >= hi {
                continue;
            }
            let reverse = q_rev != r_rev;
            let q = if reverse { n - (q as usize + K) } else { q as usize };
            hits[reverse as usize].push((q, pos));
        }
    }
    // Each strand's best chain, its score, and the score of its best chain without the blocks within a read's length
    // of that one.
    let mut found: Vec<(f64, f64, bool, Vec<(usize, u64, usize)>)> = Vec::new();
    for (strand, h) in hits.into_iter().enumerate() {
        let b = blocks(h);
        if b.is_empty() {
            continue;
        }
        let c = chain(&b);
        let s = chain_score(&c);
        let (c_lo, c_hi) = (c[0].1.saturating_sub(n as u64), c[c.len() - 1].1 + (2 * n) as u64);
        let rest: Vec<(usize, u64, usize)> = b.into_iter().filter(|x| x.1 < c_lo || x.1 > c_hi).collect();
        let s_rest = if rest.is_empty() { 0.0 } else { chain_score(&chain(&rest)) };
        found.push((s, s_rest, strand == 1, c));
    }
    found.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
    // The runner-up: the best strand's next chain, or the other strand's best.
    let second = found.iter().skip(1).map(|x| x.0).fold(found.first()?.1, f64::max);
    let (s, _, reverse, c) = found.swap_remove(0);
    if (c.iter().map(|b| b.2).sum::<usize>() as f64) < MIN_CHAINED * n as f64 {
        return None;
    }
    let letters: Vec<u8> = if reverse { read.iter().rev().map(|&x| comp(x)).collect() } else { read.to_vec() };
    let (pos, cigar) = place_read(&letters, &c, reference)?;
    // Identity over the aligned columns, coverage of the read.
    let (mut q, mut r) = (0usize, pos as usize);
    let (mut matches, mut columns, mut clipped) = (0usize, 0usize, 0usize);
    for &(op, len) in &cigar {
        let len = len as usize;
        match op {
            M => {
                matches += (0..len).filter(|&k| letters[q + k].eq_ignore_ascii_case(&reference[r + k])).count();
                columns += len;
                q += len;
                r += len;
            }
            I => {
                columns += len;
                q += len;
            }
            D => {
                columns += len;
                r += len;
            }
            _ => {
                clipped += len;
                q += len;
            }
        }
    }
    let coverage = (n - clipped) as f64 / n as f64;
    let identity = if columns == 0 { 0.0 } else { matches as f64 / columns as f64 };
    if coverage < MIN_COVERAGE || identity < MIN_IDENTITY {
        return None;
    }
    let mapq = (60.0 * (1.0 - second / s)).round().clamp(0.0, 60.0) as u8;
    Some(Rescue { pos, reverse, cigar, mapq })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A reproducible random sequence.
    #[cfg(feature = "karan")]
    fn random(n: usize, seed: u64) -> Vec<u8> {
        let mut x = seed;
        (0..n)
            .map(|_| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                b"ACGT"[(x % 4) as usize]
            })
            .collect()
    }

    #[cfg(feature = "karan")]
    fn rc(s: &[u8]) -> Vec<u8> {
        s.iter().rev().map(|&c| comp(c)).collect()
    }

    #[test]
    fn spelled_with_edits() {
        // Walk >1 <2: ACGTA then TTGC reversed-complemented (GCAA).
        let seq = |n: usize| match n {
            1 => Some(b"ACGTA".to_vec()),
            2 => Some(b"TTGC".to_vec()),
            _ => None,
        };
        let walk = [(1, false), (2, true)];
        assert_eq!(spell(&walk, &seq, 0, ":9").unwrap(), b"ACGTAGCAA");
        // From base 1: a substitution, an insertion, a deletion.
        assert_eq!(spell(&walk, &seq, 1, ":2*ta+gg:1-gc:2").unwrap(), b"CGAGGAAA");
        assert_eq!(spell(&walk, &seq, 0, ":10"), None);
        assert_eq!(spell(&[(3, false)], &seq, 0, ":1"), None);
    }

    #[test]
    fn tags_round_trip() {
        let r = Rescue { pos: 1234, reverse: true, cigar: vec![(S, 5), (M, 100), (I, 2), (M, 50), (D, 3), (M, 7)], mapq: 60 };
        assert_eq!(r.span(), 160);
        let t = r.tag("chr16");
        assert_eq!(t, "chr16:1234:-:5S100M2I50M3D7M");
        assert_eq!(parse_tag(&t), Some(("chr16", 1234, true, r.cigar.clone())));
        assert_eq!(parse_tag("chr1:5:+:10Q"), None);
        assert_eq!(parse_tag("CHM13#0#chr1:5:+:10M").map(|t| t.0), Some("CHM13#0#chr1"));
    }

    #[test]
    #[cfg(feature = "karan")]
    fn found_on_either_strand_with_edits() {
        let reference = random(200_000, 7);
        // A read from 120,000 with a substitution and a 2-letter deletion: found forward and reverse.
        let mut read = reference[120_000..135_000].to_vec();
        read[5_000] = if read[5_000] == b'A' { b'C' } else { b'A' };
        read.drain(9_000..9_002);
        let letters = vec![reference.clone()];
        let index = Index::build(&letters, &[(0, 0, 200_000)]);
        let f = rescue(&read, &index, 0, &reference, 0, 200_000).unwrap();
        assert_eq!((f.pos, f.reverse, f.mapq, f.span()), (120_000, false, 60, 15_000));
        let count = |op: u8| f.cigar.iter().filter(|o| o.0 == op).map(|o| o.1).sum::<u32>();
        assert_eq!((count(M), count(D), count(I), count(S)), (14_998, 2, 0, 0));
        let r = rescue(&rc(&read), &index, 0, &reference, 0, 200_000).unwrap();
        assert_eq!((r.pos, r.reverse, r.span()), (120_000, true, 15_000));
        // Outside the stretch searched, or too divergent: none.
        assert_eq!(rescue(&read, &index, 0, &reference, 0, 100_000), None);
        assert_eq!(rescue(&read, &index, 1, &reference, 0, 200_000), None);
        let noisy: Vec<u8> = read.iter().enumerate().map(|(i, &c)| if i % 50 == 0 { comp(c) } else { c }).collect();
        assert_eq!(rescue(&noisy, &index, 0, &reference, 0, 200_000), None);
    }

    #[test]
    #[cfg(feature = "karan")]
    fn a_second_copy_lowers_the_mapping_quality() {
        let mut reference = random(100_000, 3);
        let copy = reference[10_000..25_000].to_vec();
        reference.extend_from_slice(&random(20_000, 5));
        reference.extend_from_slice(&copy);
        let index = Index::build(&[reference.clone()], &[(0, 0, 60_000), (0, 50_000, reference.len() as u64)]);
        let r = rescue(&copy, &index, 0, &reference, 0, reference.len() as u64).unwrap();
        assert_eq!(r.mapq, 0);
    }

    #[test]
    fn full_depth_takes_no_read() {
        // Contig 0: 20 alignments over 0-30,000 and 2 over 30,000-40,000 (median 20); contig 1 uncovered.
        let mut d = Depth::default();
        for _ in 0..20 {
            d.add(0, 0, 30_000);
        }
        for _ in 0..2 {
            d.add(0, 30_000, 40_000);
        }
        d.finish();
        assert!(d.full(0, 2_000, 8_000));
        assert!(!d.full(0, 32_000, 38_000));
        // Half at 20, half at 2: 11.
        assert!(!d.full(0, 25_000, 35_000));
        assert!(!d.full(0, 50_000, 60_000));
        assert!(!d.full(1, 0, 1_000));
    }

    #[test]
    fn three_copies_depth_takes_no_read() {
        // Median 20 over 0-30,000; 2 over 30,000-40,000, filled back by 16 more (0.9) or by 24 (1.3).
        let mut d = Depth::default();
        for _ in 0..20 {
            d.add(0, 0, 30_000);
        }
        for _ in 0..2 {
            d.add(0, 30_000, 40_000);
        }
        d.finish();
        let (mut two, mut three) = (d.clone(), d.clone());
        for _ in 0..16 {
            two.add(0, 30_000, 40_000);
        }
        for _ in 0..24 {
            three.add(0, 30_000, 40_000);
        }
        assert!((two.share(0, 31_000, 39_000) - 0.9).abs() < 1e-9);
        assert!(two.share(0, 31_000, 39_000) < MAX_RESCUED_DEPTH);
        assert!(three.share(0, 31_000, 39_000) >= MAX_RESCUED_DEPTH);
        assert!(!three.full(0, 50_000, 60_000));
    }

    #[test]
    #[cfg(feature = "karan")]
    fn repeats_are_left_out() {
        // A 171-letter unit repeated 300 times: every minimizer occurs 300 times, none is indexed.
        let unit = random(171, 11);
        let array: Vec<u8> = (0..300).flat_map(|_| unit.iter().copied()).collect();
        let index = Index::build(&[array.clone()], &[(0, 0, array.len() as u64)]);
        assert_eq!(index.len(), 0);
        assert_eq!(rescue(&array[1_000..16_000], &index, 0, &array, 0, array.len() as u64), None);
    }
}
