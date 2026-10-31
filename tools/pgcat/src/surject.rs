//! `pgcat surject` (experimental): the reads' graph alignments (GAF) as alignments to one reference contig (SAM),
//! from the graph alone -- pgphase's `--bam` without vg surject.
//!
//! Each read is rebuilt from its walk's node sequences and its `cs` string (qualities from `bq`). A read base aligned
//! to a node the reference path visits once takes that node's position; the read's orientation is that of most such
//! bases. Placed stretches (blocks) that continue in both read and reference are joined; between two blocks, a gap
//! as long in the read as on the reference (a substitution bubble) is aligned base for base, one with nothing on one
//! side is an insertion or a deletion, and any other (an indel bubble) is aligned globally, read against reference
//! (affine gaps, in a band around the diagonal as wide as the two stretches' length difference and BAND_SLACK more;
//! past MAX_CELLS band cells, an insertion then a deletion). The blocks kept are the best chain (`chain`):
//! the most placed bases less a concave cost of the difference between the read's and the reference's step from one
//! block to the next, as minimap2 chains seeds -- a walk touching reference nodes far apart (satellites) is not
//! stitched across the jump; the blocks left out fall into gaps. The read's bases before the first block and after
//! the last are aligned to the reference beyond them, the far end on the reference free (as vg surject aligns read
//! ends off the reference), and kept if at least END_MATCHES of them match; else they are soft-clipped. Reads with
//! no base placed on the contig are left out.
//!
//! Each end is then clipped where the read diverges from its own graph alignment (BWA-MEM's local scoring: match 1,
//! mismatch -4, gap -6 less one a base, clip penalty 5; the cut before the matching base where the prefix score is
//! lowest, if that gains more than the penalty). Scored on the read's cs against its walk, not against the reference:
//! another haplotype's nodes realigned to CHM13 look like divergence. A gap longer than SMALL_BUBBLE scores 0 (an
//! allele of the pangenome's).
//!
//! With `graph_placement` (--graph-placement) a read base is placed only where the graph puts it: a stretch of more
//! than SMALL_BUBBLE read or reference bases between two blocks (sequence on another haplotype's or copy's nodes) is
//! not aligned to the reference -- the read's letters are an insertion and the reference's a deletion -- and a read
//! end longer than SMALL_BUBBLE is clipped. Substitution and small indel bubbles (the pangenome's own alleles) are
//! placed as before.
//!
//! With `catalog` (--catalog, the contig's .pgsites) a read spanning one of the catalog's snarl windows (pgcat
//! pileup's: a snarl's reference bases between its boundary nodes, padded, those in a repeat merged) is aligned there
//! as its allele: its own letters between the window's two anchors against the reference's, at ASM_COSTS (the truth
//! sets' spelling) -- whatever walk and edits the graph alignment used. Its bases are placed by the graph elsewhere.

use crate::reference::ReferencePath;
use gbz::GBZ;
use rayon::prelude::*;
use std::collections::HashMap;
use std::io::{BufRead, BufReader, Read, Write};

pub struct SurjectOptions {
    pub graph: String,
    pub gaf: String,
    pub contig: String,
    pub reference_sample: String,
    /// The contig's name in the SAM (the graph's, CHM13#0#chr20, unless given).
    pub sam_contig: Option<String>,
    /// Place read bases only where the graph puts them (no alignment of long stretches off the reference).
    pub graph_placement: bool,
    /// The contig's catalog (.pgsites): reads aligned at its snarl windows as their alleles.
    pub catalog: Option<String>,
}

/// Band cells of an alignment (an indel bubble's read and reference stretches, or a read end), at most; and the band
/// beyond the two stretches' length difference.
const MAX_CELLS: usize = 60_000_000;
const BAND_SLACK: usize = 128;
/// The band of a read end's alignment beyond BAND_SLACK, at most: its indels against the reference beyond its block.
const END_BAND: usize = 2_000;
/// Alignment costs (a match costs nothing).
#[derive(Clone, Copy)]
pub struct Costs {
    pub mismatch: i32,
    pub open: i32,
    pub extend: i32,
}

/// A read against the reference, gaps cheap enough to follow HiFi indels.
pub const READ_COSTS: Costs = Costs { mismatch: 5, open: 6, extend: 3 };
/// An allele against the reference as the truth sets spell it: minimap2's asm5 scores (match 1, mismatch 19, gap
/// 39 + 3 per base) as costs at equal lengths (twice: mismatch 40, open 78, extend 7) -- mismatches before gaps.
pub const ASM_COSTS: Costs = Costs { mismatch: 40, open: 78, extend: 7 };
/// The longest stretch off the reference (read or reference side) graph placement still aligns: a small bubble.
const SMALL_BUBBLE: usize = 50;
/// The share of a read end's bases an end alignment must match for the end to be aligned rather than clipped: a
/// read end of the same locus matches nearly all of it; another copy of a satellite array or duplication, at 70-85%,
/// is clipped as vg surject clips it (aligned, at 70%, it gave reads of chr20's 32 Mb array twice vg's mismatches,
/// past pgphase's 5% per read, and the reads were left out of the calls beside it).
const END_MATCHES: f64 = 0.9;

pub const M: u8 = b'M';
pub const I: u8 = b'I';
pub const D: u8 = b'D';
pub const S: u8 = b'S';

fn push(cigar: &mut Vec<(u8, u32)>, op: u8, n: usize) {
    if n == 0 {
        return;
    }
    match cigar.last_mut() {
        Some(last) if last.0 == op => last.1 += n as u32,
        _ => cigar.push((op, n as u32)),
    }
}

/// The affine-gap alignment of `q` (read) against `r` (reference) from both starts, in a band around the diagonal
/// (|j - i| <= half), as one traceback byte per cell (`(n + 1) * (2 half + 1)`, cell (i, j) at `i * w + j - i +
/// half`): bits 0-1 where the best score came from (0 a match or mismatch, 1 a deletion, 2 an insertion), bit 2 a
/// deletion opened there, bit 3 an insertion opened there. Scores are kept two rows at a time; the last row's best
/// scores are kept for a free reference end.
struct Banded {
    half: usize,
    w: usize,
    tb: Vec<u8>,
    last: Vec<i32>,
}

const INF: i32 = i32::MAX / 4;

impl Banded {
    fn at(&self, i: usize, j: usize) -> Option<usize> {
        let d = j as i64 - i as i64 + self.half as i64;
        (d >= 0 && d < self.w as i64).then(|| i * self.w + d as usize)
    }

    fn fill(q: &[u8], r: &[u8], half: usize, c: Costs) -> Option<Banded> {
        let (n, m) = (q.len(), r.len());
        let w = 2 * half + 1;
        if (n + 1) * w > MAX_CELLS {
            return None;
        }
        let mut tb = vec![0u8; (n + 1) * w];
        // Rows by band offset d = j - i + half: (i, j - 1) is d - 1 in the same row, (i - 1, j) is d + 1 in the row
        // before, (i - 1, j - 1) is d in the row before.
        let (mut h_prev, mut f_prev) = (vec![INF; w], vec![INF; w]);
        let (mut h_cur, mut e_cur, mut f_cur) = (vec![INF; w], vec![INF; w], vec![INF; w]);
        for i in 0..=n {
            h_cur.iter_mut().for_each(|x| *x = INF);
            e_cur.iter_mut().for_each(|x| *x = INF);
            f_cur.iter_mut().for_each(|x| *x = INF);
            let lo = i.saturating_sub(half);
            let hi = (i + half).min(m);
            for j in lo..=hi {
                let d = j + half - i;
                if i == 0 && j == 0 {
                    h_cur[d] = 0;
                    continue;
                }
                let mut bits = 0u8;
                let mut e = INF;
                if j > 0 && d > 0 {
                    let (open, ext) = (h_cur[d - 1].saturating_add(c.open + c.extend), e_cur[d - 1].saturating_add(c.extend));
                    e = open.min(ext);
                    if open <= ext {
                        bits |= 4;
                    }
                }
                let mut f = INF;
                if i > 0 && d + 1 < w {
                    let (open, ext) = (h_prev[d + 1].saturating_add(c.open + c.extend), f_prev[d + 1].saturating_add(c.extend));
                    f = open.min(ext);
                    if open <= ext {
                        bits |= 8;
                    }
                }
                let diag = if i > 0 && j > 0 {
                    h_prev[d].saturating_add(if q[i - 1].eq_ignore_ascii_case(&r[j - 1]) { 0 } else { c.mismatch })
                } else {
                    INF
                };
                let h = if diag <= e && diag <= f {
                    diag
                } else if e <= f {
                    bits |= 1;
                    e
                } else {
                    bits |= 2;
                    f
                };
                h_cur[d] = h;
                e_cur[d] = e;
                f_cur[d] = f;
                tb[i * w + d] = bits;
            }
            std::mem::swap(&mut h_prev, &mut h_cur);
            std::mem::swap(&mut f_prev, &mut f_cur);
        }
        Some(Banded { half, w, tb, last: h_prev })
    }

    /// The best score at (n, j) of the last row, if j is in the band.
    fn last_score(&self, n: usize, j: usize) -> Option<i32> {
        let d = j as i64 - n as i64 + self.half as i64;
        (d >= 0 && d < self.w as i64).then(|| self.last[d as usize])
    }

    /// The operations from (0, 0) to (i, j); on ties a match or mismatch, so gaps go as far left as the costs allow.
    fn trace(&self, mut i: usize, mut j: usize) -> Vec<(u8, u32)> {
        let mut state = 0u8; // 0 best, 1 in a deletion, 2 in an insertion
        let mut rev: Vec<u8> = Vec::with_capacity(i + j);
        while i > 0 || j > 0 {
            let bits = self.at(i, j).map_or(0, |k| self.tb[k]);
            match state {
                0 => match bits & 3 {
                    0 => {
                        rev.push(M);
                        i -= 1;
                        j -= 1;
                    }
                    1 => state = 1,
                    _ => state = 2,
                },
                1 => {
                    rev.push(D);
                    j -= 1;
                    if bits & 4 != 0 {
                        state = 0;
                    }
                }
                _ => {
                    rev.push(I);
                    i -= 1;
                    if bits & 8 != 0 {
                        state = 0;
                    }
                }
            }
        }
        let mut out = Vec::new();
        for &op in rev.iter().rev() {
            push(&mut out, op, 1);
        }
        out
    }
}

/// The band half-width an alignment of an `n`-letter and an `m`-letter stretch needs: their difference, and
/// BAND_SLACK more for indels that cancel out.
fn half_band(n: usize, m: usize) -> usize {
    n.abs_diff(m) + BAND_SLACK
}

/// The global alignment of `q` (read) against `r` (reference) with affine gaps, as CIGAR operations; past MAX_CELLS
/// band cells, an insertion then a deletion.
pub fn align_global(q: &[u8], r: &[u8]) -> Vec<(u8, u32)> {
    align_global_with(q, r, READ_COSTS)
}

/// align_global at the given costs.
pub fn align_global_with(q: &[u8], r: &[u8], c: Costs) -> Vec<(u8, u32)> {
    let (n, m) = (q.len(), r.len());
    let fallback = || {
        let mut out = Vec::new();
        push(&mut out, I, n);
        push(&mut out, D, m);
        out
    };
    if n == 0 || m == 0 {
        return fallback();
    }
    match Banded::fill(q, r, half_band(n, m), c) {
        Some(b) => b.trace(n, m),
        None => fallback(),
    }
}

/// `q` aligned from its start against `r` from its start, `r`'s far end free: the CIGAR and how many letters of `r`
/// it uses, or None if it matches fewer than END_MATCHES of `q`'s letters (or is too large to align).
pub fn align_end(q: &[u8], r: &[u8]) -> Option<(Vec<(u8, u32)>, usize)> {
    let (n, m) = (q.len(), r.len());
    if n == 0 {
        return Some((Vec::new(), 0));
    }
    if m == 0 {
        return None;
    }
    // The far end is free, so the band need only hold the indels a read end may have, not the reference's slack.
    let b = Banded::fill(q, r, n.abs_diff(m).min(END_BAND) + BAND_SLACK, READ_COSTS)?;
    let end = (1..=m).filter_map(|j| b.last_score(n, j).map(|s| (s, j))).min()?.1;
    let ops = b.trace(n, end);
    let mut matches = 0usize;
    let (mut i, mut j) = (0usize, 0usize);
    for &(op, len) in &ops {
        for _ in 0..len {
            match op {
                M => {
                    matches += q[i].eq_ignore_ascii_case(&r[j]) as usize;
                    i += 1;
                    j += 1;
                }
                I => i += 1,
                _ => j += 1,
            }
        }
    }
    (matches as f64 >= END_MATCHES * n as f64).then_some((ops, end))
}

/// Predecessors a block's chaining looks back over, and the farthest step in read or reference it may join across.
const CHAIN_LOOKBACK: usize = 64;
const CHAIN_MAX_STEP: u64 = 100_000;

/// The cost of joining two blocks whose read and reference steps differ by `d` bases: minimap2's 0.01 w |d| +
/// 0.5 log2 |d| at w 5 (an indel of d bases joins blocks of more than about d / 20 bases).
fn join_cost(d: u64) -> f64 {
    if d == 0 { 0.0 } else { 0.05 * d as f64 + 0.5 * (d as f64).log2() }
}

/// The best chain of blocks (read offset, reference position, length; sorted by read offset): each later than the
/// one before in both read and reference, with the most bases less the join costs. The chain's blocks, in order.
pub fn chain(blocks: &[(usize, u64, usize)]) -> Vec<(usize, u64, usize)> {
    let n = blocks.len();
    if n <= 1 {
        return blocks.to_vec();
    }
    let mut score = vec![0f64; n];
    let mut prev = vec![usize::MAX; n];
    for j in 0..n {
        let (qj, rj, lj) = blocks[j];
        score[j] = lj as f64;
        for i in j.saturating_sub(CHAIN_LOOKBACK)..j {
            let (qi, ri, li) = blocks[i];
            let (qe, re) = (qi + li, ri + li as u64);
            if qe > qj || re > rj {
                continue;
            }
            let (dq, dr) = ((qj - qe) as u64, rj - re);
            if dq > CHAIN_MAX_STEP || dr > CHAIN_MAX_STEP {
                continue;
            }
            let s = score[i] + lj as f64 - join_cost(dq.abs_diff(dr));
            if s > score[j] {
                score[j] = s;
                prev[j] = i;
            }
        }
    }
    let mut j = (0..n).max_by(|&a, &b| score[a].partial_cmp(&score[b]).unwrap().then(b.cmp(&a))).unwrap();
    let mut out = vec![blocks[j]];
    while prev[j] != usize::MAX {
        j = prev[j];
        out.push(blocks[j]);
    }
    out.reverse();
    out
}

/// A read's placed blocks (read offset, 0-based reference position, length; the BAM's orientation, sorted by read
/// offset) as an alignment: its first reference position and CIGAR. `read` is the read in that orientation;
/// `reference(a, b)` the reference's letters a..b. None if no block. `graph`: graph placement (see the module notes).
// End clipping, as BWA-MEM scores a local alignment, on the read's alignment to the graph (its cs edits), not on the
// alignment to the reference rebuilt from it: a base matching the walk +1, a mismatch -4, an insertion or deletion
// -(6 + its length) (free beyond SMALL_BUBBLE: a structural variant), and an end clipped only where that gains more than
// CLIP_PENALTY. A read's end in sequence of another copy -- a satellite array beside the region, its walk a fifth
// mismatched -- is clipped, as vg's surjection clips it; aligned, it left the read over the 5% difference pgphase skips
// reads at (chr20:32.2 Mb, the alt reads of 150 truth hets). Scored on the reference instead, a VNTR allele the read
// walks (other haplotypes' nodes, realigned to CHM13 as scattered mismatches) read as divergence: 43 of 75 reads were
// clipped inside chr20:49,574,839-49,582,911 and 21 truth hets lost their MSA's reads.
const CLIP_MATCH: i64 = 1;
const CLIP_MISMATCH: i64 = -4;
const CLIP_GAP_OPEN: i64 = -6;
const CLIP_PENALTY: i64 = 5;

/// The score of an indel of `n` bases at the clipping scores.
fn clip_gap(n: usize) -> i64 {
    if n > SMALL_BUBBLE { 0 } else { CLIP_GAP_OPEN - n as i64 }
}

/// Read bases to clip at each end, from per-base scores (an indel's score on the base before it): the cut before a
/// matching base with the lowest prefix score, and after one with the lowest suffix score, where that gains more than
/// CLIP_PENALTY.
pub fn clip_lengths(score: &[i64]) -> (usize, usize) {
    let n = score.len();
    let (mut sum, mut best, mut left) = (0i64, 0i64, 0usize);
    for (i, &v) in score.iter().enumerate() {
        if v == CLIP_MATCH && sum < best {
            best = sum;
            left = i;
        }
        sum += v;
    }
    if -best <= CLIP_PENALTY {
        left = 0;
    }
    let (mut sum, mut best, mut right) = (0i64, 0i64, n);
    for i in (left..n).rev() {
        if score[i] == CLIP_MATCH && sum < best {
            best = sum;
            right = i + 1;
        }
        sum += score[i];
    }
    if -best <= CLIP_PENALTY {
        right = n;
    }
    (left, n - right)
}

/// The alignment (start, CIGAR) with `head` read bases soft-clipped at its start and `tail` at its end (soft clips
/// already there included in the counts); a gap left at either end goes with the clip.
pub fn clip_cigar(start: u64, cigar: &[(u8, u32)], head: usize, tail: usize) -> (u64, Vec<(u8, u32)>) {
    if head == 0 && tail == 0 {
        return (start, cigar.to_vec());
    }
    let qlen: usize = cigar.iter().filter(|o| o.0 != D).map(|o| o.1 as usize).sum();
    if head + tail >= qlen {
        return (start, cigar.to_vec());
    }
    // Expand to one op per read base or deleted base, cut, then re-run-length.
    let mut ops: Vec<u8> = Vec::new();
    for &(op, n) in cigar {
        ops.extend(std::iter::repeat(op).take(n as usize));
    }
    let (mut q, mut pos) = (0usize, start);
    let mut out: Vec<u8> = Vec::new();
    for &op in &ops {
        let read_base = op != D;
        let in_head = q < head;
        let in_tail = q >= qlen - tail;
        if in_head || (in_tail && read_base) {
            if read_base {
                out.push(S);
            }
            if in_head && (op == M || op == D) {
                pos += 1;
            }
        } else if in_tail && op == D {
            // a deletion after the last kept base: dropped
        } else {
            out.push(if op == S { S } else { op });
        }
        if read_base {
            q += 1;
        }
    }
    // A clip ends at an aligned base: gaps beside a soft clip join it (insertions) or move the start (deletions).
    let mut runs: Vec<(u8, u32)> = Vec::new();
    for op in out {
        push(&mut runs, op, 1);
    }
    loop {
        let k = runs.iter().position(|o| o.0 != S).unwrap_or(runs.len());
        match runs.get(k) {
            Some(&(I, n)) => { runs.remove(k); if k > 0 { runs[k - 1].1 += n } else { runs.insert(0, (S, n)) } }
            Some(&(D, n)) => { runs.remove(k); pos += n as u64; }
            _ => break,
        }
        let mut merged: Vec<(u8, u32)> = Vec::new();
        for &(op, n) in &runs { push(&mut merged, op, n as usize); }
        runs = merged;
    }
    loop {
        let k = runs.iter().rposition(|o| o.0 != S);
        match k.map(|k| (k, runs[k])) {
            Some((k, (I, n))) => { runs.remove(k); if k < runs.len() { runs[k].1 += n } else { runs.push((S, n)) } }
            Some((k, (D, _))) => { runs.remove(k); }
            _ => break,
        }
        let mut merged: Vec<(u8, u32)> = Vec::new();
        for &(op, n) in &runs { push(&mut merged, op, n as usize); }
        runs = merged;
    }
    (pos, runs)
}

pub fn linearize(read: &[u8], blocks: &[(usize, u64, usize)], reference: &dyn Fn(u64, u64) -> Vec<u8>, ref_len: u64, graph: bool) -> Option<(u64, Vec<(u8, u32)>)> {
    linearize_windows(read, blocks, reference, ref_len, graph, &[])
}

/// linearize, the read aligned at each window (sorted, disjoint; 0-based half-open reference ranges) it holds both
/// anchors of (the reference bases just before and just after it) as its allele, at ASM_COSTS.
pub fn linearize_windows(read: &[u8], blocks: &[(usize, u64, usize)], reference: &dyn Fn(u64, u64) -> Vec<u8>, ref_len: u64, graph: bool, windows: &[(u64, u64)]) -> Option<(u64, Vec<(u8, u32)>)> {
    let mut kept = chain(blocks);
    // The windows whose two anchors the chain places: their bases taken out of the blocks, so each is one gap.
    let placed = |r: u64| kept.iter().any(|&(_, r0, n)| r0 <= r && r < r0 + n as u64);
    let mut spanned: Vec<(u64, u64)> = Vec::new();
    if let (Some(f), Some(l)) = (kept.first(), kept.last()) {
        let (lo, hi) = (f.1, l.1 + l.2 as u64);
        let first = windows.partition_point(|w| w.1 < lo);
        for &(a, b) in windows[first..].iter().take_while(|w| w.0 <= hi) {
            if a > 0 && placed(a - 1) && placed(b) {
                spanned.push((a, b));
            }
        }
    }
    if !spanned.is_empty() {
        let mut cut: Vec<(usize, u64, usize)> = Vec::with_capacity(kept.len() + spanned.len());
        for &(q0, r0, n) in &kept {
            let mut parts = vec![(q0, r0, n)];
            for &(a, b) in &spanned {
                let mut next = Vec::new();
                for (q, r, len) in parts {
                    let e = r + len as u64;
                    if e <= a || r >= b {
                        next.push((q, r, len));
                        continue;
                    }
                    if r < a {
                        next.push((q, r, (a - r) as usize));
                    }
                    if e > b {
                        next.push((q + (b - r) as usize, b, (e - b) as usize));
                    }
                }
                parts = next;
            }
            cut.extend(parts);
        }
        kept = cut;
    }
    let is_window = |ra: u64, rb: u64| spanned.binary_search(&(ra, rb)).is_ok();
    let first = *kept.first()?;
    let mut cigar = Vec::new();
    let mut start = first.1;
    // The head: the read's first bases, reversed, against the reference before the first block, reversed.
    let head = first.0;
    let reach = (head + (head / 5).max(50)) as u64;
    let r0 = first.1.saturating_sub(reach);
    let aligned_head = (head > 0 && first.1 > 0 && !(graph && head > SMALL_BUBBLE))
        .then(|| {
            let q: Vec<u8> = read[..head].iter().rev().copied().collect();
            let r: Vec<u8> = reference(r0, first.1).into_iter().rev().collect();
            align_end(&q, &r)
        })
        .flatten();
    match aligned_head {
        Some((ops, used)) => {
            start = first.1 - used as u64;
            for &(op, n) in ops.iter().rev() {
                push(&mut cigar, op, n as usize);
            }
        }
        None => push(&mut cigar, S, head),
    }
    push(&mut cigar, M, first.2);
    for w in kept.windows(2) {
        let (a, b) = (w[0], w[1]);
        let (qa, ra) = (a.0 + a.2, a.1 + a.2 as u64);
        let (gq, gr) = (b.0 - qa, (b.1 - ra) as usize);
        if is_window(ra, b.1) {
            for (op, n) in align_global_with(&read[qa..b.0], &reference(ra, b.1), ASM_COSTS) {
                push(&mut cigar, op, n as usize);
            }
        } else if graph && gq.max(gr) > SMALL_BUBBLE {
            push(&mut cigar, I, gq);
            push(&mut cigar, D, gr);
        } else if gq == gr {
            push(&mut cigar, M, gq);
        } else if gr == 0 {
            push(&mut cigar, I, gq);
        } else if gq == 0 {
            push(&mut cigar, D, gr);
        } else {
            for (op, n) in align_global(&read[qa..b.0], &reference(ra, b.1)) {
                push(&mut cigar, op, n as usize);
            }
        }
        push(&mut cigar, M, b.2);
    }
    let last = *kept.last().unwrap();
    let (tq, tr) = (last.0 + last.2, last.1 + last.2 as u64);
    let tail = read.len() - tq;
    let reach = ((tail + (tail / 5).max(50)) as u64).min(ref_len.saturating_sub(tr));
    let aligned_tail = (tail > 0 && reach > 0 && !(graph && tail > SMALL_BUBBLE)).then(|| align_end(&read[tq..], &reference(tr, tr + reach))).flatten();
    match aligned_tail {
        Some((ops, _)) => {
            for (op, n) in ops {
                push(&mut cigar, op, n as usize);
            }
        }
        None => push(&mut cigar, S, tail),
    }
    // An end alignment may start or finish with a deletion: not a CIGAR's first or last operation.
    while let Some(&(D, n)) = cigar.first() {
        start += n as u64;
        cigar.remove(0);
    }
    while let Some(&(D, _)) = cigar.last() {
        cigar.pop();
    }
    Some((start, cigar))
}

pub(crate) fn comp(b: u8) -> u8 {
    match b {
        b'A' | b'a' => b'T',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        b'T' | b't' => b'A',
        _ => b'N',
    }
}

pub(crate) fn parse_walk(s: &str) -> Option<Vec<(usize, bool)>> {
    let b = s.as_bytes();
    let mut out = Vec::new();
    let mut i = 0;
    while i < b.len() {
        let rev = match b[i] {
            b'>' => false,
            b'<' => true,
            _ => return None,
        };
        i += 1;
        let st = i;
        while i < b.len() && b[i].is_ascii_digit() {
            i += 1;
        }
        out.push((s[st..i].parse().ok()?, rev));
    }
    Some(out)
}

/// The contig as the projection needs it: its letters, and each node the path visits once, with where and whether
/// reversed.
pub(crate) struct Contig {
    pub(crate) seq: Vec<u8>,
    pub(crate) place: HashMap<usize, (u64, bool)>,
}

pub(crate) fn load_contig(g: &GBZ, sample: &str, contig: &str) -> Result<(Contig, String), String> {
    let r = ReferencePath::load(g, sample, contig)?;
    let mut seq = vec![b'N'; r.length];
    let mut count: HashMap<usize, u32> = HashMap::new();
    for (k, h) in r.steps.iter().enumerate() {
        *count.entry(h.0).or_default() += 1;
        let s = g.sequence(h.0).ok_or("a reference node without sequence")?;
        let at = r.pos[k];
        if h.1 == gbz::support::Orientation::Forward {
            seq[at..at + s.len()].copy_from_slice(s);
        } else {
            for (i, &c) in s.iter().rev().enumerate() {
                seq[at + i] = comp(c);
            }
        }
    }
    let mut place = HashMap::new();
    for (k, h) in r.steps.iter().enumerate() {
        if count[&h.0] == 1 {
            place.insert(h.0, (r.pos[k] as u64, h.1 == gbz::support::Orientation::Reverse));
        }
    }
    Ok((Contig { seq, place }, r.chrom))
}

/// One GAF line as a SAM record (without the newline), or None if it places no base on the contig.
fn surject_line(g: &GBZ, c: &Contig, sam_contig: &str, graph: bool, windows: &[(u64, u64)], line: &str) -> Option<String> {
    let f: Vec<&str> = line.split('\t').collect();
    let k = if f.len() > 4 && (f[4] == "+" || f[4] == "-") { 0 } else { 3 };
    if f.len() < k + 12 || f[k + 5] == "*" {
        return None;
    }
    let name = f[k];
    let minus = f[k + 4] == "-";
    let (qlen, qstart, qend): (usize, usize, usize) = (f[k + 1].parse().ok()?, f[k + 2].parse().ok()?, f[k + 3].parse().ok()?);
    let pstart: usize = f[k + 7].parse().ok()?;
    let mapq: u32 = f[k + 11].parse().unwrap_or(255);
    let cs = f.iter().find_map(|t| t.strip_prefix("cs:Z:"))?;
    let bq = f.iter().find_map(|t| t.strip_prefix("bq:Z:"));
    let walk = parse_walk(f[k + 5])?;
    let lens: Vec<usize> = walk.iter().map(|&(n, _)| g.sequence_len(n).unwrap_or(0)).collect();
    let (mut w, mut off) = (0usize, pstart);
    while w < walk.len() && off >= lens[w] {
        off -= lens[w];
        w += 1;
    }
    // Read letters in walk orientation, and placed runs (walk offset, length, reference position of the first
    // letter, +1 or -1).
    let mut read: Vec<u8> = Vec::new();
    let mut runs: Vec<(usize, usize, u64, i64)> = Vec::new();
    // Each read letter's clipping score against the walk (clip_lengths); an indel's on the letter before it.
    let mut score: Vec<i64> = Vec::new();
    let place = |w: usize, off: usize| -> Option<(u64, i64)> {
        let (node, rev) = walk[w];
        let &(start, node_rev) = c.place.get(&node)?;
        let len = lens[w];
        let fo = if rev { len - 1 - off } else { off };
        let r = if node_rev { start + (len - 1 - fo) as u64 } else { start + fo as u64 };
        Some((r, if rev == node_rev { 1 } else { -1 }))
    };
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
                let mut left: usize = cs[st..i].parse().ok()?;
                while left > 0 {
                    if w >= walk.len() {
                        return None;
                    }
                    let (node, rev) = walk[w];
                    let s = g.sequence(node)?;
                    let take = left.min(lens[w] - off);
                    if rev {
                        read.extend(s[lens[w] - off - take..lens[w] - off].iter().rev().map(|&x| comp(x)));
                    } else {
                        read.extend(s[off..off + take].iter().map(|x| x.to_ascii_uppercase()));
                    }
                    score.extend(std::iter::repeat(CLIP_MATCH).take(take));
                    if let Some((r, d)) = place(w, off) {
                        runs.push((read.len() - take, take, r, d));
                    }
                    off += take;
                    left -= take;
                    if off == lens[w] {
                        w += 1;
                        off = 0;
                    }
                }
            }
            b'*' => {
                if i + 2 > b.len() || w >= walk.len() {
                    return None;
                }
                read.push(b[i + 1].to_ascii_uppercase());
                i += 2;
                score.push(CLIP_MISMATCH);
                if let Some((r, d)) = place(w, off) {
                    runs.push((read.len() - 1, 1, r, d));
                }
                off += 1;
                if off == lens[w] {
                    w += 1;
                    off = 0;
                }
            }
            b'+' => {
                while i < b.len() && b[i].is_ascii_alphabetic() {
                    i += 1;
                }
                read.extend(b[st..i].iter().map(|x| x.to_ascii_uppercase()));
                let n = i - st;
                score.extend(std::iter::repeat(0).take(n));
                if let Some(v) = score.len().checked_sub(n + 1).and_then(|k| score.get_mut(k)) {
                    *v += clip_gap(n);
                }
            }
            b'-' => {
                while i < b.len() && b[i].is_ascii_alphabetic() {
                    i += 1;
                }
                let mut left = i - st;
                if let Some(v) = score.last_mut() {
                    *v += clip_gap(left);
                }
                while left > 0 {
                    if w >= walk.len() {
                        return None;
                    }
                    let take = left.min(lens[w] - off);
                    off += take;
                    left -= take;
                    if off == lens[w] {
                        w += 1;
                        off = 0;
                    }
                }
            }
            _ => return None,
        }
    }
    let fwd_bases: usize = runs.iter().filter(|r| r.3 == 1).map(|r| r.1).sum();
    let rev_bases: usize = runs.iter().filter(|r| r.3 == -1).map(|r| r.1).sum();
    if fwd_bases + rev_bases == 0 {
        return None;
    }
    let fwd = fwd_bases >= rev_bases;
    let rl = read.len();
    let mut blocks: Vec<(usize, u64, usize)> = runs
        .iter()
        .filter(|r| (r.3 == 1) == fwd)
        .map(|&(q, n, r, _)| if fwd { (q, r, n) } else { (rl - (q + n), r + 1 - n as u64, n) })
        .collect();
    blocks.sort_unstable();
    let mut merged: Vec<(usize, u64, usize)> = Vec::with_capacity(blocks.len());
    for blk in blocks {
        if let Some(last) = merged.last_mut() {
            if last.0 + last.2 == blk.0 && last.1 + last.2 as u64 == blk.1 {
                last.2 += blk.2;
                continue;
            }
        }
        merged.push(blk);
    }
    // The read in the BAM's orientation: a GAF read on `-` is the path's letters reverse-complemented.
    let seq: Vec<u8> = if fwd { read } else { read.iter().rev().map(|&x| comp(x)).collect() };
    // bq covers the whole read; the aligned part (qstart..qend, on the read as sequenced) is what was rebuilt. The
    // GAF holds no letters for unaligned read ends, so they are not in the record.
    let mut qual: Vec<u8> = match bq {
        Some(q) if q.len() == qlen && qend <= qlen && qend - qstart == rl => q.as_bytes()[qstart..qend].to_vec(),
        _ => vec![b'I'; rl],
    };
    if minus != !fwd {
        // bq runs along the read as sequenced; the walk's letters along the read reversed when it is on `-`.
        qual.reverse();
    }
    let reference = |a: u64, b: u64| c.seq[a as usize..b as usize].to_vec();
    let (pos, cigar) = linearize_windows(&seq, &merged, &reference, c.seq.len() as u64, graph, windows)?;
    // Clipping on the read's alignment to the walk, in the BAM's orientation.
    if !fwd {
        score.reverse();
    }
    let (head, tail) = clip_lengths(&score);
    let lead = match cigar.first() { Some(&(S, n)) => n as usize, _ => 0 };
    let trail = match cigar.last() { Some(&(S, n)) if cigar.len() > 1 => n as usize, _ => 0 };
    let (pos, cigar) = clip_cigar(pos, &cigar, head.max(lead), tail.max(trail));
    let cigar_s: String = cigar.iter().map(|&(op, n)| format!("{n}{}", op as char)).collect();
    Some(format!(
        "{name}\t{}\t{sam_contig}\t{}\t{mapq}\t{cigar_s}\t*\t0\t0\t{}\t{}",
        if fwd { 0 } else { 16 },
        pos + 1,
        String::from_utf8_lossy(&seq),
        String::from_utf8_lossy(&qual)
    ))
}

pub fn run(o: &SurjectOptions) -> Result<(), String> {
    let g: GBZ = simple_sds::serialize::load_from(&o.graph).map_err(|e| format!("{}: {e}", o.graph))?;
    let (contig, chrom) = load_contig(&g, &o.reference_sample, &o.contig)?;
    let sam_contig = o.sam_contig.clone().unwrap_or(chrom);
    let windows = match &o.catalog {
        Some(p) => crate::pileup::load_windows(&g, &contig, p, None)?,
        None => Vec::new(),
    };
    if o.catalog.is_some() {
        eprintln!("pgcat surject: {} snarl windows", windows.len());
    }
    let file = std::fs::File::open(&o.gaf).map_err(|e| format!("{}: {e}", o.gaf))?;
    let input: Box<dyn Read> = if o.gaf.ends_with(".gz") { Box::new(flate2::read::MultiGzDecoder::new(file)) } else { Box::new(file) };
    let mut lines = BufReader::with_capacity(1 << 22, input).lines();
    let stdout = std::io::stdout();
    let mut w = std::io::BufWriter::with_capacity(1 << 24, stdout.lock());
    let io = |e: std::io::Error| e.to_string();
    writeln!(w, "@HD\tVN:1.6\tSO:unsorted").map_err(io)?;
    writeln!(w, "@SQ\tSN:{sam_contig}\tLN:{}", contig.seq.len()).map_err(io)?;
    writeln!(w, "@PG\tID:pgcat-surject\tPN:pgcat\tVN:{}", env!("CARGO_PKG_VERSION")).map_err(io)?;
    let (mut n, mut written) = (0u64, 0u64);
    loop {
        let batch: Vec<String> = lines.by_ref().take(20_000).filter_map(|l| l.ok()).filter(|l| !l.starts_with('#') && !l.is_empty()).collect();
        if batch.is_empty() {
            break;
        }
        let out: Vec<Option<String>> = batch.par_iter().map(|l| surject_line(&g, &contig, &sam_contig, o.graph_placement, &windows, l)).collect();
        for rec in out {
            n += 1;
            if let Some(r) = rec {
                writeln!(w, "{r}").map_err(io)?;
                written += 1;
            }
        }
    }
    w.flush().map_err(io)?;
    eprintln!("pgcat surject: {n} alignments, {written} written, {} without a base on {}", n - written, o.contig);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cigar(c: &[(u8, u32)]) -> String {
        c.iter().map(|&(op, n)| format!("{n}{}", op as char)).collect()
    }

    #[test]
    fn global_alignment() {
        assert_eq!(cigar(&align_global(b"ACGT", b"ACGT")), "4M");
        assert_eq!(cigar(&align_global(b"ACTT", b"ACGT")), "4M");
        // A deletion in a run goes to its left end.
        assert_eq!(cigar(&align_global(b"GAAC", b"GAAAC")), "1M1D3M");
        assert_eq!(cigar(&align_global(b"GAAAAC", b"GAAC")), "1M2I3M");
        assert_eq!(cigar(&align_global(b"", b"ACG")), "3D");
        assert_eq!(cigar(&align_global(b"ACG", b"")), "3I");
    }

    #[test]
    fn long_gap_is_aligned_in_a_band() {
        // A 3 kb stretch with two mismatches and a 2-letter deletion (once past the old 4M-cell cap: written as a
        // 3 kb insertion and a 3 kb deletion).
        let mut x: u64 = 7;
        let r: Vec<u8> = (0..3000)
            .map(|_| {
                x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                b"ACGT"[(x >> 33) as usize % 4]
            })
            .collect();
        let mut q = r.clone();
        q[500] = if q[500] == b'A' { b'C' } else { b'A' };
        q[2500] = if q[2500] == b'A' { b'C' } else { b'A' };
        q.drain(1500..1502);
        let got = align_global(&q, &r);
        let d: u32 = got.iter().filter(|o| o.0 == D).map(|o| o.1).sum();
        let i: u32 = got.iter().filter(|o| o.0 == I).map(|o| o.1).sum();
        let m: u32 = got.iter().filter(|o| o.0 == M).map(|o| o.1).sum();
        assert_eq!((m, i, d), (2998, 0, 2));
        // A read end of 3 kb against the reference beyond it: aligned, not clipped.
        let mut rr = r.clone();
        rr.extend_from_slice(b"TTTTTTTTTT");
        let (ops, used) = align_end(&q, &rr).unwrap();
        assert_eq!(used, 3000);
        assert_eq!(ops.iter().filter(|o| o.0 == M).map(|o| o.1).sum::<u32>(), 2998);
    }

    #[test]
    fn chaining() {
        // A block far off the others' diagonal (a walk touching a reference node 50 kb away) is left out.
        let got = chain(&[(0, 1000, 500), (500, 51_500, 20), (520, 1520, 500)]);
        assert_eq!(got, vec![(0, 1000, 500), (520, 1520, 500)]);
        // A real 3 kb deletion between long blocks is kept.
        let got = chain(&[(0, 1000, 5000), (5000, 9000, 5000)]);
        assert_eq!(got.len(), 2);
    }

    #[test]
    fn linearized_blocks() {
        let reference = b"AAAACCCCGGGGTTTTAAAACCCC".to_vec();
        let at = |a: u64, b: u64| reference[a as usize..b as usize].to_vec();
        // Two blocks with a substitution bubble between (equal gaps): one match run.
        let read = b"AAAACCCCGGGG";
        let len = reference.len() as u64;
        assert_eq!(linearize(read, &[(0, 0, 4), (6, 6, 6)], &at, len, false).map(|(p, c)| (p, cigar(&c))), Some((0, "12M".into())));
        // An insertion, then a deletion, between blocks; read ends that match nothing beyond the blocks clipped.
        let read = b"TTAAAAXXCCCCGGNN";
        let got = linearize(read, &[(2, 0, 4), (8, 4, 4), (12, 10, 2)], &at, len, false).map(|(p, c)| (p, cigar(&c)));
        assert_eq!(got, Some((0, "2S4M2I4M2D2M2S".into())));
        // A block running backwards on the reference is dropped into the gap.
        let got = linearize(b"AAAACCCCGGGG", &[(0, 0, 4), (4, 2, 2), (6, 6, 6)], &at, len, false).map(|(p, c)| (p, cigar(&c)));
        assert_eq!(got, Some((0, "12M".into())));
        // An indel bubble: aligned against the reference between the blocks.
        let read = b"AAAACCCGGGGTTTT";
        let got = linearize(read, &[(0, 0, 4), (11, 12, 4)], &at, len, false).map(|(p, c)| (p, cigar(&c)));
        assert_eq!(got, Some((0, "4M1D11M".into())));
        // Read ends that match the reference beyond the blocks are aligned, not clipped.
        let got = linearize(b"CCCGGGGTTTTAAA", &[(3, 8, 8)], &at, len, false).map(|(p, c)| (p, cigar(&c)));
        assert_eq!(got, Some((5, "14M".into())));
        assert_eq!(linearize(b"ACGT", &[], &at, len, false), None);
    }

    #[test]
    fn window_aligned_as_an_allele() {
        // Reference: AAAA GAGAGAGA CCCC; a read's graph alignment put a 2-letter deletion in one place and the same
        // 2 letters back as an insertion further on (blocks either side of a window 4..12); aligned as the window's
        // allele (the read's letters between the anchors, equal length, one letter different) it is matches and one
        // mismatch, not an insertion and a deletion.
        let r = b"AAAAGAGAGAGACCCC".to_vec();
        let at = |a: u64, b: u64| r[a as usize..b as usize].to_vec();
        let read = b"AAAAGAGAGTGACCCC";
        let blocks = [(0usize, 0u64, 6usize), (6, 8, 2), (10, 10, 6)];
        let got = linearize_windows(read, &blocks, &at, r.len() as u64, false, &[(4, 12)]).map(|(p, c)| (p, cigar(&c)));
        assert_eq!(got, Some((0, "16M".into())));
    }

    #[test]
    fn divergent_ends_clipped() {
        // 40 bases mismatched every third, 100 matching, 40 mismatched every third: both ends clipped.
        let mut score = vec![CLIP_MATCH; 180];
        for i in (0..40).chain(140..180) {
            if i % 3 == 0 {
                score[i] = CLIP_MISMATCH;
            }
        }
        assert_eq!(clip_lengths(&score), (40, 39));  // base 140 matches
        // One mismatch near an end is HiFi's; a long insertion (a structural variant) is free.
        let mut one = vec![CLIP_MATCH; 180];
        one[3] = CLIP_MISMATCH;
        assert_eq!(clip_lengths(&one), (0, 0));
        let mut sv = vec![CLIP_MATCH; 180];
        sv[19] += clip_gap(80);
        assert_eq!(clip_lengths(&sv), (0, 0));
        // The CIGAR: clipped bases soft-clipped, the start past the clipped reference, a gap at the cut joins the clip.
        assert_eq!(clip_cigar(100, &[(M, 10), (I, 2), (M, 10), (D, 3), (M, 20)], 10, 5), (110, vec![(S, 12), (M, 10), (D, 3), (M, 15), (S, 5)]));
        assert_eq!(clip_cigar(100, &[(S, 4), (M, 30)], 6, 0), (102, vec![(S, 6), (M, 28)]));
    }

    #[test]
    fn graph_placement() {
        let r: Vec<u8> = (0..400).map(|i| b"ACGT"[(i * 7 + i / 3) % 4]).collect();
        let at = |a: u64, b: u64| r[a as usize..b as usize].to_vec();
        let len = r.len() as u64;
        // Read = reference 0..100, 80 letters of another copy in place of reference 100..180, then 180..300.
        let mut read = r[..100].to_vec();
        read.extend(std::iter::repeat(b'T').take(80));
        read.extend_from_slice(&r[180..300]);
        let blocks = [(0usize, 0u64, 100usize), (180, 180, 120)];
        let got = linearize(&read, &blocks, &at, len, true).map(|(p, c)| (p, cigar(&c)));
        assert_eq!(got, Some((0, "100M80I80D120M".into())));
        // Without it the equal-length stretch is aligned base for base.
        let got = linearize(&read, &blocks, &at, len, false).map(|(p, c)| (p, cigar(&c)));
        assert_eq!(got, Some((0, "300M".into())));
        // A read end of 60 letters beyond the last block is clipped, even where it matches the reference.
        let got = linearize(&r[..260], &[(0, 0, 200)], &at, len, true).map(|(p, c)| (p, cigar(&c)));
        assert_eq!(got, Some((0, "200M60S".into())));
        let got = linearize(&r[..260], &[(0, 0, 200)], &at, len, false).map(|(p, c)| (p, cigar(&c)));
        assert_eq!(got, Some((0, "260M".into())));
    }
}
