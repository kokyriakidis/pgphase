//! A tabix (.tbi) index, built the way htslib's `tabix` builds one for
//! `-s 1 -b 2 -e 3` (hts_idx_push, hts_idx_finish, compress_binning), so
//! htslib reads it exactly as its own. Bins are written in numeric order
//! (htslib writes them in hash order; readers do not depend on it).

use crate::bgzf::{compress_block, BLOCK, EOF_BLOCK};
use std::collections::{BTreeMap, HashMap};
use std::io::Write;

const MIN_SHIFT: u32 = 14;
const N_LVLS: u32 = 5;
const N_BINS: u32 = ((1 << (3 * N_LVLS + 3)) - 1) / 7;
/// The pseudo-bin holding each contig's offsets and record counts.
const META_BIN: u32 = N_BINS + 1;
const MIN_MARKER_DIST: u64 = 0x10000;
const NONE: u64 = u64::MAX;

fn reg2bin(beg: i64, end: i64) -> u32 {
    let end = end - 1;
    let mut s = MIN_SHIFT;
    let mut t: i64 = ((1 << (3 * N_LVLS)) - 1) / 7;
    for l in (1..=N_LVLS).rev() {
        if beg >> s == end >> s {
            return (t + (beg >> s)) as u32;
        }
        // htslib's loop decrements the level before this step.
        s += 3;
        t -= 1 << (3 * (l - 1));
    }
    0
}

fn bin_first(level: u32) -> u32 {
    ((1 << (3 * level)) - 1) / 7
}

fn bin_parent(bin: u32) -> u32 {
    (bin - 1) >> 3
}

fn bin_level(bin: u32) -> u32 {
    let (mut l, mut b) = (0, bin);
    while b != 0 {
        b = bin_parent(b);
        l += 1;
    }
    l
}

#[derive(Default)]
struct Contig {
    bins: BTreeMap<u32, Vec<(u64, u64)>>,
    linear: Vec<u64>,
}

pub struct TbiBuilder {
    names: Vec<String>,
    ids: HashMap<String, usize>,
    contigs: Vec<Contig>,
    // hts_idx_t.z
    save_tid: i64,
    last_tid: i64,
    save_bin: u32,
    last_bin: u32,
    save_off: u64,
    last_off: u64,
    off_beg: u64,
    last_coor: i64,
    n_mapped: u64,
    n_unmapped: u64,
}

impl Default for TbiBuilder {
    fn default() -> Self {
        Self::new()
    }
}

impl TbiBuilder {
    pub fn new() -> Self {
        TbiBuilder {
            names: Vec::new(),
            ids: HashMap::new(),
            contigs: Vec::new(),
            save_tid: -1,
            last_tid: -1,
            save_bin: u32::MAX,
            last_bin: u32::MAX,
            save_off: 0,
            last_off: 0,
            off_beg: 0,
            last_coor: -1,
            n_mapped: 0,
            n_unmapped: 0,
        }
    }

    /// The virtual offset where the first record starts (after header lines).
    pub fn set_first_offset(&mut self, offset: u64) {
        self.save_off = offset;
        self.last_off = offset;
        self.off_beg = offset;
    }

    fn insert_bin(&mut self, tid: i64, bin: u32, a: u64, b: u64) {
        self.contigs[tid as usize].bins.entry(bin).or_default().push((a, b));
    }

    /// A record on `name` covering [beg, end) (0-based), whose line ends at
    /// virtual offset `offset`. Records must come contig by contig, sorted.
    pub fn push(&mut self, name: &str, beg: i64, end: i64, offset: u64) -> Result<(), String> {
        let tid = match self.ids.get(name) {
            Some(&t) => t as i64,
            None => {
                self.ids.insert(name.to_string(), self.names.len());
                self.names.push(name.to_string());
                self.contigs.push(Contig::default());
                (self.names.len() - 1) as i64
            }
        };
        if end > 1i64 << (MIN_SHIFT + 3 * N_LVLS) {
            return Err(format!("{name}:{end} is beyond what a .tbi index can hold"));
        }
        if self.last_tid != tid {
            if !self.contigs[tid as usize].bins.is_empty() {
                return Err(format!("records of {name} are not together"));
            }
            self.last_tid = tid;
            self.last_bin = u32::MAX;
        } else if self.last_coor > beg {
            return Err(format!("records of {name} are not sorted"));
        }
        if end < beg {
            return Err(format!("record on {name} ends before it begins"));
        }
        let (beg, end) = (beg.max(0), end.max(1));
        // Linear index: the first record overlapping each 16 kb window.
        let l = &mut self.contigs[tid as usize].linear;
        let (wb, we) = ((beg >> MIN_SHIFT) as usize, ((end - 1) >> MIN_SHIFT) as usize);
        if l.len() < we + 1 {
            l.resize(we + 1, NONE);
        }
        for w in &mut l[wb..=we] {
            if *w == NONE {
                *w = self.last_off;
            }
        }
        let bin = reg2bin(beg, end);
        if self.last_bin != bin {
            if self.save_bin != u32::MAX {
                self.insert_bin(self.save_tid, self.save_bin, self.save_off, self.last_off);
            }
            if self.last_bin == u32::MAX && self.save_bin != u32::MAX {
                // A new contig: close the previous one's metadata.
                let off_end = self.last_off;
                self.insert_bin(self.save_tid, META_BIN, self.off_beg, off_end);
                self.insert_bin(self.save_tid, META_BIN, self.n_mapped, self.n_unmapped);
                self.n_mapped = 0;
                self.n_unmapped = 0;
                self.off_beg = off_end;
            }
            self.save_off = self.last_off;
            self.save_bin = bin;
            self.last_bin = bin;
            self.save_tid = tid;
        }
        self.n_mapped += 1;
        self.last_off = offset;
        self.last_coor = beg;
        Ok(())
    }

    /// `final_offset`: the virtual offset after the last record.
    pub fn finish(&mut self, final_offset: u64) {
        if self.save_tid >= 0 {
            self.insert_bin(self.save_tid, self.save_bin, self.save_off, final_offset);
            self.insert_bin(self.save_tid, META_BIN, self.off_beg, final_offset);
            self.insert_bin(self.save_tid, META_BIN, self.n_mapped, self.n_unmapped);
        }
        for c in &mut self.contigs {
            // Windows no record starts in point to the next one that does.
            for i in (0..c.linear.len().saturating_sub(1)).rev() {
                if c.linear[i] == NONE {
                    c.linear[i] = c.linear[i + 1];
                }
            }
            compress_binning(&mut c.bins);
        }
    }

    /// The index as a BGZF file.
    pub fn write(&self, path: &str) -> std::io::Result<()> {
        let mut raw: Vec<u8> = Vec::new();
        let i32le = |v: usize| (v as i32).to_le_bytes();
        raw.extend_from_slice(b"TBI\x01");
        raw.extend_from_slice(&i32le(self.names.len()));
        // tbx_conf_t: TBX_GENERIC, sequence column 1, begin 2, end 3, meta '#', skip 0.
        for v in [0i32, 1, 2, 3, b'#' as i32, 0] {
            raw.extend_from_slice(&v.to_le_bytes());
        }
        let names_len: usize = self.names.iter().map(|n| n.len() + 1).sum();
        raw.extend_from_slice(&i32le(names_len));
        for n in &self.names {
            raw.extend_from_slice(n.as_bytes());
            raw.push(0);
        }
        for c in &self.contigs {
            raw.extend_from_slice(&i32le(c.bins.len()));
            for (bin, chunks) in &c.bins {
                raw.extend_from_slice(&bin.to_le_bytes());
                raw.extend_from_slice(&i32le(chunks.len()));
                for (a, b) in chunks {
                    raw.extend_from_slice(&a.to_le_bytes());
                    raw.extend_from_slice(&b.to_le_bytes());
                }
            }
            raw.extend_from_slice(&i32le(c.linear.len()));
            for o in &c.linear {
                raw.extend_from_slice(&o.to_le_bytes());
            }
        }
        raw.extend_from_slice(&0u64.to_le_bytes()); // records without coordinates
        let tmp = format!("{path}.tmp");
        let mut f = std::io::BufWriter::new(std::fs::File::create(&tmp)?);
        for block in raw.chunks(BLOCK) {
            f.write_all(&compress_block(block, flate2::Compression::default()))?;
        }
        f.write_all(&EOF_BLOCK)?;
        f.into_inner().map_err(|e| e.into_error())?.sync_all()?;
        std::fs::rename(&tmp, path)
    }
}

/// htslib's compress_binning: a bin whose chunks span under 64 kb of
/// compressed file joins its parent, deepest level first; then chunks that
/// start in the same BGZF block as the previous one ends are merged.
fn compress_binning(bins: &mut BTreeMap<u32, Vec<(u64, u64)>>) {
    for level in (1..=N_LVLS).rev() {
        let start = bin_first(level);
        let at_level: Vec<u32> = bins.keys().copied().filter(|&b| b >= start && b < N_BINS && bin_level(b) == level).collect();
        for b in at_level {
            let list = bins.get_mut(&b).unwrap();
            list.sort_by_key(|c| c.0);
            let span = (list[list.len() - 1].1 >> 16) - (list[0].0 >> 16);
            let parent = bin_parent(b);
            if span < MIN_MARKER_DIST && bins.contains_key(&parent) {
                let moved = bins.remove(&b).unwrap();
                bins.get_mut(&parent).unwrap().extend(moved);
            }
        }
    }
    for (&bin, list) in bins.iter_mut() {
        if bin >= N_BINS {
            continue;
        }
        list.sort_by_key(|c| c.0);
        let mut m = 0;
        for l in 1..list.len() {
            if list[m].1 >> 16 >= list[l].0 >> 16 {
                if list[m].1 < list[l].1 {
                    list[m].1 = list[l].1;
                }
            } else {
                m += 1;
                list[m] = list[l];
            }
        }
        list.truncate(m + 1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bins_as_htslib() {
        assert_eq!(N_BINS, 37449);
        assert_eq!(reg2bin(0, 1), 4681);
        assert_eq!(reg2bin(0, 1 << 14), 4681);
        assert_eq!(reg2bin(0, (1 << 14) + 1), 585);
        assert_eq!(reg2bin(0, 1 << 29), 0);
        assert_eq!(bin_level(0), 0);
        assert_eq!(bin_level(4681), 5);
        assert_eq!(bin_parent(4681), 585);
    }
}
