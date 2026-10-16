//! The `.pgsites` site catalog: VCF lines in independently zstd-compressed
//! position blocks, so a region query decompresses only the blocks it needs.
//!
//! Layout (integers little-endian):
//!
//! ```text
//! "PGSITES2"
//! u64 header length, zstd(VCF header text)
//! zstd(block 0) zstd(block 1) ...
//! index: u64 count, then per block u64 first_pos, last_pos, max_end, offset,
//!        compressed length, raw length, records, min_allele, max_walk
//! u64 offset of the index
//! ```
//!
//! A block holds the records whose POS falls in one BLOCK_BP window, split
//! every BLOCK_BYTES of text; a record over OWN_BLOCK_BYTES has a block of its
//! own. max_end is the largest POS + len(REF) - 1, the span tabix uses for a
//! VCF record, so a query for [beg, end] reads the blocks with first_pos <= end
//! and max_end >= beg and keeps the records that overlap. min_allele is the
//! fewest bases any allele of the block's records spans (padding base not
//! counted): a reader whose reads align along fewer graph bases can skip the
//! block unread, as no read can traverse any of its alleles whole. max_walk is
//! the most nodes any allele walk of its records visits.
//! A block of one giant record is written with a long window (log <= 31);
//! readers raise zstd's window limit to match.

use rayon::prelude::*;
use std::io::{Read, Write};

pub const MAGIC: &[u8; 8] = b"PGSITES2";
const BLOCK_BP: usize = 100_000;
const BLOCK_BYTES: usize = 4 << 20;
const OWN_BLOCK_BYTES: usize = 1 << 20;
const INDEX_FIELDS: usize = 9;

/// Blocks holding a giant record are compressed at this level at most: with the
/// long window it is within 1.5x of level 19's size at a fifth of its time
/// (chr20's 7 giant records: 6.7 MB in 5.6 s against 4.5 MB in 27 s).
const GIANT_LEVEL: i32 = 15;

/// One block as one zstd frame. A block over 16 MB is a giant record: hundreds
/// of near-identical haplotype alleles, which only a window as long as the block
/// finds; zstd's own threads keep it from holding up the others.
/// zstd parameters for a block of `len` raw bytes: past 16 MB it is a giant
/// record (hundreds of near-identical haplotype alleles), which only a window
/// as long as the block finds; zstd's own threads keep it from holding up the
/// others.
fn block_params(len: usize, level: i32, threads: usize) -> (i32, Vec<zstd::zstd_safe::CParameter>) {
    use zstd::zstd_safe::CParameter;
    if len <= 16 << 20 {
        return (level, vec![]);
    }
    let window_log = (usize::BITS - (len - 1).leading_zeros()).clamp(24, 31);
    (
        level.min(GIANT_LEVEL),
        vec![
            CParameter::EnableLongDistanceMatching(true),
            CParameter::WindowLog(window_log),
            CParameter::NbWorkers(threads as u32),
        ],
    )
}

/// One block as one zstd frame.
fn compress_block(raw: &[u8], level: i32, threads: usize) -> Vec<u8> {
    let (level, params) = block_params(raw.len(), level, threads);
    let mut c = zstd::bulk::Compressor::new(level).unwrap();
    c.include_contentsize(true).unwrap();
    for p in params {
        c.set_parameter(p).unwrap();
    }
    c.compress(raw).unwrap()
}

/// Counts the bytes written through it.
struct Counting<W: Write> {
    inner: W,
    n: u64,
}

impl<W: Write> Write for Counting<W> {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        let k = self.inner.write(b)?;
        self.n += k as u64;
        Ok(k)
    }
    fn flush(&mut self) -> std::io::Result<()> {
        self.inner.flush()
    }
}

struct Pending {
    window: usize,
    raw: Vec<u8>,
    first: usize,
    last: usize,
    max_end: usize,
    records: usize,
    min_allele: usize,
    max_walk: usize,
    /// Takes no further records.
    sealed: bool,
}

/// Writes `<path>.tmp` and renames it to `path` once complete, so a catalog
/// is never seen half written (an interrupted build leaves only the .tmp).
pub struct Writer {
    path: String,
    out: std::io::BufWriter<std::fs::File>,
    offset: u64,
    level: i32,
    threads: usize,
    pending: Vec<Pending>,
    index: Vec<[u64; INDEX_FIELDS]>,
}

impl Writer {
    pub fn create(path: &str, header: &str, level: i32, threads: usize) -> std::io::Result<Self> {
        let mut out = std::io::BufWriter::new(std::fs::File::create(format!("{path}.tmp"))?);
        out.write_all(MAGIC)?;
        let h = zstd::bulk::compress(header.as_bytes(), level)?;
        out.write_all(&(h.len() as u64).to_le_bytes())?;
        out.write_all(&h)?;
        Ok(Writer { path: path.to_string(), out, offset: 16 + h.len() as u64, level, threads, pending: Vec::new(), index: Vec::new() })
    }

    /// Records must arrive in POS order; `line` without its newline.
    pub fn add(&mut self, pos: usize, max_end: usize, min_allele: usize, max_walk: usize, line: &[u8]) {
        let window = pos / BLOCK_BP;
        let own = line.len() > OWN_BLOCK_BYTES;
        let new_block = own
            || self.pending.last().map(|b| b.sealed || b.window != window || b.raw.len() >= BLOCK_BYTES).unwrap_or(true);
        if new_block {
            self.pending.push(Pending {
                window,
                raw: Vec::new(),
                first: pos,
                last: pos,
                max_end: 0,
                records: 0,
                min_allele: usize::MAX,
                max_walk: 0,
                sealed: own,
            });
        }
        let b = self.pending.last_mut().unwrap();
        b.raw.extend_from_slice(line);
        b.raw.push(b'\n');
        b.last = pos;
        b.max_end = b.max_end.max(max_end);
        b.min_allele = b.min_allele.min(min_allele);
        b.max_walk = b.max_walk.max(max_walk);
        b.records += 1;
    }

    /// One record of `len` bytes (newline excluded) in a block of its own,
    /// written by `write` straight into the compressor: the line never exists
    /// in memory whole.
    pub fn add_streamed(
        &mut self,
        pos: usize,
        max_end: usize,
        min_allele: usize,
        max_walk: usize,
        len: usize,
        write: impl FnOnce(&mut dyn Write) -> std::io::Result<()>,
    ) -> std::io::Result<()> {
        self.write_blocks(self.pending.len())?;
        let raw_len = len + 1;
        let (level, params) = block_params(raw_len, self.level, self.threads);
        let mut counting = Counting { inner: &mut self.out, n: 0 };
        {
            let mut enc = zstd::stream::write::Encoder::new(&mut counting, level)?;
            enc.include_contentsize(true)?;
            enc.set_pledged_src_size(Some(raw_len as u64))?;
            for p in params {
                enc.set_parameter(p)?;
            }
            let mut buf = std::io::BufWriter::with_capacity(1 << 20, Counting { inner: enc, n: 0 });
            write(&mut buf)?;
            buf.write_all(b"\n")?;
            let inner = buf.into_inner().map_err(|e| e.into_error())?;
            if inner.n != raw_len as u64 {
                return Err(std::io::Error::other(format!("record at {pos}: wrote {} bytes, expected {raw_len}", inner.n)));
            }
            inner.inner.finish()?;
        }
        let comp = counting.n;
        self.index.push([pos as u64, pos as u64, max_end as u64, self.offset, comp, raw_len as u64, 1, min_allele as u64, max_walk as u64]);
        self.offset += comp;
        Ok(())
    }

    /// Compress and write every block but the one still filling, in parallel.
    pub fn flush(&mut self) -> std::io::Result<()> {
        self.write_blocks(self.pending.len().saturating_sub(1))
    }

    fn write_blocks(&mut self, n: usize) -> std::io::Result<()> {
        let blocks: Vec<Pending> = self.pending.drain(..n).collect();
        let (level, threads) = (self.level, self.threads);
        let comp: Vec<Vec<u8>> = blocks.par_iter().map(|b| compress_block(&b.raw, level, threads)).collect();
        for (b, c) in blocks.iter().zip(comp) {
            self.out.write_all(&c)?;
            self.index.push([
                b.first as u64, b.last as u64, b.max_end as u64, self.offset,
                c.len() as u64, b.raw.len() as u64, b.records as u64, b.min_allele as u64, b.max_walk as u64,
            ]);
            self.offset += c.len() as u64;
        }
        Ok(())
    }

    /// Writes the remaining blocks and the index; returns the block count.
    pub fn finish(mut self) -> std::io::Result<usize> {
        self.write_blocks(self.pending.len())?;
        let index_offset = self.offset;
        self.out.write_all(&(self.index.len() as u64).to_le_bytes())?;
        for e in &self.index {
            for v in e {
                self.out.write_all(&v.to_le_bytes())?;
            }
        }
        self.out.write_all(&index_offset.to_le_bytes())?;
        self.out.into_inner().map_err(|e| e.into_error())?.sync_all()?;
        std::fs::rename(format!("{}.tmp", self.path), &self.path)?;
        Ok(self.index.len())
    }
}

pub struct Reader {
    raw: Vec<u8>,
    pub header: String,
    /// The CHROM of every record (from the header's ##contig line).
    pub chrom: String,
    /// (first_pos, max_end, offset, compressed length, raw length) per block.
    blocks: Vec<[usize; 5]>,
}

impl Reader {
    pub fn open(path: &str) -> std::io::Result<Self> {
        let mut raw = Vec::new();
        std::fs::File::open(path)?.read_to_end(&mut raw)?;
        let bad = |m: &str| std::io::Error::new(std::io::ErrorKind::InvalidData, format!("{path}: {m}"));
        if raw.starts_with(b"PGSITES") && !raw.starts_with(MAGIC) {
            return Err(bad("an older .pgsites format; rebuild it with this pgcat"));
        }
        if !raw.starts_with(MAGIC) || raw.len() < 24 {
            return Err(bad("not a .pgsites file"));
        }
        let u64_at = |o: usize| u64::from_le_bytes(raw[o..o + 8].try_into().unwrap()) as usize;
        let header_len = u64_at(8);
        let header = String::from_utf8(decompress(&raw[16..16 + header_len], None)?).map_err(|_| bad("header"))?;
        let at = header.find("##contig=<ID=").ok_or_else(|| bad("no ##contig line"))? + 13;
        let chrom = header[at..at + header[at..].find([',', '>']).unwrap_or(0)].to_string();
        let index_offset = u64_at(raw.len() - 8);
        let n = u64_at(index_offset);
        let blocks = (0..n)
            .map(|k| {
                let e = |j: usize| u64_at(index_offset + 8 + (k * INDEX_FIELDS + j) * 8);
                [e(0), e(2), e(3), e(4), e(5)]
            })
            .collect();
        Ok(Reader { raw, header, chrom, blocks })
    }

    /// Calls `f` with each record line overlapping [beg, end] (1-based,
    /// inclusive), in order, as `tabix` would return them from the VCF.
    pub fn for_each_record(&self, beg: usize, end: usize, mut f: impl FnMut(&[u8])) -> std::io::Result<()> {
        for &[first, max_end, offset, comp_len, raw_len] in &self.blocks {
            if first > end || max_end < beg {
                continue;
            }
            let block = decompress(&self.raw[offset..offset + comp_len], Some(raw_len))?;
            for line in block.split(|&c| c == b'\n').filter(|l| !l.is_empty()) {
                let mut it = line.split(|&c| c == b'\t');
                it.next();
                let pos: usize = std::str::from_utf8(it.next().unwrap_or_default()).ok()
                    .and_then(|p| p.parse().ok()).unwrap_or(0);
                it.next();
                let ref_len = it.next().map(|r| r.len()).unwrap_or(1).max(1);
                if pos <= end && pos + ref_len > beg {
                    f(line);
                }
            }
        }
        Ok(())
    }
}

fn decompress(frame: &[u8], raw_len: Option<usize>) -> std::io::Result<Vec<u8>> {
    let mut d = zstd::bulk::Decompressor::new()?;
    d.set_parameter(zstd::zstd_safe::DParameter::WindowLogMax(31))?;
    d.decompress(frame, raw_len.unwrap_or(1 << 24))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Records across several windows, one with a REF reaching far past its
    /// window: every region query must return exactly the overlapping records.
    #[test]
    fn region_queries_return_the_overlapping_records() {
        let path = std::env::temp_dir().join(format!("pgcat_test_{}.pgsites", std::process::id()));
        let path = path.to_str().unwrap();
        let header = "##fileformat=VCFv4.2\n##contig=<ID=CHM13#0#chrT,length=1000000>\n#CHROM\tPOS\tID\tREF\tALT\n";
        let mut records: Vec<(usize, String)> = Vec::new();
        for k in 0..400 {
            let pos = 1 + k * 1_337;
            let reference = if k == 7 { "A".repeat(250_000) } else { "AC".to_string() };
            records.push((pos, reference));
        }
        let mut w = Writer::create(path, header, 3, 2).unwrap();
        for (k, (pos, reference)) in records.iter().enumerate() {
            let line = format!("CHM13#0#chrT\t{pos}\t>{k}>{}\t{reference}\tG", k + 1);
            w.add(*pos, pos + reference.len() - 1, reference.len(), 3, line.as_bytes());
            if k % 50 == 0 {
                w.flush().unwrap();
            }
        }
        assert!(w.finish().unwrap() > 1);

        let r = Reader::open(path).unwrap();
        assert_eq!(r.chrom, "CHM13#0#chrT");
        assert_eq!(r.header, header);
        for &(beg, end) in &[(1, 1), (2, 2), (3, 3), (100_000, 100_001), (200_000, 260_000), (1, usize::MAX), (600_000, 700_000)] {
            let mut got = Vec::new();
            r.for_each_record(beg, end, |l| got.push(String::from_utf8(l.to_vec()).unwrap())).unwrap();
            let want: Vec<usize> = records.iter().filter(|(p, s)| *p <= end && p + s.len() > beg).map(|(p, _)| *p).collect();
            let got_pos: Vec<usize> = got.iter().map(|l| l.split('\t').nth(1).unwrap().parse().unwrap()).collect();
            assert_eq!(got_pos, want, "region {beg}-{end}");
        }
        std::fs::remove_file(path).unwrap();
    }
}
