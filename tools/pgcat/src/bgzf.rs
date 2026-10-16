//! BGZF output compressed in parallel, with a tabix index built as it goes.
//!
//! The stream is cut into blocks of exactly BLOCK bytes, so the virtual offset
//! of any byte follows from its stream offset and the blocks' compressed
//! addresses. Records are indexed once the block holding their end is written.

use crate::tabix::TbiBuilder;
use flate2::{Compress, Compression, Crc, FlushCompress, Status};
use rayon::prelude::*;
use std::collections::VecDeque;
use std::io::{BufWriter, Write};

/// Uncompressed bytes per block, as bgzip and htslib use.
pub const BLOCK: usize = 0xff00;
/// Blocks compressed per parallel batch.
const BATCH_BLOCKS: usize = 256;
const MAX_BLOCK_SIZE: usize = 0x10000;

/// The empty block that ends every BGZF file.
pub const EOF_BLOCK: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0, 0x42, 0x43, 0x02, 0, 0x1b, 0, 0x03, 0, 0, 0, 0, 0, 0, 0, 0, 0,
];

/// One BGZF block: gzip header with the BC extra field, raw deflate, CRC32, ISIZE.
pub fn compress_block(data: &[u8], level: Compression) -> Vec<u8> {
    let deflate = |level: Compression| -> Option<Vec<u8>> {
        let mut c = Compress::new(level, false);
        let mut out = Vec::with_capacity(MAX_BLOCK_SIZE);
        match c.compress_vec(data, &mut out, FlushCompress::Finish) {
            Ok(Status::StreamEnd) if out.len() + 26 <= MAX_BLOCK_SIZE => Some(out),
            _ => None,
        }
    };
    // Incompressible data is stored, as htslib does, to stay within a block.
    let body = deflate(level).or_else(|| deflate(Compression::none())).expect("BGZF block overflow");
    let mut crc = Crc::new();
    crc.update(data);
    let size = body.len() + 26;
    let mut block = Vec::with_capacity(size);
    block.extend_from_slice(&[0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0, 0x42, 0x43, 0x02, 0]);
    block.extend_from_slice(&((size - 1) as u16).to_le_bytes());
    block.extend_from_slice(&body);
    block.extend_from_slice(&crc.sum().to_le_bytes());
    block.extend_from_slice(&(data.len() as u32).to_le_bytes());
    block
}

/// A record waiting for the block holding its end.
struct Pending {
    contig: u32,
    beg: u64,
    end: u64,
    stream_end: u64,
}

pub struct IndexedWriter {
    out: BufWriter<std::fs::File>,
    level: Compression,
    buf: Vec<u8>,
    /// Stream offset of buf[0]: always a block boundary.
    buf_start: u64,
    /// Compressed bytes written so far: the address of the next block.
    written: u64,
    /// Compressed address of each block from the one holding `pending[0]`'s end.
    addrs: VecDeque<(u64, u64)>, // (block number, address)
    pending: VecDeque<Pending>,
    index: TbiBuilder,
    contigs: Vec<String>,
}

impl IndexedWriter {
    /// `contigs` names the records' contig numbers; `level` 0-9.
    pub fn create(path: &str, contigs: Vec<String>, level: u32) -> std::io::Result<Self> {
        Ok(IndexedWriter {
            out: BufWriter::with_capacity(1 << 22, std::fs::File::create(path)?),
            level: Compression::new(level),
            buf: Vec::with_capacity(BLOCK * BATCH_BLOCKS * 2),
            buf_start: 0,
            written: 0,
            addrs: VecDeque::new(),
            pending: VecDeque::new(),
            index: TbiBuilder::new(),
            contigs,
        })
    }

    /// Writes header lines (each starting with '#', with its newline), before
    /// any record: tabix skips them, and the index starts after them.
    pub fn header(&mut self, text: &[u8]) -> std::io::Result<()> {
        assert!(self.pending.is_empty() && self.buf_start == 0, "header after records");
        self.buf.extend_from_slice(text);
        let start = self.buf.len() as u64;
        if start >= BLOCK as u64 {
            return Err(std::io::Error::other("header longer than a BGZF block"));
        }
        // Block 0 starts at address 0: the first record starts at (0, start).
        self.index.set_first_offset(start);
        Ok(())
    }

    /// Appends one line (with its newline) covering [beg, end) of `contig`.
    pub fn add(&mut self, line: &[u8], contig: u32, beg: u64, end: u64) -> std::io::Result<()> {
        self.buf.extend_from_slice(line);
        let stream_end = self.buf_start + self.buf.len() as u64;
        self.pending.push_back(Pending { contig, beg, end, stream_end });
        if self.buf.len() >= BLOCK * BATCH_BLOCKS {
            self.write_blocks(self.buf.len() / BLOCK)?;
        }
        Ok(())
    }

    /// Compresses and writes the first `n` blocks of buf (all, if the last is
    /// partial, when finishing), then indexes the records they complete.
    fn write_blocks(&mut self, n: usize) -> std::io::Result<()> {
        let take = (n * BLOCK).min(self.buf.len());
        let level = self.level;
        let blocks: Vec<Vec<u8>> = self.buf[..take].par_chunks(BLOCK).map(|b| compress_block(b, level)).collect();
        let first_block = self.buf_start / BLOCK as u64;
        for (k, b) in blocks.iter().enumerate() {
            self.addrs.push_back((first_block + k as u64, self.written));
            self.out.write_all(b)?;
            self.written += b.len() as u64;
        }
        self.buf.drain(..take);
        self.buf_start += take as u64;
        self.index_written(false)
    }

    /// The virtual offset of stream offset `s`, which must lie in a written
    /// block or at the start of the next. A block's end reads as the start of
    /// the next block, as htslib reports it after reading a line.
    fn voffset(&self, s: u64) -> u64 {
        if s == self.buf_start {
            return self.written << 16;
        }
        let block = s / BLOCK as u64;
        let within = s % BLOCK as u64;
        let first = self.addrs.front().expect("no block for the offset").0;
        let (_, addr) = self.addrs[(block - first) as usize];
        (addr << 16) | within
    }

    fn index_written(&mut self, all: bool) -> std::io::Result<()> {
        while let Some(p) = self.pending.front() {
            if !all && p.stream_end > self.buf_start {
                break;
            }
            let v = self.voffset(p.stream_end);
            let name = &self.contigs[p.contig as usize];
            self.index.push(name, p.beg as i64, p.end as i64, v).map_err(std::io::Error::other)?;
            self.pending.pop_front();
        }
        // Keep only the block addresses later records may need.
        let keep_from = self.pending.front().map(|p| p.stream_end / BLOCK as u64).unwrap_or(u64::MAX);
        while self.addrs.len() > 1 && self.addrs.front().map(|a| a.0 < keep_from).unwrap_or(false) {
            self.addrs.pop_front();
        }
        Ok(())
    }

    /// Writes the remaining blocks and the EOF block, and the index to `tbi`.
    pub fn finish(mut self, tbi: &str) -> std::io::Result<()> {
        let n = self.buf.len().div_ceil(BLOCK);
        if n > 0 {
            self.write_blocks(n)?;
        }
        self.index_written(true)?;
        let final_offset = self.written << 16;
        self.out.write_all(&EOF_BLOCK)?;
        self.out.into_inner().map_err(|e| e.into_error())?.sync_all()?;
        self.index.finish(final_offset);
        self.index.write(tbi)
    }
}
