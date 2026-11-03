//! Snarls files written by `vg snarls`: a libvgio stream of `vg.Snarl` messages.
//!
//! The stream is optionally BGZF/gzip compressed and is a sequence of groups,
//! each a varint count, a varint-length type tag ("SNARL") and that many minus
//! one varint-length messages. Only the fields needed here are decoded: the
//! start and end visits and the parent's start and end.

use gbz::support::Orientation;
use std::io::Read;

/// A node visit: node id and orientation.
pub type Handle = (usize, Orientation);

pub fn flip(h: Handle) -> Handle {
    (h.0, h.1.flip())
}

pub struct Snarl {
    pub start: Handle,
    pub end: Handle,
    /// The parent snarl's (start, end), as the file records it.
    pub parent: Option<(Handle, Handle)>,
}

/// A snarl read in either direction; files often record a parent reversed.
pub fn canonical(k: (Handle, Handle)) -> (Handle, Handle) {
    let r = (flip(k.1), flip(k.0));
    if r < k { r } else { k }
}

/// `>12>34`: the snarl's id as vg writes it.
pub fn snarl_name(start: Handle, end: Handle) -> String {
    let o = |h: Handle| if h.1 == Orientation::Reverse { '<' } else { '>' };
    format!("{}{}{}{}", o(start), start.0, o(end), end.0)
}

fn bad(what: &str) -> std::io::Error {
    std::io::Error::new(std::io::ErrorKind::InvalidData, format!("malformed snarls file: {what}"))
}

fn varint(b: &[u8], i: &mut usize) -> std::io::Result<u64> {
    let (mut shift, mut v) = (0, 0u64);
    loop {
        let &c = b.get(*i).ok_or_else(|| bad("truncated varint"))?;
        *i += 1;
        v |= ((c & 0x7f) as u64) << shift;
        if c < 0x80 {
            return Ok(v);
        }
        shift += 7;
        if shift > 63 {
            return Err(bad("varint too long"));
        }
    }
}

/// The `len` bytes at `*i`, advancing past them.
fn take<'a>(b: &'a [u8], i: &mut usize, len: usize) -> std::io::Result<&'a [u8]> {
    let s = b.get(*i..*i + len).ok_or_else(|| bad("truncated message"))?;
    *i += len;
    Ok(s)
}

/// (field number, varint value, bytes) of each field of a message.
fn fields(m: &[u8]) -> std::io::Result<Vec<(u64, u64, &[u8])>> {
    let mut out = Vec::new();
    let mut i = 0;
    while i < m.len() {
        let key = varint(m, &mut i)?;
        match key & 7 {
            0 => out.push((key >> 3, varint(m, &mut i)?, &m[..0])),
            2 => {
                let len = varint(m, &mut i)? as usize;
                out.push((key >> 3, 0, take(m, &mut i, len)?));
            }
            1 => drop(take(m, &mut i, 8)?),
            5 => drop(take(m, &mut i, 4)?),
            wt => return Err(bad(&format!("protobuf wire type {wt}"))),
        }
    }
    Ok(out)
}

fn visit(m: &[u8]) -> std::io::Result<Handle> {
    let (mut id, mut back) = (0usize, false);
    for (f, v, _) in fields(m)? {
        match f {
            1 => id = v as usize,
            3 => back = v != 0,
            _ => {}
        }
    }
    Ok((id, if back { Orientation::Reverse } else { Orientation::Forward }))
}

/// A Snarl message's start and end visits (fields 2 and 3) and its parent's
/// (field 4, itself a Snarl holding only those).
fn snarl_message(m: &[u8]) -> std::io::Result<Snarl> {
    let (mut start, mut end, mut parent) = (None, None, None);
    for (f, _, b) in fields(m)? {
        match f {
            2 => start = Some(visit(b)?),
            3 => end = Some(visit(b)?),
            4 => {
                let p = snarl_message(b)?;
                parent = Some((p.start, p.end));
            }
            _ => {}
        }
    }
    Ok(Snarl { start: start.ok_or_else(|| bad("snarl without start"))?, end: end.ok_or_else(|| bad("snarl without end"))?, parent })
}

pub fn read_snarls(path: &str) -> std::io::Result<Vec<Snarl>> {
    let mut raw = Vec::new();
    std::fs::File::open(path)?.read_to_end(&mut raw)?;
    if raw.starts_with(&[0x1f, 0x8b]) {
        let mut v = Vec::new();
        flate2::read::MultiGzDecoder::new(&raw[..]).read_to_end(&mut v)?;
        raw = v;
    }
    let mut out = Vec::new();
    let mut i = 0;
    while i < raw.len() {
        let count = varint(&raw, &mut i)?;
        let tag_len = varint(&raw, &mut i)? as usize;
        let is_snarl = take(&raw, &mut i, tag_len)?.eq_ignore_ascii_case(b"SNARL");
        for _ in 1..count {
            let len = varint(&raw, &mut i)? as usize;
            let m = take(&raw, &mut i, len)?;
            if is_snarl {
                out.push(snarl_message(m)?);
            }
        }
    }
    Ok(out)
}
