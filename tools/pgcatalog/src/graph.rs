//! What can be learned about a GBZ without loading it.

use gbz::headers::{GBWTPayload, GBZPayload, Header};
use gbz::support::Tags;
use gbz::Metadata;
use simple_sds::serialize::Serialize;
use simple_sds::sparse_vector::SparseVector;
use std::collections::BTreeMap;
use std::io::{BufReader, Read, Seek, SeekFrom};

fn open(path: &str) -> Result<BufReader<std::fs::File>, String> {
    let f = std::fs::File::open(path).map_err(|e| format!("{path}: {e}"))?;
    Ok(BufReader::with_capacity(1 << 20, f))
}

/// Reads the GBZ up to the end of the GBWT's tags; returns them.
fn gbwt_tags(r: &mut BufReader<std::fs::File>, path: &str) -> Result<Tags, String> {
    let err = |what: &str, e: std::io::Error| format!("{path}: {what}: {e}");
    Header::<GBZPayload>::load(r).map_err(|e| err("not a GBZ file", e))?;
    Tags::load(r).map_err(|e| err("GBZ tags", e))?;
    Header::<GBWTPayload>::load(r).map_err(|e| err("GBWT header", e))?;
    Tags::load(r).map_err(|e| err("GBWT tags", e))
}

/// The reference samples a GBZ declares (the GBWT's `reference_samples` tag),
/// read from the start of the file.
pub fn reference_samples(path: &str) -> Result<Vec<String>, String> {
    let tags = gbwt_tags(&mut open(path)?, path)?;
    Ok(tags
        .get(gbz::REFERENCE_SAMPLES_KEY)
        .map(|v| v.split(' ').filter(|x| !x.is_empty()).map(String::from).collect())
        .unwrap_or_default())
}

/// The sample's phase-0 paths by contig, named as vg names them
/// (`CHM13#0#chr6[31350872]` for a fragment starting at 31350872), from the
/// GBWT's metadata. The GBWT body before it (the BWT index, then its records
/// and the document array samples, two size-prefixed arrays) is read past
/// rather than loaded: seconds, against a whole-graph load for `vg paths`.
pub fn reference_paths(path: &str, sample: &str) -> Result<BTreeMap<String, Vec<String>>, String> {
    let mut r = open(path)?;
    gbwt_tags(&mut r, path)?;
    let err = |what: &str, e: std::io::Error| format!("{path}: {what}: {e}");
    SparseVector::load(&mut r).map_err(|e| err("BWT index", e))?;
    let bytes = usize::load(&mut r).map_err(|e| err("BWT", e))?;
    r.seek(SeekFrom::Current(bytes.div_ceil(8) as i64 * 8)).map_err(|e| err("BWT", e))?;
    let words = usize::load(&mut r).map_err(|e| err("DA samples", e))?;
    r.seek(SeekFrom::Current(words as i64 * 8)).map_err(|e| err("DA samples", e))?;
    let md = Option::<Metadata>::load(&mut r)
        .map_err(|e| err("metadata", e))?
        .ok_or_else(|| format!("{path}: the GBZ has no path metadata"))?;
    let Some(sample_id) = md.sample_id(sample) else { return Ok(BTreeMap::new()) };
    let mut by_contig: BTreeMap<String, Vec<(usize, String)>> = BTreeMap::new();
    for id in 0..md.paths() {
        let Some(p) = md.path(id) else { continue };
        if p.sample() != sample_id || p.phase() != 0 {
            continue;
        }
        let contig = md.contig_name(p.contig());
        let mut name = format!("{sample}#0#{contig}");
        if p.fragment != 0 {
            name.push_str(&format!("[{}]", p.fragment));
        }
        by_contig.entry(contig).or_default().push((p.fragment as usize, name));
    }
    Ok(by_contig
        .into_iter()
        .map(|(c, mut v)| {
            v.sort();
            (c, v.into_iter().map(|(_, n)| n).collect())
        })
        .collect())
}

/// The reference sample: the one given, or the graph's only declared one.
pub fn resolve_reference_sample(path: &str, given: Option<&str>) -> Result<String, String> {
    if let Some(s) = given {
        return Ok(s.to_string());
    }
    let found = reference_samples(path)?;
    match found.len() {
        1 => Ok(found[0].clone()),
        0 => Err(format!("{path} declares no reference sample; give one with --reference-sample")),
        _ => Err(format!(
            "{path} declares several reference samples ({}); choose one with --reference-sample",
            found.join(", ")
        )),
    }
}

/// An identity for the graph's content, cheap to compute: the file size and a
/// hash of its first and last 4 MB (header, tags and index structure).
pub fn fingerprint(path: &str) -> Result<String, String> {
    const SPAN: u64 = 4 << 20;
    let mut f = std::fs::File::open(path).map_err(|e| format!("{path}: {e}"))?;
    let size = f.metadata().map_err(|e| e.to_string())?.len();
    let mut h: u64 = 0xcbf29ce484222325; // FNV-1a
    let mut feed = |buf: &[u8]| {
        for &b in buf {
            h ^= b as u64;
            h = h.wrapping_mul(0x100000001b3);
        }
    };
    let mut buf = vec![0u8; SPAN.min(size) as usize];
    f.read_exact(&mut buf).map_err(|e| e.to_string())?;
    feed(&buf);
    f.seek(SeekFrom::Start(size.saturating_sub(SPAN))).map_err(|e| e.to_string())?;
    let mut tail = vec![0u8; SPAN.min(size) as usize];
    f.read_exact(&mut tail).map_err(|e| e.to_string())?;
    feed(&tail);
    Ok(format!("size={size} fnv1a={h:016x}"))
}
