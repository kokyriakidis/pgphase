//! `pgcatalog node-table`: a GBZ graph's node lengths and letters as one flat file, the node table dromos
//! (github.com/kokyriakidis/dromos) reads GAM against without the graph. Little-endian, each section at a multiple
//! of 8 bytes: "DROMOSN1"; u64 count (node ids 0 .. count - 1, an id the graph lacks of length 0); u32 lengths[count]
//! and zero padding; u64 1 if the letters follow, else 0; then u64 offsets[count + 1] and the letters, node i's at
//! offsets[i] .. offsets[i + 1] on its forward strand.

use gbz::GBZ;
use std::io::{BufWriter, Write};

pub fn write(g: &GBZ, path: &str, letters: bool) -> Result<(), String> {
    let io = |e: std::io::Error| format!("{path}: {e}");
    let count = g.max_node() + 1;
    let length = |id: usize| if g.has_node(id) { g.sequence_len(id).unwrap_or(0) } else { 0 };
    let mut w = BufWriter::with_capacity(1 << 22, std::fs::File::create(path).map_err(io)?);
    w.write_all(b"DROMOSN1").map_err(io)?;
    w.write_all(&(count as u64).to_le_bytes()).map_err(io)?;
    for id in 0..count {
        let n = length(id);
        if n > u32::MAX as usize {
            return Err(format!("node {id}: {n} letters, more than a node table holds"));
        }
        w.write_all(&(n as u32).to_le_bytes()).map_err(io)?;
    }
    if count % 2 == 1 {
        w.write_all(&[0u8; 4]).map_err(io)?;
    }
    w.write_all(&(letters as u64).to_le_bytes()).map_err(io)?;
    if letters {
        let mut at = 0u64;
        for id in 0..count {
            w.write_all(&at.to_le_bytes()).map_err(io)?;
            at += length(id) as u64;
        }
        w.write_all(&at.to_le_bytes()).map_err(io)?;
        for id in 0..count {
            if let Some(s) = g.has_node(id).then(|| g.sequence(id)).flatten() {
                w.write_all(s).map_err(io)?;
            }
        }
    }
    w.flush().map_err(io)
}
