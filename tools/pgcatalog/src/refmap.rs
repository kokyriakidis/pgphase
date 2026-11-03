//! The reference map of a graph: where each node lies on the reference sample's paths (each placement: the path,
//! the step along it, the reference position and the strand), with each node's length and the contigs' names. Saved
//! once per graph by `pgcatalog ref-map`, for dromos (github.com/kokyriakidis/dromos), which places alignments on
//! the reference by it (`dromos index`) without loading the graph.

use gbz::support::Orientation;
use gbz::GBZ;
use std::collections::HashMap;
use std::io::BufWriter;

/// `pgcatalog ref-map`: the reference map of a graph's reference sample, saved for dromos.
pub fn save_ref_map(graph: &str, sample: Option<&str>, out: &str) -> Result<(), String> {
    let g: GBZ = simple_sds::serialize::load_from(graph).map_err(|e| format!("{graph}: {e}"))?;
    let sample = crate::graph::resolve_reference_sample(graph, sample)?;
    RefMap::build(&g, &sample)?.save(&sample, out)
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

    /// Writes the map to `path` (a ref map: `pgcatalog ref-map`), for dromos to load instead of the graph. Little-endian:
    /// "PGCATRM1"; the reference sample; node_len and first (u32 each, per node id); the placements (path, step,
    /// start, reverse); the further placements of nodes visited twice (node, then its list); path_contig; contigs.
    fn save(&self, sample: &str, path: &str) -> Result<(), String> {
        use std::io::Write;
        let io = |e: std::io::Error| format!("{path}: {e}");
        let mut w = BufWriter::with_capacity(1 << 22, std::fs::File::create(path).map_err(io)?);
        let mut put = |b: &[u8]| w.write_all(b).map_err(io);
        let text = |s: &str| [(s.len() as u64).to_le_bytes().to_vec(), s.as_bytes().to_vec()].concat();
        put(b"PGCATRM1")?;
        put(&text(sample))?;
        put(&(self.node_len.len() as u64).to_le_bytes())?;
        for v in [&self.node_len, &self.first] {
            let mut buf = Vec::with_capacity(4 * v.len());
            for x in v.iter() {
                buf.extend_from_slice(&x.to_le_bytes());
            }
            put(&buf)?;
        }
        put(&(self.placements.len() as u64).to_le_bytes())?;
        let mut buf = Vec::with_capacity(17 * self.placements.len());
        for p in &self.placements {
            buf.extend_from_slice(&p.path.to_le_bytes());
            buf.extend_from_slice(&p.step.to_le_bytes());
            buf.extend_from_slice(&p.start.to_le_bytes());
            buf.push(p.reverse as u8);
        }
        put(&buf)?;
        let mut more: Vec<(&usize, &Vec<u32>)> = self.more.iter().collect();
        more.sort_unstable();
        put(&(more.len() as u64).to_le_bytes())?;
        for (node, ks) in more {
            put(&(*node as u64).to_le_bytes())?;
            put(&(ks.len() as u64).to_le_bytes())?;
            for k in ks {
                put(&k.to_le_bytes())?;
            }
        }
        put(&(self.path_contig.len() as u64).to_le_bytes())?;
        for c in &self.path_contig {
            put(&c.to_le_bytes())?;
        }
        put(&(self.contigs.len() as u64).to_le_bytes())?;
        for c in &self.contigs {
            put(&text(c))?;
        }
        drop(put);
        w.flush().map_err(io)
    }

}
