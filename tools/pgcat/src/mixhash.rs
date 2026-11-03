//! A fast hasher for integer keys (node handles, minimizer values): a multiply and rotate per word, finished by
//! MurmurHash3's mix. The standard hasher (SipHash, keyed against attacker-chosen keys) costs more than the lookups
//! it serves here, and these keys are not chosen by anyone.

use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};

#[derive(Default)]
pub struct MixHasher(u64);

const K: u64 = 0x9E37_79B9_7F4A_7C15;

impl Hasher for MixHasher {
    fn write(&mut self, bytes: &[u8]) {
        for &b in bytes {
            self.write_u64(b as u64);
        }
    }
    fn write_u8(&mut self, x: u8) {
        self.write_u64(x as u64);
    }
    fn write_u32(&mut self, x: u32) {
        self.write_u64(x as u64);
    }
    fn write_usize(&mut self, x: usize) {
        self.write_u64(x as u64);
    }
    fn write_u64(&mut self, x: u64) {
        self.0 = (self.0.rotate_left(5) ^ x).wrapping_mul(K);
    }
    fn finish(&self) -> u64 {
        let mut h = self.0;
        h ^= h >> 33;
        h = h.wrapping_mul(0xff51_afd7_ed55_8ccd);
        h ^= h >> 33;
        h = h.wrapping_mul(0xc4ce_b9fe_1a85_ec53);
        h ^ (h >> 33)
    }
}

/// A HashMap keyed by integers through MixHasher.
pub type MixMap<K, V> = HashMap<K, V, BuildHasherDefault<MixHasher>>;
