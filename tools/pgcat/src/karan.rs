//! karan (github.com/kokyriakidis/karan) through its C API: exact pairwise alignment. Linked when pgcat is built with
//! the `karan` feature, from the folder KARAN_LIB names (karan's `pixi run build-c` leaves libkaran.so and the
//! runtime libraries it loads in build/c; build.rs links it and sets the rpath).
//!
//! Each thread aligns through an aligner of its own, which keeps its memory from call to call.

use std::cell::RefCell;
use std::os::raw::c_char;

#[repr(C)]
struct Costs {
    mismatch: i64,
    opening: i64,
    extension: i64,
    opening2: i64,
    extension2: i64,
    deletion_opening: i64,
    deletion_extension: i64,
    deletion_opening2: i64,
    deletion_extension2: i64,
}

#[repr(C)]
struct Mode {
    kind: i64,
    reference_start: i64,
    reference_end: i64,
    query_start: i64,
    query_end: i64,
    match_score: i64,
    anchor: i64,
    zdrop: i64,
    end_bonus: i64,
}

#[repr(C)]
struct Options {
    band_low: i64,
    band_high: i64,
    max_cost: i64,
    eqx: i64,
    right_ties: i64,
    max_memory: i64,
}

#[repr(C)]
struct Alignment {
    cost: i64,
    score: i64,
    reference_start: i64,
    reference_end: i64,
    query_start: i64,
    query_end: i64,
    cigar: *mut c_char,
    cigar_length: i64,
}

#[repr(C)]
struct RawAligner {
    _private: [u8; 0],
}

extern "C" {
    fn karan_aligner_new() -> *mut RawAligner;
    fn karan_aligner_free(aligner: *mut RawAligner);
    fn karan_aligner_align(
        aligner: *mut RawAligner,
        reference: *const c_char,
        reference_length: i64,
        query: *const c_char,
        query_length: i64,
        costs: *const Costs,
        mode: *const Mode,
        options: *const Options,
        alignment: *mut Alignment,
    ) -> i64;
    fn karan_free(cigar: *mut c_char);
}

/// KARAN_ENDS_FREE, and KARAN_ALL: as many free letters as a sequence has.
const ENDS_FREE: i64 = 0;
const ALL: i64 = i64::MAX;

/// A thread's aligner, freed with the thread.
struct Aligner(*mut RawAligner);

impl Drop for Aligner {
    fn drop(&mut self) {
        // SAFETY: the pointer came from karan_aligner_new and is freed once, here.
        unsafe { karan_aligner_free(self.0) }
    }
}

thread_local! {
    static ALIGNER: RefCell<Option<Aligner>> = const { RefCell::new(None) };
}

/// A read's alignment inside a stretch of reference: the stretch's and the read's aligned spans, its cost, and its
/// CIGAR over them (`=`, `X`, `I`, `D`).
#[derive(Debug, PartialEq)]
pub struct Placed {
    pub reference_start: usize,
    pub reference_end: usize,
    pub query_start: usize,
    pub query_end: usize,
    pub cost: i64,
    pub cigar: Vec<(u8, u32)>,
}

/// The least-cost alignment of all of `query` inside `reference` -- the read anywhere in the stretch, every letter of
/// it aligned -- at minimap2's HiFi costs (mismatch 4, a gap of k letters 6 + 2k or 26 + k, the lesser), its moves on
/// the diagonals `band.0 ..= band.1` (reference letters less query letters from the stretch's start). Exact within the
/// band: karan finds the cheapest alignment there, ties placing indels left. None if no alignment fits the band.
pub fn place(reference: &[u8], query: &[u8], band: (i64, i64)) -> Option<Placed> {
    let costs = Costs {
        mismatch: 4,
        opening: 6,
        extension: 2,
        opening2: 26,
        extension2: 1,
        deletion_opening: 0,
        deletion_extension: 0,
        deletion_opening2: -1,
        deletion_extension2: 0,
    };
    let mode = Mode {
        kind: ENDS_FREE,
        reference_start: ALL,
        reference_end: ALL,
        query_start: 0,
        query_end: 0,
        match_score: 0,
        anchor: 0,
        zdrop: 0,
        end_bonus: 0,
    };
    let options = Options { band_low: band.0, band_high: band.1, max_cost: -1, eqx: 1, right_ties: 0, max_memory: 0 };
    let mut a = Alignment {
        cost: 0,
        score: 0,
        reference_start: 0,
        reference_end: 0,
        query_start: 0,
        query_end: 0,
        cigar: std::ptr::null_mut(),
        cigar_length: 0,
    };
    let status = ALIGNER.with(|cell| {
        let mut slot = cell.borrow_mut();
        if slot.is_none() {
            // SAFETY: a new aligner; null only when memory runs out.
            let raw = unsafe { karan_aligner_new() };
            if raw.is_null() {
                return None;
            }
            *slot = Some(Aligner(raw));
        }
        let aligner = slot.as_ref().unwrap().0;
        // SAFETY: the sequences are valid for their lengths, the structs live for the call, and this thread alone
        // uses its aligner.
        Some(unsafe {
            karan_aligner_align(
                aligner,
                reference.as_ptr() as *const c_char,
                reference.len() as i64,
                query.as_ptr() as *const c_char,
                query.len() as i64,
                &costs,
                &mode,
                &options,
                &mut a,
            )
        })
    })?;
    if status != 0 {
        return None;
    }
    // SAFETY: on success the CIGAR is cigar_length bytes, freed below with karan_free.
    let text = unsafe { std::slice::from_raw_parts(a.cigar as *const u8, a.cigar_length as usize) }.to_vec();
    unsafe { karan_free(a.cigar) };
    let mut cigar: Vec<(u8, u32)> = Vec::new();
    let mut n: u32 = 0;
    for c in text {
        if c.is_ascii_digit() {
            n = n * 10 + (c - b'0') as u32;
        } else {
            cigar.push((c, n));
            n = 0;
        }
    }
    Some(Placed {
        reference_start: a.reference_start as usize,
        reference_end: a.reference_end as usize,
        query_start: a.query_start as usize,
        query_end: a.query_end as usize,
        cost: a.cost,
        cigar,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_read_inside_a_stretch() {
        // The read is the stretch's letters 10..40 with one substitution and one deleted letter.
        let reference = b"TTTTTTTTTTACGTACGGATCCATGCAAGTTCGATCGGACTAGGGGGGGGGG".to_vec();
        let mut query = reference[10..40].to_vec();
        query[5] = b'T';
        query.remove(20);
        let p = place(&reference, &query, (i64::MIN, i64::MAX)).unwrap();
        assert_eq!((p.reference_start, p.reference_end, p.query_start, p.query_end), (10, 40, 0, 29));
        let (mut r, mut q) = (0u32, 0u32);
        for &(op, n) in &p.cigar {
            match op {
                b'=' | b'X' => {
                    r += n;
                    q += n;
                }
                b'D' => r += n,
                b'I' => q += n,
                _ => panic!("op {}", op as char),
            }
        }
        assert_eq!((r, q), (30, 29));
        assert!(p.cigar.contains(&(b'X', 1)) && p.cigar.contains(&(b'D', 1)));
        // A band that the read's diagonal (10) lies outside: no alignment.
        assert_eq!(place(&reference, &query, (-5, 5)).map(|p| p.reference_start).filter(|&s| s == 10), None);
    }
}
