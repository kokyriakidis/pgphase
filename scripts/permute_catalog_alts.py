#!/usr/bin/env python3
"""A dromos site catalog (DRMSITES, format version 1) with each record's ALT alleles in another order.

The catalog's records, positions, walks, node tables and header are kept; only the order of each record's walks after
the reference's (its ALTs, and their AT entries with them) changes: reversed (--order reverse, the default) or rotated
one place (--order rotate: the first ALT last). For the gates (scripts/validate_collect_gates.sh): pgphase's output
must not depend on the order a catalog lists a site's ALTs in.

The format is dromos's docs/catalog-format.md; zstd comes from the libzstd pgphase links (ctypes), so nothing beyond
the Python standard library is needed.
"""

import argparse
import ctypes
import ctypes.util
import struct
import sys
import zlib

MAGIC = b"DRMSITES"
K = 13
COL_NWALK, COL_WLEN, COL_WALK = 4, 5, 6


def load_zstd():
    for name in (ctypes.util.find_library("zstd"), "libzstd.so.1", "libzstd.dylib"):
        if not name:
            continue
        try:
            lib = ctypes.CDLL(name)
        except OSError:
            continue
        lib.ZSTD_getFrameContentSize.restype = ctypes.c_ulonglong
        lib.ZSTD_getFrameContentSize.argtypes = [ctypes.c_char_p, ctypes.c_size_t]
        lib.ZSTD_decompress.restype = ctypes.c_size_t
        lib.ZSTD_decompress.argtypes = [ctypes.c_char_p, ctypes.c_size_t, ctypes.c_char_p, ctypes.c_size_t]
        lib.ZSTD_isError.restype = ctypes.c_uint
        lib.ZSTD_isError.argtypes = [ctypes.c_size_t]
        lib.ZSTD_compressBound.restype = ctypes.c_size_t
        lib.ZSTD_compressBound.argtypes = [ctypes.c_size_t]
        lib.ZSTD_createCCtx.restype = ctypes.c_void_p
        lib.ZSTD_CCtx_setParameter.restype = ctypes.c_size_t
        lib.ZSTD_CCtx_setParameter.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int]
        lib.ZSTD_compress2.restype = ctypes.c_size_t
        lib.ZSTD_compress2.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_size_t, ctypes.c_char_p,
                                       ctypes.c_size_t]
        return lib
    sys.exit("error: libzstd not found")


ZSTD = load_zstd()
CCTX = ZSTD.ZSTD_createCCtx()
ZSTD.ZSTD_CCtx_setParameter(CCTX, 100, 9)   # ZSTD_c_compressionLevel
ZSTD.ZSTD_CCtx_setParameter(CCTX, 200, 1)   # ZSTD_c_contentSizeFlag
ZSTD.ZSTD_CCtx_setParameter(CCTX, 201, 1)   # ZSTD_c_checksumFlag


def decompress(frame, size=None):
    if size is None:
        size = ZSTD.ZSTD_getFrameContentSize(frame, len(frame))
        if size >= (1 << 63):
            sys.exit("error: a zstd frame without its content size")
    out = ctypes.create_string_buffer(max(size, 1))
    n = ZSTD.ZSTD_decompress(out, size, frame, len(frame))
    if ZSTD.ZSTD_isError(n) or n != size:
        sys.exit("error: a damaged zstd frame")
    return out.raw[:size]


def compress(data):
    cap = ZSTD.ZSTD_compressBound(len(data))
    out = ctypes.create_string_buffer(cap)
    n = ZSTD.ZSTD_compress2(CCTX, out, cap, data, len(data))
    if ZSTD.ZSTD_isError(n):
        sys.exit("error: zstd compression failed")
    return out.raw[:n]


def fnv1a(data):
    h = 0xCBF29CE484222325
    for b in data:
        h = ((h ^ b) * 0x100000001B3) & 0xFFFFFFFFFFFFFFFF
    return h


class Reader:
    def __init__(self, data):
        self.data, self.at = data, 0

    def u64(self):
        v = struct.unpack_from("<Q", self.data, self.at)[0]
        self.at += 8
        return v

    def raw(self, n):
        v = self.data[self.at:self.at + n]
        self.at += n
        return v

    def string(self):
        return self.raw(self.u64())

    def varint(self):
        v, shift = 0, 0
        while True:
            b = self.data[self.at]
            self.at += 1
            v |= (b & 0x7F) << shift
            shift += 7
            if b < 0x80:
                return v

    def done(self):
        return self.at == len(self.data)


def varint(v):
    out = bytearray()
    while True:
        b = v & 0x7F
        v >>= 7
        if v:
            out.append(b | 0x80)
        else:
            out.append(b)
            return bytes(out)


def zigzag(d):
    return (d << 1) ^ (d >> 63) if d < 0 else d << 1


def unzigzag(z):
    return (z >> 1) ^ -(z & 1)


def string(b):
    return struct.pack("<Q", len(b)) + b


def permute_block(cols, records, order):
    """Columns 5 and 6 of a block with each record's ALT walks reordered; the others as they are."""
    nwalk_r, wlen_r, walk_r = Reader(cols[COL_NWALK]), Reader(cols[COL_WLEN]), Reader(cols[COL_WALK])
    prev = 0
    new_wlen, new_walk, prev_out = bytearray(), bytearray(), 0
    for _ in range(records):
        walks = []
        for _ in range(nwalk_r.varint()):
            visits = []
            for _ in range(wlen_r.varint()):
                v = walk_r.varint()
                prev += unzigzag(v >> 1)
                visits.append((prev, v & 1))
            walks.append(visits)
        alts = walks[1:]
        alts = alts[::-1] if order == "reverse" else alts[1:] + alts[:1]
        for visits in [walks[0]] + alts:
            new_wlen += varint(len(visits))
            for i, rev in visits:
                new_walk += varint(zigzag(i - prev_out) << 1 | rev)
                prev_out = i
    if not (nwalk_r.done() and wlen_r.done() and walk_r.done()):
        sys.exit("error: a block's walk columns are not consumed exactly")
    out = list(cols)
    out[COL_WLEN], out[COL_WALK] = bytes(new_wlen), bytes(new_walk)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("catalog")
    ap.add_argument("out")
    ap.add_argument("--order", choices=("reverse", "rotate"), default="reverse")
    args = ap.parse_args()
    data = open(args.catalog, "rb").read()
    if data[:8] != MAGIC or len(data) < 64:
        sys.exit("error: not a dromos site catalog")
    version, _, size, header_off, header_size, index_off, index_size, hsh = struct.unpack_from("<IIQQQQQQ", data, 8)
    if version != 1 or hsh != fnv1a(data[:56]) or size != len(data):
        sys.exit("error: not a version 1 catalog, or damaged")

    header = Reader(decompress(data[header_off:header_off + header_size]))
    if header.u64() != version:
        sys.exit("error: the header's version differs from the preamble's")
    props = [(header.string(), header.string()) for _ in range(header.u64())]
    contigs = []
    for _ in range(header.u64()):
        c = {"name": header.string(), "chrom": header.string(), "length": header.u64(), "md5": header.raw(16),
             "blocks": header.u64(), "records": header.u64(), "index_offset": header.u64(), "index_size": header.u64()}
        c["regions"] = [(header.u64(), header.u64()) for _ in range(header.u64())]
        contigs.append(c)
    if not header.done():
        sys.exit("error: the header is not consumed exactly")

    body = bytearray()
    indexes = []
    for c in contigs:
        idx = Reader(decompress(data[c["index_offset"]:c["index_offset"] + c["index_size"]]))
        b, k = idx.u64(), idx.u64()
        if k != K or b != c["blocks"]:
            sys.exit("error: an index with unexpected columns or blocks")
        arrays = [[idx.u64() for _ in range(b)] for _ in range(9)]
        comp = [idx.u64() for _ in range(b * K)]
        raw = [idx.u64() for _ in range(b * K)]
        first_pos, last_pos, max_end, reach, offset, records, min_allele, max_walk, crc = arrays
        new_offset, new_crc, new_comp, new_raw = [], [], [], []
        for bi in range(b):
            at = offset[bi]
            block = data[at:at + sum(comp[bi * K:(bi + 1) * K])]
            if zlib.crc32(block) != crc[bi]:
                sys.exit("error: a block's CRC does not match")
            cols, pos = [], 0
            for ci in range(K):
                cs, rs = comp[bi * K + ci], raw[bi * K + ci]
                stored = block[pos:pos + cs]
                pos += cs
                cols.append(stored if cs == rs else decompress(stored, rs))
            cols = permute_block(cols, records[bi], args.order)
            stored_block = bytearray()
            for ci in range(K):
                if ci in (COL_WLEN, COL_WALK):
                    col = cols[ci]
                    z = compress(col) if col else b""
                    stored = z if col and len(z) < len(col) else col
                    new_comp.append(len(stored))
                    new_raw.append(len(col))
                else:
                    cs = comp[bi * K + ci]
                    start = sum(comp[bi * K:bi * K + ci])
                    stored = block[start:start + cs]
                    new_comp.append(cs)
                    new_raw.append(raw[bi * K + ci])
                stored_block += stored
            new_offset.append(64 + len(body))
            new_crc.append(zlib.crc32(stored_block))
            body += stored_block
        arrays = [first_pos, last_pos, max_end, reach, new_offset, records, min_allele, max_walk, new_crc]
        frame = struct.pack("<QQ", b, K) + b"".join(struct.pack("<%dQ" % b, *a) for a in arrays)
        frame += struct.pack("<%dQ" % (b * K), *new_comp) + struct.pack("<%dQ" % (b * K), *new_raw)
        indexes.append(compress(frame))

    index_off = 64 + len(body)
    index_blob = bytearray()
    for c, frame in zip(contigs, indexes):
        c["index_offset"], c["index_size"] = index_off + len(index_blob), len(frame)
        index_blob += frame
    hdr = struct.pack("<QQ", version, len(props)) + b"".join(string(k) + string(v) for k, v in props)
    hdr += struct.pack("<Q", len(contigs))
    for c in contigs:
        hdr += string(c["name"]) + string(c["chrom"]) + struct.pack("<Q", c["length"]) + c["md5"]
        hdr += struct.pack("<QQQQQ", c["blocks"], c["records"], c["index_offset"], c["index_size"], len(c["regions"]))
        hdr += b"".join(struct.pack("<QQ", s, e) for s, e in c["regions"])
    header_frame = compress(hdr)
    header_off = index_off + len(index_blob)
    total = header_off + len(header_frame)
    pre = MAGIC + struct.pack("<IIQQQQQ", version, 0, total, header_off, len(header_frame), index_off, len(index_blob))
    pre += struct.pack("<Q", fnv1a(pre))
    with open(args.out, "wb") as f:
        f.write(pre + bytes(body) + bytes(index_blob) + header_frame)


if __name__ == "__main__":
    main()
