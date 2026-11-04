#include "site_catalog.hpp"

#include <algorithm>
#include <cerrno>
#include <cstring>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <unordered_set>

#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

#include <zlib.h>
#include <zstd.h>

namespace pgphase_collect {

namespace {

constexpr char kMagic[8] = {'D', 'R', 'M', 'S', 'I', 'T', 'E', 'S'};
constexpr uint32_t kVersion = 1;
constexpr uint64_t kPreamble = 64;
constexpr size_t K = SiteCatalogIndex::kColumns;
constexpr size_t kIndexArrays = 9;

// Columns (docs/catalog-format.md in dromos).
enum Column { POS, FLAGS, PAD, BOUND, NWALK, WLEN, WALK, LV, PS, NODEID, NODELEN, SEQ, EXC };
constexpr uint8_t kFlagReversed = 1, kFlagPad = 2, kFlagPs = 4;

uint64_t u64_le(const unsigned char* p) {
    uint64_t v = 0;
    for (int i = 7; i >= 0; --i) v = (v << 8) | p[i];
    return v;
}

uint64_t fnv1a(const unsigned char* p, size_t n) {
    uint64_t h = 0xCBF29CE484222325ULL;
    for (size_t i = 0; i < n; ++i) {
        h ^= p[i];
        h *= 0x100000001B3ULL;
    }
    return h;
}

int64_t unzigzag(uint64_t z) {
    return static_cast<int64_t>(z >> 1) ^ -static_cast<int64_t>(z & 1);
}

[[noreturn]] void damaged(const std::string& what) {
    throw std::runtime_error(what);
}

// One zstd frame of exactly `comp`, decompressed whole (its content size
// recorded; zstd checks its checksum when it has one).
std::string decompress_frame(ZSTD_DCtx* dctx, const std::string& comp, const std::string& what) {
    const unsigned long long raw = ZSTD_getFrameContentSize(comp.data(), comp.size());
    if (raw == ZSTD_CONTENTSIZE_ERROR || raw == ZSTD_CONTENTSIZE_UNKNOWN)
        damaged(what + ": damaged (not a zstd frame with its size)");
    if (ZSTD_findFrameCompressedSize(comp.data(), comp.size()) != comp.size())
        damaged(what + ": damaged (its frame's size)");
    std::string out(static_cast<size_t>(raw), '\0');
    const size_t n = ZSTD_decompressDCtx(dctx, out.data(), out.size(), comp.data(), comp.size());
    if (ZSTD_isError(n) || n != out.size()) damaged(what + ": damaged (" + (ZSTD_isError(n) ? ZSTD_getErrorName(n) : "short") + ")");
    return out;
}

ZSTD_DCtx* new_dctx() {
    ZSTD_DCtx* d = ZSTD_createDCtx();
    if (d == nullptr) throw std::runtime_error("cannot allocate a zstd context");
    // Long columns are compressed with windows up to 2^31.
    ZSTD_DCtx_setParameter(d, ZSTD_d_windowLogMax, 31);
    return d;
}

struct DCtxOwner {
    ZSTD_DCtx* d = new_dctx();
    ~DCtxOwner() { ZSTD_freeDCtx(d); }
};

// The header's fields, read in turn.
struct ByteReader {
    const std::string& data;
    size_t at = 0;
    const std::string& what;

    uint64_t u64() {
        if (data.size() - at < 8) damaged(what + ": damaged (it ends early)");
        const uint64_t v = u64_le(reinterpret_cast<const unsigned char*>(data.data()) + at);
        at += 8;
        return v;
    }
    std::string bytes(uint64_t n) {
        if (data.size() - at < n) damaged(what + ": damaged (it ends early)");
        std::string s = data.substr(at, static_cast<size_t>(n));
        at += static_cast<size_t>(n);
        return s;
    }
    std::string string() { return bytes(u64()); }
};

// A column read in record order.
struct Cursor {
    const unsigned char* p = nullptr;
    const unsigned char* e = nullptr;

    explicit Cursor(std::string_view s)
        : p(reinterpret_cast<const unsigned char*>(s.data())), e(p + s.size()) {}
    bool done() const { return p == e; }
    uint64_t varint() {
        uint64_t v = 0;
        for (int shift = 0; shift < 64; shift += 7) {
            if (p == e) damaged("a column ends early");
            const unsigned char b = *p++;
            v |= static_cast<uint64_t>(b & 0x7F) << shift;
            if ((b & 0x80) == 0) return v;
        }
        damaged("a varint over 10 bytes");
    }
    unsigned char byte() {
        if (p == e) damaged("a column ends early");
        return *p++;
    }
};

std::string hex(const unsigned char* d, size_t n) {
    static const char* digits = "0123456789abcdef";
    std::string s;
    for (size_t i = 0; i < n; ++i) {
        s.push_back(digits[d[i] >> 4]);
        s.push_back(digits[d[i] & 15]);
    }
    return s;
}

void put_visit(std::string& out, uint64_t v) {
    out.push_back((v & 1) ? '<' : '>');
    out += std::to_string(v >> 1);
}

char complement(char c) {
    switch (c) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        case 'a': return 'T';
        case 'c': return 'G';
        case 'g': return 'C';
        case 't': return 'A';
        default: return 'N';
    }
}

// Says something about a contig once per run (handles are per thread).
void warn_once(const std::string& key, const std::string& message) {
    static std::mutex mu;
    static std::unordered_set<std::string> said;
    std::lock_guard<std::mutex> lock(mu);
    if (said.insert(key).second) std::cerr << "warning: " << message << "\n";
}

}  // namespace

std::string SiteCatalogContig::md5_hex() const { return hex(md5.data(), md5.size()); }

std::pair<size_t, size_t> SiteCatalogIndex::overlapping(uint64_t beg, uint64_t end) const {
    const size_t first = static_cast<size_t>(std::lower_bound(reach.begin(), reach.end(), beg) - reach.begin());
    const size_t last = static_cast<size_t>(std::upper_bound(first_pos.begin() + static_cast<std::ptrdiff_t>(first),
                                                             first_pos.end(), end) - first_pos.begin());
    return {first, std::max(first, last)};
}

bool is_site_catalog_file(const std::string& path) {
    struct stat st;
    if (stat(path.c_str(), &st) != 0 || !S_ISREG(st.st_mode)) return false;
    const int fd = ::open(path.c_str(), O_RDONLY);
    if (fd < 0) return false;
    char magic[8] = {};
    const bool ok = ::pread(fd, magic, 8, 0) == 8 && std::memcmp(magic, kMagic, 8) == 0;
    ::close(fd);
    return ok;
}

void refuse_pgcatalog_sites(const std::string& path) {
    struct stat st;
    if (stat(path.c_str(), &st) != 0) return;
    bool pgsites = S_ISDIR(st.st_mode);
    if (S_ISREG(st.st_mode)) {
        const int fd = ::open(path.c_str(), O_RDONLY);
        char magic[7] = {};
        pgsites = fd >= 0 && ::pread(fd, magic, 7, 0) == 7 && std::memcmp(magic, "PGSITES", 7) == 0;
        if (fd >= 0) ::close(fd);
    }
    if (pgsites) SiteCatalogFile::open(path);  // throws, saying what the file is
}

std::shared_ptr<const SiteCatalogFile> SiteCatalogFile::open(const std::string& path) {
    // Kept for the process: its header is read once, however often it is opened.
    static std::mutex mu;
    static std::map<std::string, std::shared_ptr<const SiteCatalogFile>> opened;
    std::lock_guard<std::mutex> lock(mu);
    std::shared_ptr<const SiteCatalogFile>& slot = opened[path];
    if (!slot) slot = std::make_shared<const SiteCatalogFile>(path);
    return slot;
}

SiteCatalogFile::SiteCatalogFile(const std::string& path) : path_(path) {
    fd_ = ::open(path.c_str(), O_RDONLY);
    if (fd_ < 0) throw std::runtime_error("cannot open site catalog " + path + ": " + std::strerror(errno));
    try {
        struct stat st;
        if (fstat(fd_, &st) != 0) throw std::runtime_error("cannot stat site catalog " + path);
        if (S_ISDIR(st.st_mode))
            throw std::runtime_error(path + ": a folder; pgphase reads dromos's one-file site catalog, not pgcatalog's "
                                     "folders of <contig>.pgsites: make it with dromos catalog");
        size_ = static_cast<uint64_t>(st.st_size);
        // The preamble, checked in the specification's order.
        unsigned char pre[kPreamble] = {};
        const ssize_t got = ::pread(fd_, pre, kPreamble, 0);
        const size_t have = got < 0 ? 0 : static_cast<size_t>(got);
        if (have < 8 || std::memcmp(pre, kMagic, 8) != 0) {
            if (have >= 7 && std::memcmp(pre, "PGSITES", 7) == 0)
                throw std::runtime_error(path + ": a pgcatalog .pgsites file (one contig's catalog), which pgphase "
                                         "no longer reads: make the graph's catalog with dromos catalog");
            throw std::runtime_error(path + ": not a dromos site catalog (no DRMSITES magic)");
        }
        if (have < kPreamble)
            throw std::runtime_error(path + ": truncated (a catalog's first 64 bytes are its preamble; it has " +
                                     std::to_string(have) + ")");
        uint32_t version = 0;
        for (int i = 3; i >= 0; --i) version = (version << 8) | pre[8 + i];
        if (version != kVersion)
            throw std::runtime_error(path + ": a site catalog of format version " + std::to_string(version) +
                                     "; pgphase reads version " + std::to_string(kVersion) +
                                     ": rebuild it with a dromos that writes version 1");
        if (u64_le(pre + 56) != fnv1a(pre, 56)) damaged(path + ": damaged (its preamble's hash differs)");
        const uint64_t recorded = u64_le(pre + 16);
        if (recorded != size_) {
            if (size_ < recorded)
                damaged(path + ": truncated (" + std::to_string(size_) + " of its " + std::to_string(recorded) +
                        " bytes)");
            damaged(path + ": " + std::to_string(size_) + " bytes where its preamble says " +
                    std::to_string(recorded) + ": damaged");
        }
        const uint64_t header_at = u64_le(pre + 24), header_size = u64_le(pre + 32);
        const uint64_t index_at = u64_le(pre + 40), index_size = u64_le(pre + 48);
        if (header_at < kPreamble || header_size == 0 || header_size > size_ || header_at > size_ - header_size ||
            index_at < kPreamble || index_size > header_at || index_at > header_at - index_size)
            damaged(path + ": damaged (its preamble points outside the file)");
        blocks_end_ = index_at;

        DCtxOwner z;
        std::string comp;
        read_at(header_at, header_size, comp);
        const std::string what = path + ": its header";
        const std::string header = decompress_frame(z.d, comp, what);
        ByteReader h{header, 0, what};
        if (h.u64() != kVersion) damaged(path + ": damaged (its header's version differs from its preamble's)");
        const uint64_t nprops = h.u64();
        for (uint64_t i = 0; i < nprops; ++i) {
            std::string key = h.string();
            props_.emplace_back(std::move(key), h.string());
        }
        const uint64_t ncontigs = h.u64();
        for (uint64_t i = 0; i < ncontigs; ++i) {
            SiteCatalogContig c;
            c.name = h.string();
            c.chrom = h.string();
            c.length = h.u64();
            const std::string md5 = h.bytes(16);
            std::memcpy(c.md5.data(), md5.data(), 16);
            c.blocks = h.u64();
            c.records = h.u64();
            c.index_offset = h.u64();
            c.index_size = h.u64();
            if (c.index_offset < index_at || c.index_size == 0 || c.index_size > index_size ||
                c.index_offset - index_at > index_size - c.index_size)
                damaged(path + ": damaged (contig " + c.name + "'s index lies outside the index)");
            const uint64_t nregions = h.u64();
            for (uint64_t r = 0; r < nregions; ++r) {
                const uint64_t a = h.u64();
                c.regions.emplace_back(a, h.u64());
            }
            contigs_.push_back(std::move(c));
        }
        if (h.at != header.size()) damaged(path + ": damaged (its header has bytes left over)");
    } catch (...) {
        ::close(fd_);
        fd_ = -1;
        throw;
    }
    index_.resize(contigs_.size());
    index_once_ = std::make_unique<std::once_flag[]>(contigs_.size());
}

SiteCatalogFile::~SiteCatalogFile() {
    if (fd_ >= 0) ::close(fd_);
}

void SiteCatalogFile::read_at(uint64_t offset, uint64_t size, std::string& out) const {
    out.resize(static_cast<size_t>(size));
    size_t done = 0;
    while (done < size) {
        const ssize_t n = ::pread(fd_, out.data() + done, static_cast<size_t>(size) - done,
                                  static_cast<off_t>(offset + done));
        if (n < 0 && errno == EINTR) continue;
        if (n <= 0) throw std::runtime_error(path_ + ": cannot read " + std::to_string(size) + " bytes at " +
                                             std::to_string(offset) + (n < 0 ? std::string(": ") + std::strerror(errno) : ": truncated"));
        done += static_cast<size_t>(n);
    }
}

std::string SiteCatalogFile::property(const std::string& key) const {
    for (const auto& kv : props_)
        if (kv.first == key) return kv.second;
    return {};
}

int SiteCatalogFile::find_contig(const std::string& name) const {
    for (size_t k = 0; k < contigs_.size(); ++k)
        if (graph_chrom_names_match(name, contigs_[k].chrom)) return static_cast<int>(k);
    return -1;
}

const SiteCatalogIndex& SiteCatalogFile::index(size_t k) const {
    std::call_once(index_once_[k], [&]() {
        const SiteCatalogContig& c = contigs_[k];
        const std::string what = path_ + ": " + c.name + "'s index";
        DCtxOwner z;
        std::string comp;
        read_at(c.index_offset, c.index_size, comp);
        const std::string data = decompress_frame(z.d, comp, what);
        const uint64_t n = c.blocks;
        if (n > data.size() / 8 || data.size() != 16 + 8 * n * (kIndexArrays + 2 * K) ||
            u64_le(reinterpret_cast<const unsigned char*>(data.data())) != n)
            damaged(what + ": damaged (the wrong size)");
        const auto* p = reinterpret_cast<const unsigned char*>(data.data());
        if (u64_le(p + 8) != K) damaged(what + ": damaged (the wrong number of columns)");
        auto x = std::make_unique<SiteCatalogIndex>();
        std::vector<uint64_t>* arrays[kIndexArrays] = {&x->first_pos, &x->last_pos, &x->max_end,
                                                       &x->reach,     &x->offset,   &x->records,
                                                       &x->min_allele, &x->max_walk, &x->crc32};
        size_t at = 16;
        for (std::vector<uint64_t>* a : arrays) {
            a->resize(n);
            for (uint64_t b = 0; b < n; ++b, at += 8) (*a)[b] = u64_le(p + at);
        }
        x->comp.resize(n * K);
        x->raw.resize(n * K);
        for (uint64_t i = 0; i < n * K; ++i, at += 8) x->comp[i] = u64_le(p + at);
        for (uint64_t i = 0; i < n * K; ++i, at += 8) x->raw[i] = u64_le(p + at);
        x->bytes.assign(n, 0);
        uint64_t records = 0;
        for (uint64_t b = 0; b < n; ++b) {
            uint64_t total = 0;
            for (size_t col = 0; col < K; ++col) {
                const uint64_t cs = x->comp[b * K + col], rs = x->raw[b * K + col];
                if ((cs == 0) != (rs == 0) || cs > rs || cs > blocks_end_)
                    damaged(what + ": damaged (block " + std::to_string(b) + "'s column sizes)");
                total += cs;
            }
            x->bytes[b] = total;
            const uint64_t off = x->offset[b];
            if (off < kPreamble || off > blocks_end_ || total > blocks_end_ - off)
                damaged(what + ": damaged (block " + std::to_string(b) + " lies outside the blocks)");
            if (b > 0 && (x->first_pos[b] < x->first_pos[b - 1] || x->reach[b] < x->reach[b - 1]))
                damaged(what + ": damaged (its blocks are out of order)");
            if (x->records[b] == 0) damaged(what + ": damaged (block " + std::to_string(b) + " holds no records)");
            records += x->records[b];
        }
        if (records != c.records) damaged(what + ": damaged (its blocks' records differ from its header's)");
        index_[k] = std::move(x);
    });
    return *index_[k];
}

// --- SiteCatalogReader -------------------------------------------------------

SiteCatalogReader::SiteCatalogReader(std::shared_ptr<const SiteCatalogFile> file)
    : file_(std::move(file)), dctx_(new_dctx()) {}

SiteCatalogReader::~SiteCatalogReader() { ZSTD_freeDCtx(dctx_); }

void SiteCatalogReader::read_block(size_t k, size_t b) {
    const SiteCatalogIndex& ix = file_->index(k);
    block_k_ = k;
    block_b_ = b;
    file_->read_at(ix.offset[b], ix.bytes[b], block_);
    uLong crc = crc32(0L, Z_NULL, 0);
    for (size_t at = 0; at < block_.size();) {
        const size_t n = std::min<size_t>(block_.size() - at, 1u << 30);
        crc = crc32(crc, reinterpret_cast<const Bytef*>(block_.data() + at), static_cast<uInt>(n));
        at += n;
    }
    if (static_cast<uint64_t>(crc) != ix.crc32[b]) damaged(block_damage("its CRC differs"));
    // Every column but the letters (SEQ, EXC), which only a block that keeps
    // a record needs (letters()).
    size_t at = 0;
    for (size_t c = 0; c < K; ++c) {
        const uint64_t cs = ix.comp[b * K + c];
        col_at_[c] = at;
        at += static_cast<size_t>(cs);
        if (c != SEQ && c != EXC) column(c);
    }
    letters_ready_ = false;

    // The node table.
    node_ids_.clear();
    node_starts_.assign(1, 0);
    {
        Cursor ids(col_[NODEID]), lens(col_[NODELEN]);
        uint64_t id = 0;
        size_t total = 0;
        while (!ids.done()) {
            const uint64_t d = ids.varint();
            if (!node_ids_.empty() && d == 0) damaged(block_damage("node ids do not increase"));
            id += d;
            node_ids_.push_back(id);
            total += static_cast<size_t>(lens.varint());
            node_starts_.push_back(total);
        }
        if (!lens.done()) damaged(block_damage("node lengths left over"));
    }

    // The records.
    records_.clear();
    walks_.clear();
    visits_.clear();
    Cursor pos(col_[POS]), flags(col_[FLAGS]), pad(col_[PAD]), bound(col_[BOUND]), nwalk(col_[NWALK]),
        wlen(col_[WLEN]), walk(col_[WALK]), lv(col_[LV]), ps(col_[PS]);
    const int64_t nodes = static_cast<int64_t>(node_ids_.size());
    hts_pos_t at_pos = static_cast<hts_pos_t>(ix.first_pos[b]);
    uint64_t prev_first = 0;
    int64_t widx = 0;
    for (uint64_t i = 0; i < ix.records[b]; ++i) {
        Record r;
        at_pos += static_cast<hts_pos_t>(pos.varint());
        r.pos = at_pos;
        r.flags = flags.byte();
        if (r.flags & ~(kFlagReversed | kFlagPad | kFlagPs)) damaged(block_damage("unknown record flags"));
        r.pad = (r.flags & kFlagPad) ? static_cast<char>(pad.byte()) : 0;
        r.first = prev_first + static_cast<uint64_t>(unzigzag(bound.varint()));
        r.last = r.first + static_cast<uint64_t>(unzigzag(bound.varint()));
        prev_first = r.first;
        r.nwalk = static_cast<size_t>(nwalk.varint());
        if (r.nwalk < 1) damaged(block_damage("a record without walks"));
        r.walk0 = walks_.size();
        for (size_t w = 0; w < r.nwalk; ++w) {
            walks_.push_back(visits_.size());
            const uint64_t m = wlen.varint();
            for (uint64_t j = 0; j < m; ++j) {
                const uint64_t x = walk.varint();
                widx += unzigzag(x >> 1);
                if (widx < 0 || widx >= nodes) damaged(block_damage("a walk visits a node outside its table"));
                visits_.push_back(static_cast<uint64_t>(widx) << 1 | (x & 1));
            }
        }
        r.level = lv.varint();
        if (r.flags & kFlagPs) {
            const uint64_t base = (r.flags & kFlagReversed) ? r.last ^ 1 : r.first;
            r.ps_start = base + static_cast<uint64_t>(unzigzag(ps.varint()));
            r.ps_end = r.ps_start + static_cast<uint64_t>(unzigzag(ps.varint()));
        }
        records_.push_back(r);
    }
    walks_.push_back(visits_.size());
    if (!pos.done() || !flags.done() || !pad.done() || !bound.done() || !nwalk.done() || !wlen.done() ||
        !walk.done() || !lv.done() || !ps.done())
        damaged(block_damage("columns left over"));
}

std::string SiteCatalogReader::block_damage(const std::string& why) const {
    return file_->path() + ": damaged (" + file_->contigs()[block_k_].name + "'s block " + std::to_string(block_b_) +
           ": " + why + ")";
}

void SiteCatalogReader::column(size_t c) {
    const SiteCatalogIndex& ix = file_->index(block_k_);
    const uint64_t cs = ix.comp[block_b_ * K + c], rs = ix.raw[block_b_ * K + c];
    const std::string_view stored(block_.data() + col_at_[c], static_cast<size_t>(cs));
    if (cs == rs) {
        col_[c] = stored;  // empty, or stored as it is
        return;
    }
    std::string& dst = cols_[c];
    dst.resize(static_cast<size_t>(rs));
    const size_t n = ZSTD_decompressDCtx(dctx_, dst.data(), dst.size(), stored.data(), stored.size());
    if (ZSTD_isError(n) || n != dst.size()) damaged(block_damage("a column does not decompress"));
    col_[c] = dst;
}

void SiteCatalogReader::letters() {
    if (letters_ready_) return;
    column(SEQ);
    column(EXC);
    const size_t total = node_starts_.back();
    const std::string_view seq = col_[SEQ];
    if (seq.size() != (total + 3) / 4) damaged(block_damage("the letters' column has the wrong size"));
    static const char acgt[4] = {'A', 'C', 'G', 'T'};
    letters_.resize(total);
    for (size_t i = 0; i < total; ++i)
        letters_[i] = acgt[(static_cast<unsigned char>(seq[i >> 2]) >> ((i & 3) * 2)) & 3];
    Cursor exc(col_[EXC]);
    size_t e = 0;
    while (!exc.done()) {
        e += static_cast<size_t>(exc.varint());
        if (e >= total) damaged(block_damage("an exception past the letters"));
        letters_[e] = static_cast<char>(exc.byte());
    }
    letters_ready_ = true;
}

size_t SiteCatalogReader::allele_length(const Record& r, size_t w) const {
    size_t n = (r.flags & kFlagPad) ? 1 : 0;
    for (size_t j = walks_[r.walk0 + w]; j < walks_[r.walk0 + w + 1]; ++j) {
        const size_t node = static_cast<size_t>(visits_[j] >> 1);
        n += node_starts_[node + 1] - node_starts_[node];
    }
    return n;
}

void SiteCatalogReader::allele(const Record& r, size_t w, std::string& out) const {
    out.clear();
    out.reserve(allele_length(r, w));
    if (r.flags & kFlagPad) out.push_back(r.pad);
    for (size_t j = walks_[r.walk0 + w]; j < walks_[r.walk0 + w + 1]; ++j) {
        const size_t node = static_cast<size_t>(visits_[j] >> 1);
        const char* a = letters_.data() + node_starts_[node];
        const size_t n = node_starts_[node + 1] - node_starts_[node];
        if (visits_[j] & 1) {
            for (size_t i = n; i > 0; --i) out.push_back(complement(a[i - 1]));
        } else {
            out.append(a, n);
        }
    }
}

// The GraphSite the PGSITES2 reader parsed from the record's VCF line
// (append_graph_site_from_vcf_data_line without the allele traversal text).
GraphSite SiteCatalogReader::site(const SiteCatalogContig& contig, const Record& r) const {
    GraphSite s;
    s.chrom = contig.chrom;
    const size_t hash = contig.chrom.rfind('#');
    s.ref_contig = hash == std::string::npos ? contig.chrom : contig.chrom.substr(hash + 1);
    s.pos = r.pos;
    const bool rev = (r.flags & kFlagReversed) != 0;
    put_visit(s.id, rev ? r.last ^ 1 : r.first);
    put_visit(s.id, rev ? r.first ^ 1 : r.last);
    allele(r, 0, s.ref);
    // The ALT field is the other alleles joined by ','; the VCF reader splits
    // it unless it is empty or '.', so a lone empty ALT gives no alleles.
    if (r.nwalk == 2) {
        std::string a;
        allele(r, 1, a);
        if (!a.empty() && a != ".") s.alts.push_back(std::move(a));
    } else {
        for (size_t w = 1; w < r.nwalk; ++w) {
            s.alts.emplace_back();
            allele(r, w, s.alts.back());
        }
    }
    return s;
}

void SiteCatalogReader::append_sites(const std::string& contig_name, hts_pos_t beg, hts_pos_t end,
                                     hts_pos_t max_span, GraphSiteCatalog& out) {
    const int k = file_->find_contig(contig_name);
    if (k < 0) {
        warn_once("site catalog contig " + contig_name,
                  "no sites on " + contig_name + ": the site catalog " + file_->path() + " has no such contig");
        return;
    }
    const SiteCatalogContig& contig = file_->contigs()[static_cast<size_t>(k)];
    if (!contig.regions.empty()) {
        bool inside = false;
        for (const auto& [a, b] : contig.regions)
            inside = inside || (static_cast<hts_pos_t>(a) <= beg && end <= static_cast<hts_pos_t>(b));
        if (!inside)
            warn_once("site catalog regions " + contig_name,
                      "the site catalog holds " + contig.name +
                          " in regions alone: outside them it has no sites (a query of " + contig_name + ":" +
                          std::to_string(beg) + "-" + std::to_string(end) + ")");
    }
    if (contig.blocks == 0 || end < 1 || beg > end) return;
    const SiteCatalogIndex& ix = file_->index(static_cast<size_t>(k));
    const uint64_t qbeg = static_cast<uint64_t>(std::max<hts_pos_t>(beg, 0));
    const auto [first, last] = ix.overlapping(qbeg, static_cast<uint64_t>(end));
    GraphSiteLoadStats& stats = out.stats;
    for (size_t b = first; b < last; ++b) {
        if (ix.max_end[b] < qbeg) continue;
        if (static_cast<hts_pos_t>(ix.min_allele[b]) > max_span) {
            stats.beyond_read_span += ix.records[b];
            continue;
        }
        if (ix.records[b] == 1 && ix.max_walk[b] > kMaxSiteWalkSteps) {
            ++stats.walk_too_long;
            continue;
        }
        read_block(static_cast<size_t>(k), b);
        for (const Record& r : records_) {
            const hts_pos_t ref_len = static_cast<hts_pos_t>(allele_length(r, 0));
            if (r.pos > end || r.pos + std::max<hts_pos_t>(1, ref_len) - 1 < beg) continue;
            // unmatchable_record: the shortest allele (REF, then each ALT; an
            // empty ALT field one of no bases), less the padding base when
            // the alleles' lengths differ, longer than any alignment spans...
            hts_pos_t shortest = ref_len;
            bool equal = true;
            if (r.nwalk == 1) {
                equal = ref_len == 0;
                shortest = 0;
            }
            for (size_t w = 1; w < r.nwalk; ++w) {
                const hts_pos_t l = static_cast<hts_pos_t>(allele_length(r, w));
                equal = equal && l == ref_len;
                shortest = std::min(shortest, l);
            }
            if ((equal ? shortest : shortest - 1) > max_span) {
                ++stats.beyond_read_span;
                continue;
            }
            // ...or a walk (its boundary visits counted) longer than the matcher takes.
            size_t longest = 0;
            for (size_t w = 0; w < r.nwalk; ++w)
                longest = std::max(longest, walks_[r.walk0 + w + 1] - walks_[r.walk0 + w] + 2);
            if (longest > kMaxSiteWalkSteps) {
                ++stats.walk_too_long;
                continue;
            }
            letters();
            ++stats.data_lines;
            if (r.pos <= 0) {
                ++stats.bad_position;
                continue;
            }
            GraphSite s = site(contig, r);
            if (!graph_allele_is_representable(s.ref)) {
                ++stats.unsupported_allele;
                continue;
            }
            bool representable = true;
            for (const std::string& a : s.alts) representable = representable && graph_allele_is_representable(a);
            if (!representable) {
                ++stats.unsupported_allele;
                continue;
            }
            const Record& rr = r;
            s.allele_walks.resize(rr.nwalk);
            for (size_t w = 0; w < rr.nwalk; ++w) {
                GraphWalk& gw = s.allele_walks[w];
                gw.reserve(walks_[rr.walk0 + w + 1] - walks_[rr.walk0 + w] + 2);
                gw.push_back(GraphWalkStep{rr.first >> 1, (rr.first & 1) != 0});
                for (size_t j = walks_[rr.walk0 + w]; j < walks_[rr.walk0 + w + 1]; ++j)
                    gw.push_back(GraphWalkStep{node_ids_[static_cast<size_t>(visits_[j] >> 1)], (visits_[j] & 1) != 0});
                gw.push_back(GraphWalkStep{rr.last >> 1, (rr.last & 1) != 0});
            }
            s.level = static_cast<int>(rr.level);
            if (rr.flags & kFlagPs) {
                put_visit(s.parent, rr.ps_start);
                put_visit(s.parent, rr.ps_end);
            }
            s.ref_beg = s.pos;
            s.ref_end = s.pos;
            s.has_spanning_deletion = std::find(s.alts.begin(), s.alts.end(), "*") != s.alts.end();
            s.skip_reason = graph_site_validation_skip_reason(s);
            s.eligible = s.skip_reason.empty();
            ++stats.sites_parsed;
            out.sites.push_back(std::move(s));
        }
    }
}

} // namespace pgphase_collect
