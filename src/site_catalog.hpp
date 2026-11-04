#ifndef PGPHASE_SITE_CATALOG_HPP
#define PGPHASE_SITE_CATALOG_HPP

// dromos's site catalog (`dromos catalog`; the DRMSITES format, version 1, of
// dromos's docs/catalog-format.md): one file per graph and reference sample,
// each contig's records in position blocks of zstd-framed columns, an index
// frame per contig, and a header naming the graph it was built from.
//
// SiteCatalogFile is the file opened once per process: its preamble checked,
// its header read, each contig's index read the first time it is queried, and
// shared by every thread. SiteCatalogReader is one thread's: its zstd context
// and block buffers. A region query binary-searches the contig's index and
// decodes the blocks it selects into the GraphSites pgphase's old PGSITES2
// reader made of the same records' VCF lines.

#include "graph_sites.hpp"

#include <array>
#include <cstdint>
#include <memory>
#include <mutex>
#include <string>
#include <utility>
#include <vector>

struct ZSTD_DCtx_s;

namespace pgphase_collect {

/// One contig of a catalog, as its header gives it.
struct SiteCatalogContig {
    std::string name;   ///< chr20
    std::string chrom;  ///< CHM13#0#chr20: the records' CHROM
    uint64_t length = 0;
    std::array<unsigned char, 16> md5{};
    uint64_t blocks = 0;
    uint64_t records = 0;
    uint64_t index_offset = 0;
    uint64_t index_size = 0;
    /// 1-based, inclusive; empty: the whole contig.
    std::vector<std::pair<uint64_t, uint64_t>> regions;

    std::string md5_hex() const;
};

/// A contig's index: per block its position bounds, where it lies, what it
/// holds, its CRC, and each column's compressed and raw size.
struct SiteCatalogIndex {
    static constexpr size_t kColumns = 13;
    std::vector<uint64_t> first_pos, last_pos, max_end, reach, offset, records, min_allele, max_walk, crc32;
    std::vector<uint64_t> comp;  ///< block b's column c at b * kColumns + c
    std::vector<uint64_t> raw;
    std::vector<uint64_t> bytes;  ///< each block's stored bytes (its columns' compressed sizes)

    size_t blocks() const { return first_pos.size(); }
    /// Blocks [first, last) that may hold records overlapping [beg, end]
    /// (1-based, inclusive): two binary searches, on reach and first_pos.
    std::pair<size_t, size_t> overlapping(uint64_t beg, uint64_t end) const;
};

class SiteCatalogFile {
public:
    /// The catalog at `path`, opened once per process and kept: its preamble
    /// and header are checked and read on the first open. Throws when it is
    /// not a version 1 catalog or is damaged.
    static std::shared_ptr<const SiteCatalogFile> open(const std::string& path);

    explicit SiteCatalogFile(const std::string& path);
    ~SiteCatalogFile();
    SiteCatalogFile(const SiteCatalogFile&) = delete;
    SiteCatalogFile& operator=(const SiteCatalogFile&) = delete;

    const std::string& path() const { return path_; }
    const std::vector<std::pair<std::string, std::string>>& properties() const { return props_; }
    /// The property's value; empty when the header has none.
    std::string property(const std::string& key) const;
    const std::vector<SiteCatalogContig>& contigs() const { return contigs_; }
    /// The contig `name` names (chr20, or CHM13#0#chr20 of the reference
    /// sample alone); -1 when the catalog has none.
    int find_contig(const std::string& name) const;
    /// Contig k's index, read and checked the first time it is asked for.
    const SiteCatalogIndex& index(size_t k) const;
    /// Reads `size` bytes at `offset` into `out` (thread-safe: pread).
    void read_at(uint64_t offset, uint64_t size, std::string& out) const;
    uint64_t blocks_end() const { return blocks_end_; }

private:
    std::string path_;
    int fd_ = -1;
    uint64_t size_ = 0;
    uint64_t blocks_end_ = 0;
    std::vector<std::pair<std::string, std::string>> props_;
    std::vector<SiteCatalogContig> contigs_;
    mutable std::vector<std::unique_ptr<SiteCatalogIndex>> index_;
    mutable std::unique_ptr<std::once_flag[]> index_once_;
};

/// Whether `path` is a file starting with dromos's catalog magic (DRMSITES).
bool is_site_catalog_file(const std::string& path);

/// Throws, saying what to use instead, when `path` is pgcatalog's catalog: a
/// folder of <contig>.pgsites files, or one of them.
void refuse_pgcatalog_sites(const std::string& path);

/// One thread's reader of a shared catalog: its zstd context and buffers.
class SiteCatalogReader {
public:
    explicit SiteCatalogReader(std::shared_ptr<const SiteCatalogFile> file);
    ~SiteCatalogReader();
    SiteCatalogReader(const SiteCatalogReader&) = delete;
    SiteCatalogReader& operator=(const SiteCatalogReader&) = delete;

    const SiteCatalogFile& file() const { return *file_; }

    /// Appends the sites of `contig` overlapping [beg, end] (1-based,
    /// inclusive; a record spans POS to POS + max(len(REF), 1) - 1) to
    /// `out`, unsorted and not finalized, as the PGSITES2 reader appended the
    /// same records' VCF lines: a block whose every allele spans more than
    /// `max_span` bases is skipped unread, as is a record that no read can
    /// be matched to, each counted in out.stats.
    void append_sites(const std::string& contig, hts_pos_t beg, hts_pos_t end, hts_pos_t max_span,
                      GraphSiteCatalog& out);

private:
    struct Record {
        hts_pos_t pos = 0;
        uint8_t flags = 0;
        char pad = 0;
        uint64_t first = 0, last = 0;  ///< visits: id << 1 | reverse
        size_t walk0 = 0;              ///< its first walk in walks_
        size_t nwalk = 0;
        uint64_t level = 0;
        uint64_t ps_start = 0, ps_end = 0;
    };

    /// Reads contig k's block b (its CRC checked) and decodes its columns
    /// but the letters into records_, walks_, visits_ and the node table.
    void read_block(size_t k, size_t b);
    /// Column c of the block read, decompressed when stored compressed.
    void column(size_t c);
    /// The block's node letters (SEQ, EXC), decoded once.
    void letters();
    std::string block_damage(const std::string& why) const;
    /// Walk w's inner visits are visits_[walks_[r.walk0 + w] .. walks_[r.walk0 + w + 1]).
    size_t allele_length(const Record& r, size_t walk) const;
    void allele(const Record& r, size_t walk, std::string& out) const;
    GraphSite site(const SiteCatalogContig& contig, const Record& r) const;

    std::shared_ptr<const SiteCatalogFile> file_;
    ZSTD_DCtx_s* dctx_ = nullptr;
    std::string block_;
    std::array<std::string, SiteCatalogIndex::kColumns> cols_;
    std::array<std::string_view, SiteCatalogIndex::kColumns> col_;
    std::array<size_t, SiteCatalogIndex::kColumns> col_at_{};
    size_t block_k_ = 0, block_b_ = 0;
    bool letters_ready_ = false;
    std::vector<Record> records_;
    std::vector<size_t> walks_;
    std::vector<uint64_t> visits_;  ///< node-table index << 1 | reverse
    std::vector<uint64_t> node_ids_;
    std::vector<size_t> node_starts_;
    std::string letters_;
};

} // namespace pgphase_collect

#endif
