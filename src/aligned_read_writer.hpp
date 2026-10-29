#ifndef PGPHASE_ALIGNED_READ_WRITER_HPP
#define PGPHASE_ALIGNED_READ_WRITER_HPP

// Aligned read output for collect-graph-variation (--out-bam/--out-cram/--out-sam): the input records of the regions
// processed, each primary with its read's union label as HP/PS and, with --refine-aln, its refined alignment
// (refined_alignment.hpp).

#include "bam_io_types.hpp"
#include "refined_alignment.hpp"

#include <cstddef>
#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <htslib/faidx.h>
#include <htslib/sam.h>

namespace pgphase_collect {

/// A read's union label: haplotype 1 or 2 and its phase set, or haplotype 0 (unphased).
using ReadLabel = std::function<std::pair<int, hts_pos_t>(const std::string& qname)>;

/// Writes the input alignment (one BAM/CRAM) over the regions processed, sorted, then indexes it.
class AlignedReadWriter {
public:
    explicit AlignedReadWriter(const Options& opts);
    ~AlignedReadWriter();
    AlignedReadWriter(const AlignedReadWriter&) = delete;
    AlignedReadWriter& operator=(const AlignedReadWriter&) = delete;

    /// The regions (contig, 1-based inclusive span) in processing order must follow the input header's contig order
    /// and, within a contig, not overlap or go back: otherwise the output could not be sorted.
    void check_region_order(const std::vector<std::pair<std::string, std::pair<hts_pos_t, hts_pos_t>>>& regions) const;

    /// Every input record overlapping `spans` (1-based, inclusive, sorted, disjoint) on `contig` and not written for
    /// an earlier span: a primary with its refinement (`refined`, by read name) and label as HP/PS, any other record
    /// without HP/PS.
    void write_contig(const std::string& contig, const std::vector<std::pair<hts_pos_t, hts_pos_t>>& spans,
                      const std::unordered_map<std::string, const RefinedAlignment*>& refined, const ReadLabel& label);

    /// Writes what is held, closes and indexes the output.
    void finish();

    std::size_t records() const { return n_records_; }
    std::size_t refined() const { return n_refined_; }
    std::size_t left_aligned() const { return n_left_aligned_; }

private:
    struct Held {
        hts_pos_t pos;
        std::size_t order;
        std::unique_ptr<bam1_t, AlignmentDeleter> b;
    };
    void apply(bam1_t* b, const RefinedAlignment& r);
    bool left_align(bam1_t* b);
    void update_tags(bam1_t* b);
    const std::string& contig_sequence(int tid, const std::string& contig);
    void write_held_before(hts_pos_t pos);
    void emit(bam1_t* b);

    Options opts_;
    std::unique_ptr<SamFile> in_;
    std::unique_ptr<bam_hdr_t, HeaderDeleter> in_header_;
    std::unique_ptr<hts_idx_t, IndexDeleter> in_index_;
    std::unique_ptr<faidx_t, FaiDeleter> fai_;
    samFile* out_ = nullptr;
    std::unique_ptr<bam_hdr_t, HeaderDeleter> out_header_;
    // The reference of the contig being written, uppercase (0-based), for left-alignment and tags.
    std::string contig_seq_;
    int contig_seq_tid_ = -1;
    // Records held until no later input record can be written before them (a min-heap by position).
    std::vector<Held> held_;
    std::size_t order_ = 0;
    int contig_tid_ = -1;        // the contig written last
    hts_pos_t written_to_ = 0;   // records starting at or before this, on it, were written with an earlier span
    int last_tid_ = -1;
    hts_pos_t last_pos_ = -1;
    std::size_t n_records_ = 0;
    std::size_t n_refined_ = 0;
    std::size_t n_left_aligned_ = 0;
    bool finished_ = false;
};

}  // namespace pgphase_collect

#endif
