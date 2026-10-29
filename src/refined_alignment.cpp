#include "refined_alignment.hpp"
#include "aligned_read_writer.hpp"

#include "version.hpp"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>

namespace pgphase_collect {

namespace {

bool consumes_query(uint32_t op) { return (bam_cigar_type(op) & 1) != 0; }
bool consumes_reference(uint32_t op) { return (bam_cigar_type(op) & 2) != 0; }

void push_op(std::vector<uint32_t>& cigar, uint32_t op, uint32_t len) {
    if (len == 0) return;
    if (!cigar.empty() && bam_cigar_op(cigar.back()) == op)
        cigar.back() = bam_cigar_gen(bam_cigar_oplen(cigar.back()) + len, op);
    else
        cigar.push_back(bam_cigar_gen(len, op));
}

// The CIGAR with = and X written as M.
std::vector<uint32_t> as_match(const std::vector<uint32_t>& cigar) {
    std::vector<uint32_t> out;
    out.reserve(cigar.size());
    for (const uint32_t c : cigar) {
        const uint32_t op = bam_cigar_op(c);
        push_op(out, op == BAM_CEQUAL || op == BAM_CDIFF ? BAM_CMATCH : op, bam_cigar_oplen(c));
    }
    return out;
}

uint32_t cigar_op_of(DigarType type) {
    switch (type) {
        case DigarType::Equal: return BAM_CEQUAL;
        case DigarType::Snp: return BAM_CDIFF;
        case DigarType::Insertion: return BAM_CINS;
        case DigarType::Deletion: return BAM_CDEL;
        case DigarType::SoftClip: return BAM_CSOFT_CLIP;
        case DigarType::HardClip: return BAM_CHARD_CLIP;
        case DigarType::RefSkip: return BAM_CREF_SKIP;
    }
    return BAM_CEQUAL;
}

// The read's alignment as its digars now describe it, when the MSA moved it and it holds up: every = base the
// reference's, every X base not (N either way), the query consumed exactly. -1 rejected, 0 unchanged, 1 refined.
int refine_one(const ReadRecord& r, const PhasingChunk& chunk, RefinedAlignment& out) {
    const bam1_t* b = r.alignment.get();
    std::vector<uint32_t> exact;
    exact.reserve(r.digars.size());
    hts_pos_t pos = -1;
    for (const DigarOp& d : r.digars) {
        const uint32_t op = cigar_op_of(d.type);
        if (pos < 0 && consumes_reference(op)) pos = d.pos - 1;
        if (d.len < 0) return -1;
        push_op(exact, op, static_cast<uint32_t>(d.len));
    }
    if (pos < 0) return -1;
    const uint32_t* in = bam_get_cigar(b);
    const std::vector<uint32_t> input(in, in + b->core.n_cigar);
    if (pos == b->core.pos && as_match(exact) == as_match(input)) return 0;

    // A hard clip (rewritten as a soft one, longcallD would write SEQ shorter than its CIGAR), a mate whose position
    // would go stale, or a start moved further than the writer holds records back: the input kept.
    const auto hard = [](const std::vector<uint32_t>& c) {
        return std::any_of(c.begin(), c.end(), [](uint32_t x) { return bam_cigar_op(x) == BAM_CHARD_CLIP; });
    };
    if (hard(input) || hard(exact) || (b->core.flag & BAM_FPAIRED) != 0 ||
        std::llabs(static_cast<long long>(pos - b->core.pos)) > kRefineMaxShift)
        return -1;

    const auto ref_at = [&chunk](hts_pos_t pos1) -> int {
        const hts_pos_t i = pos1 - chunk.ref_beg;
        if (i < 0 || i >= static_cast<hts_pos_t>(chunk.ref_seq.size())) return -1;
        return std::toupper(static_cast<unsigned char>(chunk.ref_seq[static_cast<size_t>(i)]));
    };
    const uint8_t* seq = bam_get_seq(b);
    int qp = 0;
    hts_pos_t rp = pos + 1;
    for (const uint32_t c : exact) {
        const uint32_t op = bam_cigar_op(c);
        const int len = static_cast<int>(bam_cigar_oplen(c));
        if (op == BAM_CEQUAL || op == BAM_CDIFF) {
            if (qp + len > b->core.l_qseq) return -1;
            for (int j = 0; j < len; ++j) {
                const int rb = ref_at(rp + j);
                const int qb = seq_nt16_str[bam_seqi(seq, qp + j)];
                if (rb < 0) return -1;
                if (rb == 'N' || qb == 'N') continue;
                if ((op == BAM_CEQUAL) != (rb == qb)) return -1;
            }
        } else if (op == BAM_CDEL || op == BAM_CREF_SKIP) {
            if (ref_at(rp) < 0 || ref_at(rp + len - 1) < 0) return -1;
        }
        if (consumes_query(op)) qp += len;
        if (consumes_reference(op)) rp += len;
    }
    if (qp != b->core.l_qseq) return -1;

    const bool input_eqx = std::any_of(input.begin(), input.end(), [](uint32_t x) {
        return bam_cigar_op(x) == BAM_CEQUAL || bam_cigar_op(x) == BAM_CDIFF;
    });
    out.qname = r.qname;
    out.orig_pos = b->core.pos;
    out.pos = pos;
    out.cigar = input_eqx ? std::move(exact) : as_match(exact);
    return 1;
}

// A record's CIGAR replaced in place.
void replace_cigar(bam1_t* b, const std::vector<uint32_t>& cigar) {
    const size_t qname = b->core.l_qname;
    const size_t old_bytes = static_cast<size_t>(b->core.n_cigar) * sizeof(uint32_t);
    const size_t new_bytes = cigar.size() * sizeof(uint32_t);
    const size_t rest = static_cast<size_t>(b->l_data) - qname - old_bytes;
    const size_t new_len = static_cast<size_t>(b->l_data) - old_bytes + new_bytes;
    if (new_len > b->m_data) {
        auto* data = static_cast<uint8_t*>(std::realloc(b->data, new_len));
        if (data == nullptr) throw std::bad_alloc();
        b->data = data;
        b->m_data = static_cast<uint32_t>(new_len);
    }
    std::memmove(b->data + qname + new_bytes, b->data + qname + old_bytes, rest);
    std::memcpy(b->data + qname, cigar.data(), new_bytes);
    b->core.n_cigar = static_cast<uint32_t>(cigar.size());
    b->l_data = static_cast<int>(new_len);
}

// NM, MD and minimap2's short cs of a record's alignment against `ref` (the reference from the record's start).
// false when the alignment leaves `ref`, or (cs) skips reference bases.
bool alignment_tags(const bam1_t* b, const char* ref, hts_pos_t ref_len, int& nm, std::string& md, std::string& cs,
                    bool& cs_ok) {
    const uint8_t* seq = bam_get_seq(b);
    const uint32_t* cigar = bam_get_cigar(b);
    hts_pos_t rp = 0;
    int qp = 0, md_run = 0, cs_run = 0;
    nm = 0;
    md.clear();
    cs.clear();
    cs_ok = true;
    const auto lower = [](int c) { return static_cast<char>(std::tolower(c)); };
    const auto flush_cs = [&]() {
        if (cs_run > 0) cs += ":" + std::to_string(cs_run);
        cs_run = 0;
    };
    for (uint32_t k = 0; k < b->core.n_cigar; ++k) {
        const uint32_t op = bam_cigar_op(cigar[k]);
        const int len = static_cast<int>(bam_cigar_oplen(cigar[k]));
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            if (rp + len > ref_len) return false;
            for (int j = 0; j < len; ++j) {
                const int rb = std::toupper(static_cast<unsigned char>(ref[rp + j]));
                const int qb = seq_nt16_str[bam_seqi(seq, qp + j)];
                if (rb == qb) {
                    ++md_run;
                    ++cs_run;
                } else {
                    ++nm;
                    md += std::to_string(md_run);
                    md += static_cast<char>(rb);
                    md_run = 0;
                    flush_cs();
                    cs += '*';
                    cs += lower(rb);
                    cs += lower(qb);
                }
            }
            rp += len;
            qp += len;
        } else if (op == BAM_CINS) {
            nm += len;
            flush_cs();
            cs += '+';
            for (int j = 0; j < len; ++j) cs += lower(seq_nt16_str[bam_seqi(seq, qp + j)]);
            qp += len;
        } else if (op == BAM_CDEL) {
            if (rp + len > ref_len) return false;
            nm += len;
            md += std::to_string(md_run) + "^";
            md_run = 0;
            flush_cs();
            cs += '-';
            for (int j = 0; j < len; ++j) {
                md += static_cast<char>(std::toupper(static_cast<unsigned char>(ref[rp + j])));
                cs += lower(ref[rp + j]);
            }
            rp += len;
        } else if (op == BAM_CREF_SKIP) {
            cs_ok = false;
            rp += len;
        } else if (op == BAM_CSOFT_CLIP) {
            qp += len;
        }
    }
    md += std::to_string(md_run);
    flush_cs();
    return true;
}

// Each insertion and deletion moved left while the base before it is its last one (a deletion: the reference's; an
// insertion: the read's), so no base changes its match status (minimap2 and pbmm2 place gaps so; vg's surjected
// alignments do not, and DeepVariant, trained on pbmm2's, called 276 more easy chr20 SNPs NoCall from them). A gap
// keeps one aligned base before it: the start does not move, and no gap meets another or a clip. ref holds the
// reference from 0-based position 0.
bool left_align_cigar(std::vector<uint32_t>& cigar, const uint8_t* seq, hts_pos_t pos, const std::string& ref) {
    bool moved = false;
    hts_pos_t rp = pos;
    int qp = 0;
    const auto matching = [](uint32_t op) { return op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF; };
    for (size_t k = 0; k < cigar.size(); ++k) {
        const uint32_t op = bam_cigar_op(cigar[k]);
        const int n = static_cast<int>(bam_cigar_oplen(cigar[k]));
        if ((op == BAM_CDEL || op == BAM_CINS) && k > 0 && matching(bam_cigar_op(cigar[k - 1]))) {
            const int before = static_cast<int>(bam_cigar_oplen(cigar[k - 1]));
            int s = 0;
            while (s < before - 1) {
                bool same;
                if (op == BAM_CDEL) {
                    const hts_pos_t a = rp - 1 - s, b = rp + n - 1 - s;
                    same = a >= 0 && b < static_cast<hts_pos_t>(ref.size()) &&
                           ref[static_cast<size_t>(a)] == ref[static_cast<size_t>(b)];
                } else {
                    const int a = qp - 1 - s, b = qp + n - 1 - s;
                    same = a >= 0 && bam_seqi(seq, a) == bam_seqi(seq, b);
                }
                if (!same) break;
                ++s;
            }
            if (s > 0) {
                moved = true;
                const uint32_t prev = bam_cigar_op(cigar[k - 1]);
                cigar[k - 1] = bam_cigar_gen(static_cast<uint32_t>(before - s), prev);
                // The bases moved past the gap keep their op (= or X stay what they were).
                if (k + 1 < cigar.size() && bam_cigar_op(cigar[k + 1]) == prev)
                    cigar[k + 1] = bam_cigar_gen(bam_cigar_oplen(cigar[k + 1]) + static_cast<uint32_t>(s), prev);
                else
                    cigar.insert(cigar.begin() + static_cast<long>(k) + 1, bam_cigar_gen(static_cast<uint32_t>(s), prev));
                rp -= s;
                qp -= s;
            }
        }
        if (consumes_reference(op)) rp += n;
        if (consumes_query(op)) qp += n;
    }
    return moved;
}

// The record's HP/PS: its read's union label, or none.
void set_label_tags(bam1_t* b, int hap, hts_pos_t ps) {
    if (uint8_t* t = bam_aux_get(b, "HP")) bam_aux_del(b, t);
    if (uint8_t* t = bam_aux_get(b, "PS")) bam_aux_del(b, t);
    if ((hap != 1 && hap != 2) || ps <= 0) return;
    bam_aux_update_int(b, "HP", hap);
    bam_aux_update_int(b, "PS", ps);
}

}  // namespace

std::vector<RefinedAlignment> take_refined_alignments(PhasingChunk& chunk, RefineCounts& counts) {
    std::vector<RefinedAlignment> out;
    for (ReadRecord& r : chunk.reads) {
        const bam1_t* b = r.alignment.get();
        if (b != nullptr && !r.digars.empty() &&
            (b->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) == 0) {
            RefinedAlignment refined;
            switch (refine_one(r, chunk, refined)) {
                case 1:
                    out.push_back(std::move(refined));
                    ++counts.refined;
                    break;
                case 0: ++counts.unchanged; break;
                default: ++counts.rejected; break;
            }
        }
        // What the rewrite kept past the alignment solve (mid_free_chunk frees it otherwise).
        r.digars.clear();
        r.digars.shrink_to_fit();
        r.qual.clear();
        r.qual.shrink_to_fit();
        r.noisy_regions.clear();
        r.noisy_regions.shrink_to_fit();
    }
    return out;
}

AlignedReadWriter::AlignedReadWriter(const Options& opts) : opts_(opts) {
    if (opts_.bam_files.size() != 1)
        throw std::runtime_error("aligned output (--out-bam/--out-cram/--out-sam) takes exactly one --bam");
    const std::string& path = opts_.bam_files.front();
    in_ = std::make_unique<SamFile>(path, opts_.threads, opts_.ref_fasta);
    in_header_.reset(sam_hdr_read(in_->get()));
    if (!in_header_) throw std::runtime_error("failed to read the header of " + path);
    in_index_.reset(sam_index_load(in_->get(), path.c_str()));
    if (!in_index_) throw std::runtime_error("aligned output needs an indexed --bam: " + path);
    fai_.reset(load_reference_index(opts_.ref_fasta));

    const OutputAlignmentFormat format = opts_.output_union_aln_format;
    out_ = sam_open(opts_.output_union_aln.c_str(),
                    format == OutputAlignmentFormat::Sam ? "w" : format == OutputAlignmentFormat::Cram ? "wc" : "wb");
    if (out_ == nullptr) throw std::runtime_error("failed to open alignment output: " + opts_.output_union_aln);
    if (format == OutputAlignmentFormat::Cram && hts_set_fai_filename(out_, opts_.ref_fasta.c_str()) != 0)
        throw std::runtime_error("failed to set the reference for CRAM output: " + opts_.ref_fasta);
    if (opts_.threads > 1) hts_set_threads(out_, opts_.threads);
    out_header_.reset(sam_hdr_dup(in_header_.get()));
    if (!out_header_ || sam_hdr_add_pg(out_header_.get(), "pgphase", "VN", PGPHASE_VERSION, "CL",
                                       opts_.command_line.c_str(), nullptr) != 0 ||
        sam_hdr_write(out_, out_header_.get()) != 0)
        throw std::runtime_error("failed to write the header of " + opts_.output_union_aln);
}

AlignedReadWriter::~AlignedReadWriter() {
    if (out_ != nullptr) sam_close(out_);
}

void AlignedReadWriter::check_region_order(
    const std::vector<std::pair<std::string, std::pair<hts_pos_t, hts_pos_t>>>& regions) const {
    int last_tid = -1;
    hts_pos_t last_end = 0;
    for (const auto& [contig, span] : regions) {
        const int tid = sam_hdr_name2tid(in_header_.get(), contig.c_str());
        if (tid < 0) throw std::runtime_error("aligned output: contig " + contig + " is not in the --bam header");
        if (tid < last_tid || (tid == last_tid && span.first <= last_end))
            throw std::runtime_error("aligned output: the regions (" + contig + ":" + std::to_string(span.first) +
                                     ") are not in the --bam header's order, or overlap; the output could not be "
                                     "sorted");
        last_tid = tid;
        last_end = span.second;
    }
}

const std::string& AlignedReadWriter::contig_sequence(int tid, const std::string& contig) {
    if (tid != contig_seq_tid_) {
        hts_pos_t len = 0;
        char* seq = faidx_fetch_seq64(fai_.get(), contig.c_str(), 0, HTS_POS_MAX, &len);
        if (seq == nullptr || len <= 0) {
            std::free(seq);
            throw std::runtime_error("aligned output: cannot read " + contig + " from " + opts_.ref_fasta);
        }
        contig_seq_.assign(seq, static_cast<size_t>(len));
        std::free(seq);
        for (char& c : contig_seq_) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
        contig_seq_tid_ = tid;
    }
    return contig_seq_;
}

void AlignedReadWriter::apply(bam1_t* b, const RefinedAlignment& r) {
    replace_cigar(b, r.cigar);
    b->core.pos = r.pos;
    b->core.bin = static_cast<uint16_t>(hts_reg2bin(b->core.pos, bam_endpos(b), 14, 5));
    ++n_refined_;
}

bool AlignedReadWriter::left_align(bam1_t* b) {
    const uint32_t* in = bam_get_cigar(b);
    std::vector<uint32_t> cigar(in, in + b->core.n_cigar);
    if (!left_align_cigar(cigar, bam_get_seq(b), b->core.pos, contig_seq_)) return false;
    replace_cigar(b, cigar);
    ++n_left_aligned_;
    return true;
}

void AlignedReadWriter::update_tags(bam1_t* b) {
    uint8_t* nm_tag = bam_aux_get(b, "NM");
    uint8_t* md_tag = bam_aux_get(b, "MD");
    uint8_t* cs_tag = bam_aux_get(b, "cs");
    if (nm_tag == nullptr && md_tag == nullptr && cs_tag == nullptr) return;
    const hts_pos_t end = bam_endpos(b);
    int nm = 0;
    std::string md, cs;
    bool cs_ok = false;
    const bool ok = b->core.pos >= 0 && end <= static_cast<hts_pos_t>(contig_seq_.size()) &&
                    alignment_tags(b, contig_seq_.data() + b->core.pos, end - b->core.pos, nm, md, cs, cs_ok);
    // Tags that would no longer describe the alignment are dropped rather than left stale.
    if (nm_tag != nullptr) {
        if (ok) bam_aux_update_int(b, "NM", nm);
        else bam_aux_del(b, bam_aux_get(b, "NM"));
    }
    if (md_tag != nullptr) {
        if (ok) bam_aux_update_str(b, "MD", static_cast<int>(md.size() + 1), md.c_str());
        else bam_aux_del(b, bam_aux_get(b, "MD"));
    }
    if (cs_tag != nullptr) {
        // A long cs (with = runs spelled out) is rewritten short, as minimap2's default.
        if (ok && cs_ok) bam_aux_update_str(b, "cs", static_cast<int>(cs.size() + 1), cs.c_str());
        else bam_aux_del(b, bam_aux_get(b, "cs"));
    }
}

void AlignedReadWriter::emit(bam1_t* b) {
    const int tid = b->core.tid;
    if (tid < last_tid_ || (tid == last_tid_ && b->core.pos < last_pos_))
        throw std::logic_error("aligned output out of order at " + std::string(bam_get_qname(b)));
    last_tid_ = tid;
    last_pos_ = b->core.pos;
    if (sam_write1(out_, out_header_.get(), b) < 0)
        throw std::runtime_error("failed to write " + std::string(bam_get_qname(b)) + " to " + opts_.output_union_aln);
    ++n_records_;
}

void AlignedReadWriter::write_held_before(hts_pos_t pos) {
    const auto later = [](const Held& x, const Held& y) { return x.pos != y.pos ? x.pos > y.pos : x.order > y.order; };
    while (!held_.empty() && held_.front().pos < pos) {
        std::pop_heap(held_.begin(), held_.end(), later);
        emit(held_.back().b.get());
        held_.pop_back();
    }
}

void AlignedReadWriter::write_contig(const std::string& contig,
                                     const std::vector<std::pair<hts_pos_t, hts_pos_t>>& spans,
                                     const std::unordered_map<std::string, const RefinedAlignment*>& refined,
                                     const ReadLabel& label) {
    const int tid = sam_hdr_name2tid(in_header_.get(), contig.c_str());
    if (tid < 0) throw std::runtime_error("aligned output: contig " + contig + " is not in the --bam header");
    if (tid != contig_tid_) {
        write_held_before(std::numeric_limits<hts_pos_t>::max());
        contig_tid_ = tid;
        written_to_ = 0;
    }
    contig_sequence(tid, contig);
    const auto later = [](const Held& x, const Held& y) { return x.pos != y.pos ? x.pos > y.pos : x.order > y.order; };
    std::unique_ptr<bam1_t, AlignmentDeleter> rec(bam_init1());
    for (const auto& [beg, end] : spans) {
        std::unique_ptr<hts_itr_t, IteratorDeleter> itr(sam_itr_queryi(in_index_.get(), tid, beg - 1, end));
        if (!itr) throw std::runtime_error("aligned output: cannot query " + contig + " in the --bam");
        int ret = 0;
        while ((ret = sam_itr_next(in_->get(), itr.get(), rec.get())) >= 0) {
            if (rec->core.pos + 1 <= written_to_) continue;
            const hts_pos_t input_pos = rec->core.pos;
            std::unique_ptr<bam1_t, AlignmentDeleter> b(bam_dup1(rec.get()));
            if (!b) throw std::bad_alloc();
            bool rewritten = false;
            if ((b->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) == 0) {
                const std::string qname = bam_get_qname(b.get());
                const auto it = refined.find(qname);
                if (it != refined.end() && it->second->orig_pos == b->core.pos) {
                    apply(b.get(), *it->second);
                    rewritten = true;
                }
                const auto [hap, ps] = label(qname);
                set_label_tags(b.get(), hap, ps);
            } else {
                set_label_tags(b.get(), 0, 0);
            }
            if (opts_.left_align_aln && (b->core.flag & BAM_FUNMAP) == 0 && b->core.n_cigar > 0 && b->core.l_qseq > 0)
                rewritten = left_align(b.get()) || rewritten;
            if (rewritten) update_tags(b.get());
            const hts_pos_t pos = b->core.pos;
            held_.push_back(Held{pos, order_++, std::move(b)});
            std::push_heap(held_.begin(), held_.end(), later);
            // A refinement moves a start by at most kRefineMaxShift: no later input record can be written before
            // one held that far behind this one.
            write_held_before(input_pos - kRefineMaxShift);
        }
        if (ret < -1) throw std::runtime_error("aligned output: failed reading " + contig + " from the --bam");
        written_to_ = std::max(written_to_, end);
    }
}

void AlignedReadWriter::finish() {
    if (finished_) return;
    finished_ = true;
    write_held_before(std::numeric_limits<hts_pos_t>::max());
    if (sam_close(out_) != 0) {
        out_ = nullptr;
        throw std::runtime_error("failed to close alignment output: " + opts_.output_union_aln);
    }
    out_ = nullptr;
    if (opts_.output_union_aln_format != OutputAlignmentFormat::Sam &&
        sam_index_build3(opts_.output_union_aln.c_str(), nullptr, 0, opts_.threads) != 0)
        throw std::runtime_error("failed to index alignment output: " + opts_.output_union_aln);
}

}  // namespace pgphase_collect
