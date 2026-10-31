// Union gap phasing, stage 1b: the injection plan (union_select.cpp) applied -- the alignment solve's selected
// sites enter the graph chunk (injection, hidden complement sites), its rows flagged as the plan says.

#include "union_select.hpp"

#include "noise_filter.hpp"
#include "trace.hpp"

#include <algorithm>
#include <cstdlib>
#include <optional>
#include <map>
#include <set>

namespace pgphase_collect {

namespace {

using CandKey = CandidateIdentityKey;

CandKey cand_key_of(const CandidateVariant& cand) {
    return CandKey{cand.key.sort_pos(), static_cast<int>(cand.key.type), cand.key.ref_len,
                   cand.key.alt};
}

}  // namespace

/// VCF-form metadata for an alignment candidate the catalog never held, in the
/// convention the graph writer emits and vcf_to_variant_key reads back. An
/// empty REF makes the writer skip the record while index-parallel arrays stay
/// aligned.
static std::pair<GraphSiteMeta, std::vector<int>> synthesize_alignment_site_meta(
        const CandidateVariant& cand, const PhasingChunk& source, const char* contig_name) {
    const auto ref_at = [&source](hts_pos_t pos, int len) -> std::string {
        if (len <= 0 || source.ref_seq.empty()) return std::string();
        const hts_pos_t off = pos - source.ref_beg;
        if (off < 0 || static_cast<size_t>(off) + static_cast<size_t>(len) > source.ref_seq.size())
            return std::string();
        return source.ref_seq.substr(static_cast<size_t>(off), static_cast<size_t>(len));
    };
    const VariantKey& key = cand.key;
    const auto anchored = [&](const std::string& allele) -> std::string {
        if (key.type == VariantType::Snp) return allele;
        if (key.type == VariantType::Insertion && key.ref_len != 0) return allele;
        return ref_at(key.pos - 1, 1) + allele;
    };
    hts_pos_t vcf_pos = key.pos;
    std::string vcf_ref;
    if (key.type == VariantType::Snp) {
        vcf_ref = ref_at(key.pos, std::max(1, key.ref_len));
    } else if (key.type == VariantType::Insertion && key.ref_len == 0) {
        vcf_pos = key.pos - 1;
        vcf_ref = ref_at(vcf_pos, 1);
    } else if (key.type == VariantType::Insertion) {
        vcf_ref = ref_at(key.pos, key.ref_len);
    } else {
        vcf_pos = key.pos - 1;
        vcf_ref = ref_at(vcf_pos, 1) + ref_at(key.pos, key.ref_len);
    }
    const std::string vcf_alt = anchored(key.alt);
    bool usable = !vcf_ref.empty() && !vcf_alt.empty();
    if (usable) {
        const VariantKey round_trip = vcf_to_variant_key(key.tid, vcf_pos, vcf_ref, vcf_alt);
        usable = round_trip.type == key.type && round_trip.pos == key.pos &&
                 round_trip.ref_len == key.ref_len && round_trip.alt == key.alt;
    }
    // An MSA indel whose anchor letter is not the reference's on its haplotype (a SNP right before an insertion or a
    // deletion, aligned as the mismatch then the gap): the consensus's letter there leads the ALT, one complex record,
    // as collect-bam-variation writes it (alt_ref_base). With the reference's letter the SNP was lost: chr12:6,396,535
    // C>CG beside A>G written as A>AG, another haplotype. HG002 HiFi chr12/18/20: true het regions 29,979 -> 30,200,
    // wrong 681 -> 487; reads unchanged.
    std::string written_alt = vcf_alt;
    if (usable && key.type != VariantType::Snp && cand.msa_insertion_alts.empty() &&
        cand.alt_ref_base < 4 && !written_alt.empty()) {
        const char consensus = "ACGT"[cand.alt_ref_base];
        if (std::toupper(static_cast<unsigned char>(written_alt[0])) != consensus) written_alt[0] = consensus;
    }
    std::vector<std::string> alts;
    if (cand.msa_insertion_alts.size() >= 2) {
        for (size_t ai = 0; ai < cand.msa_insertion_alts.size(); ++ai) {
            std::string a = anchored(cand.msa_insertion_alts[ai]);
            if (a.empty()) { alts.clear(); break; }
            // Each allele's own consensus letter leads it, as for a single ALT above.
            const uint8_t anchor = msa_allele_anchor(cand, ai);
            if (key.type != VariantType::Snp && !(key.type == VariantType::Insertion && key.ref_len != 0) && anchor < 4)
                a[0] = "ACGT"[anchor];
            alts.push_back(a);
        }
    }
    if (alts.empty() && usable) alts.push_back(written_alt);
    GraphSiteMeta meta;
    meta.chrom = contig_name != nullptr ? contig_name : std::string();
    if (!usable || alts.empty()) pgphase_trace::candidate("inject", "drop", cand, "unusable_vcf_meta");
    if (usable && !alts.empty()) {
        meta.pos = vcf_pos;
        meta.ref = vcf_ref;
        meta.alts = alts;
    }
    std::vector<int> orig;
    for (int i = 0; i <= static_cast<int>(alts.size()); ++i) orig.push_back(i);
    if (orig.size() < 2) orig = std::vector<int>{0, 1};
    return {std::move(meta), std::move(orig)};
}

// Left-aligned identity of a simple insertion or deletion, so placements of
// one allele at different repeat copies compare equal. Other rows keep their key.
static CandKey left_aligned_indel_key(const CandidateVariant& c, const PhasingChunk& source) {
    VariantKey k = c.key;
    const auto base = [&source](hts_pos_t pos) -> char {
        const hts_pos_t off = pos - source.ref_beg;
        if (off < 0 || static_cast<size_t>(off) >= source.ref_seq.size()) return 'N';
        return static_cast<char>(std::toupper(static_cast<unsigned char>(source.ref_seq[static_cast<size_t>(off)])));
    };
    if (!c.msa_insertion_alts.empty()) return cand_key_of(c);
    // An indel whose anchor letter is not the reference's (alt_ref_base: a SNP right before it, written as one
    // complex record) is another allele than the plain indel it would left-align onto: keyed with its letter, it is
    // never folded into one. chr20:9,869,936: one consensus deleted a T of T8 with the T>G before it as its anchor,
    // the other the plain T; folded, the SNP was lost and AT>A written het where the truth has it homozygous beside
    // T>G.
    if (k.type != VariantType::Snp && c.alt_ref_base < 4 && "ACGT"[c.alt_ref_base] != base(k.pos - 1)) {
        CandKey key = cand_key_of(c);
        key.alt = std::string("^") + "ACGT"[c.alt_ref_base] + key.alt;
        return key;
    }
    if (k.type == VariantType::Insertion && k.ref_len == 0 && !k.alt.empty()) {
        // The inserted bases sit after pos - 1; shift while the base before
        // the insertion point equals the last inserted base.
        while (true) {
            const char prev = base(k.pos - 1);
            if (prev == 'N' || std::toupper(static_cast<unsigned char>(k.alt.back())) != prev) break;
            k.alt = std::string(1, prev) + k.alt.substr(0, k.alt.size() - 1);
            --k.pos;
        }
    } else if (k.type == VariantType::Deletion && k.ref_len > 0 && k.alt.empty()) {
        while (true) {
            const char prev = base(k.pos - 1);
            if (prev == 'N' || prev != base(k.pos + k.ref_len - 1)) break;
            --k.pos;
        }
    }
    return CandKey{k.sort_pos(), static_cast<int>(k.type), k.ref_len, k.alt};
}

namespace {

// One site per left-aligned allele: other placements become aliases whose
// read calls fold into the canonical row. Regular sites come first: an
// island-only call never absorbs one, and one at an allele a regular site
// already carries is dropped. Returns alias -> canonical candidate.
std::unordered_map<size_t, size_t> fold_aliases(Selection& sel, const PhasingChunk& bam) {
    std::unordered_map<size_t, size_t> alias_of;
    std::map<CandKey, size_t> canonical;
    std::vector<size_t> kept;
    std::stable_partition(sel.selected.begin(), sel.selected.end(),
                          [&](size_t bi) { return sel.island_only.count(bi) == 0; });
    for (const size_t bi : sel.selected) {
        const auto [it, inserted] = canonical.emplace(left_aligned_indel_key(bam.candidates[bi], bam), bi);
        if (inserted) {
            kept.push_back(bi);
        } else if (sel.island_only.count(bi) == 0) {
            alias_of.emplace(bi, it->second);
            pgphase_trace::candidate("inject", "merge", bam.candidates[bi], "alias_folded",
                                     "into=" + pgphase_trace::key_text(bam.candidates[it->second].key));
        } else {
            pgphase_trace::candidate("inject", "drop", bam.candidates[bi], "island_duplicate_dropped",
                                     "of=" + pgphase_trace::key_text(bam.candidates[it->second].key));
        }
    }
    sel.selected = std::move(kept);
    return alias_of;
}

}  // namespace

size_t inject_alignment_sites(GraphChunkBuildResult& graph_chunk, const PhasingChunk& bam,
                              const char* contig_name) {
    PhasingChunk& graph = graph_chunk.chunk;
    static const char* const diag_dir = std::getenv("PGPHASE_DIAG");
    Selection sel = plan_injection(graph_chunk, bam, diag_dir != nullptr);
    for (const size_t ci : sel.unwritten_rows)
        if (ci < graph_chunk.site_meta.size()) graph_chunk.site_meta[ci].no_record = true;
    for (const size_t ci : sel.hom_rows)
        if (ci < graph_chunk.site_meta.size()) graph_chunk.site_meta[ci].msa_hom = true;
    if (pgphase_trace::on() && graph_chunk.site_meta.size() == graph.candidates.size()) {
        const auto row = [&](size_t ci, std::string_view event, std::string_view reason) {
            const GraphSiteMeta& m = graph_chunk.site_meta[ci];
            std::string key = m.ref + ">";
            for (size_t i = 0; i < m.alts.size(); ++i) key += (i ? "," : "") + m.alts[i].substr(0, 60);
            pgphase_trace::emit("inject", event, m.pos, m.pos + std::max<hts_pos_t>(1, m.ref.size()) - 1, key, reason,
                                ci < graph_chunk.site_ids.size() ? "site=" + graph_chunk.site_ids[ci] : std::string());
        };
        for (const auto& [ci, reason] : sel.row_reason) row(ci, sel.hom_rows.count(ci) ? "regenotype" : "drop", reason);
        for (const size_t ci : sel.replaced_rows) row(ci, "drop", "row_replaced_by_msa_call");
    }
    if (diag_dir != nullptr) write_injection_diagnostics(diag_dir, graph_chunk, bam, sel);
    if (sel.selected.empty()) return 0;
    const std::unordered_map<size_t, size_t> alias_of = fold_aliases(sel, bam);

    // A catalog row its complement call describes exactly (the call's VCF allele the row's ALT) takes the call's phase
    // in its alignment solve, as an injected call keeps its own (CandidateVariant::aln_phase_set): what
    // phase_by_alignment_solve relates sites by. Catalog SNPs the EM phases at error 0.1-0.3 are written unphased
    // (chr20:1,388,819 G>T at 0.12, its complement call phased by the alignment solve beside the insertion at 816).
    std::unordered_map<size_t, std::pair<hts_pos_t, int>> row_aln_phase;
    if (graph_chunk.site_meta.size() == graph.candidates.size() &&
        graph_chunk.site_allele_orig_idx.size() == graph.candidates.size()) {
        const auto base = [&bam](hts_pos_t pos) -> char {
            const hts_pos_t off = pos - bam.ref_beg;
            if (off < 0 || static_cast<size_t>(off) >= bam.ref_seq.size()) return 'N';
            return static_cast<char>(std::toupper(static_cast<unsigned char>(bam.ref_seq[static_cast<size_t>(off)])));
        };
        for (const auto& [bi, rows] : sel.complement_rows) {
            if (rows.size() != 1) continue;
            const CandidateVariant& c = bam.candidates[bi];
            const auto& h = c.hap_to_cons_alle;
            if (c.phase_set <= 0 || !c.msa_insertion_alts.empty() || h[1] < 0 || h[2] < 0 || h[1] + h[2] != 1) continue;
            std::string call_ref;
            for (int i = 0; i < c.key.ref_len; ++i) call_ref.push_back(base(c.key.pos + i));
            const auto call = normalize_candidate_identity(c.key.pos, call_ref, c.key.alt, base);
            const size_t r = rows[0];
            const std::optional<RowAllele> allele = row_allele(graph_chunk, r);
            if (!call || !allele) continue;
            const auto row = normalize_candidate_identity(allele->pos, *allele->ref, *allele->alt, base);
            if (row && *row == *call) row_aln_phase.emplace(r, std::make_pair(c.phase_set, h[1]));
        }
    }

    // Merge in key order; every index-parallel array is rebuilt with the table.
    struct Slot { CandKey key; long graph_index; long bam_index; };
    std::vector<Slot> slots;
    slots.reserve(graph.candidates.size() + sel.selected.size());
    for (size_t ci = 0; ci < graph.candidates.size(); ++ci)
        slots.push_back(Slot{cand_key_of(graph.candidates[ci]), static_cast<long>(ci), -1});
    for (const size_t bi : sel.selected) {
        CandidateVariant relabeled = bam.candidates[bi];
        relabeled.key.tid = graph.region.tid;
        slots.push_back(Slot{cand_key_of(relabeled), -1, static_cast<long>(bi)});
    }
    std::stable_sort(slots.begin(), slots.end(),
                     [](const Slot& a, const Slot& b) { return a.key < b.key; });
    const bool had_meta = graph_chunk.site_meta.size() == graph.candidates.size();
    CandidateTable merged;
    std::vector<std::string> merged_ids;
    std::vector<GraphSiteMeta> merged_meta;
    std::vector<std::vector<int>> merged_orig;
    std::vector<long> graph_to_new(graph.candidates.size(), -1);
    std::unordered_map<size_t, size_t> bam_to_new;
    for (const Slot& slot : slots) {
        const size_t ni = merged.size();
        if (slot.graph_index >= 0) {
            const size_t ci = static_cast<size_t>(slot.graph_index);
            graph_to_new[ci] = static_cast<long>(ni);
            merged.push_back(graph.candidates[ci]);
            if (const auto it = row_aln_phase.find(ci); it != row_aln_phase.end()) {
                merged.back().aln_phase_set = it->second.first;
                merged.back().aln_hap1_allele = it->second.second;
            }
            merged_ids.push_back(ci < graph_chunk.site_ids.size() ? graph_chunk.site_ids[ci] : std::string());
            merged_meta.push_back(had_meta ? graph_chunk.site_meta[ci] : GraphSiteMeta{});
            if (sel.replaced_rows.count(ci) != 0) {
                merged_meta.back().ref.clear();  // empty REF: the writer emits no record
                merged_meta.back().alts.clear();
            }
            merged_orig.push_back(ci < graph_chunk.site_allele_orig_idx.size()
                                      ? graph_chunk.site_allele_orig_idx[ci] : std::vector<int>{});
            continue;
        }
        const size_t bi = static_cast<size_t>(slot.bam_index);
        bam_to_new.emplace(bi, ni);
        CandidateVariant cand = bam.candidates[bi];
        cand.key.tid = graph.region.tid;
        cand.bam_injected = true;
        {
            const auto& h = cand.hap_to_cons_alle;
            if (cand.phase_set > 0 && h.size() > 2 && h[1] >= 0 && h[2] >= 0 && h[1] != h[2]) {
                cand.aln_phase_set = cand.phase_set;
                cand.aln_hap1_allele = h[1];
            } else if (cand.msa_verified) {
                aln_phase_from_msa_clusters(bam, cand);
            }
            if (cand.msa_verified && h.size() > 2 && h[1] >= 0 && h[2] >= 0 && h[1] != h[2])
                for (size_t r = 0; r < bam.msa_regions.size(); ++r)
                    if (bam.msa_regions[r].het() && bam.msa_regions[r].beg < cand.key.pos &&
                        cand.key.pos < bam.msa_regions[r].end) {
                        cand.msa_region = static_cast<int>(r);
                        cand.msa_hap1_allele = h[1];
                        break;
                    }
        }
        cand.alignment_verified = true;
        // Phase comes from the union's own solve, not the alignment's gauge.
        cand.phase_set = kUnsetCandidatePhaseSet;
        cand.hap_to_cons_alle[1] = -1;
        cand.hap_to_cons_alle[2] = -1;
        cand.hap_alt = 0;
        cand.hap_ref = 0;
        if (cand.counts.alle_covs.size() < 2)
            cand.counts.alle_covs = {cand.counts.ref_cov, cand.counts.alt_cov};
        auto [meta, orig] = synthesize_alignment_site_meta(bam.candidates[bi], bam, contig_name);
        cand.low_mapq_only = sel.low_mapq_only.count(bi) != 0;
        if (sel.island_only.count(bi) != 0) {
            cand.island_only = true;
            cand.lcd_var_i_to_cate &= ~kCandGermlineVarCate;
        }
        if (sel.hom_records.count(bi) != 0) {
            cand.lcd_var_i_to_cate &= ~kCandGermlineVarCate;
            cand.hap_to_cons_alle[1] = 1;
            cand.hap_to_cons_alle[2] = 1;
            meta.msa_hom = true;
        }
        // A depth-gated MSA het: written unphased with the MSA's genotype (each consensus its own allele), if any.
        const std::array<int8_t, 2>& cons = cand.msa_cons_allele;
        const bool island_het = sel.island_written.count(bi) != 0 &&
                                (cand.msa_insertion_alts.size() == 2 ||
                                 (cons[0] >= 0 && cons[1] >= 0 && cons[0] + cons[1] == 1));
        if (island_het) {
            cand.hap_to_cons_alle[1] = cand.msa_insertion_alts.size() == 2 ? 1 : 0;
            cand.hap_to_cons_alle[2] = cand.msa_insertion_alts.size() == 2 ? 2 : 1;
            cand.vcf_unphased = true;
        }
        if (sel.complement_rows.count(bi) != 0 || (sel.island_only.count(bi) != 0 && !island_het)) {
            meta = GraphSiteMeta{};  // empty REF: the writer emits no record
            meta.chrom = contig_name != nullptr ? contig_name : std::string();
        }
        merged.push_back(std::move(cand));
        merged_ids.emplace_back();
        merged_meta.push_back(std::move(meta));
        merged_orig.push_back(std::move(orig));
    }

    // Alignment observations at the injected sites, keyed by read name. Reads
    // the GAF never had are added: they are what carries a catalog gap.
    std::unordered_map<std::string_view, size_t> graph_read_by_qname;
    graph_read_by_qname.reserve(graph.reads.size());
    for (size_t ri = 0; ri < graph.reads.size(); ++ri)
        graph_read_by_qname.try_emplace(graph.reads[ri].qname, ri);
    std::vector<std::vector<std::tuple<size_t, int, int>>> added_obs(graph.reads.size());
    const size_t old_profiles = graph.reads.size();
    std::unordered_map<size_t, const std::vector<size_t>*> complement_of_new;  // merged index -> graph rows
    for (const auto& [bi, rows] : sel.complement_rows) {
        const auto it = bam_to_new.find(bi);
        if (it != bam_to_new.end()) complement_of_new.emplace(it->second, &rows);
    }
    for (const ReadVariantProfile& profile : bam.read_var_profile) {
        if (profile.read_id < 0 || profile.start_var_idx < 0) continue;
        const ReadRecord& source_read = bam.reads[static_cast<size_t>(profile.read_id)];
        if (source_read.is_skipped) continue;
        std::vector<std::tuple<size_t, int, int>> calls;
        for (size_t k = 0; k < profile.alleles.size(); ++k) {
            if (profile.alleles[k] < 0) continue;
            size_t source_i = static_cast<size_t>(profile.start_var_idx) + k;
            const auto alias = alias_of.find(source_i);
            if (alias != alias_of.end()) source_i = alias->second;
            const auto found = bam_to_new.find(source_i);
            if (found == bam_to_new.end()) continue;
            if (sel.hom_records.count(source_i) != 0) continue;  // written 1/1, no read calls
            const int qi = k < profile.alt_qi.size() ? profile.alt_qi[k] : 0;
            // An alias and its canonical row describe one allele: ALT wins.
            const auto prior = std::find_if(calls.begin(), calls.end(),
                [&](const auto& call) { return std::get<0>(call) == found->second; });
            if (prior == calls.end())
                calls.emplace_back(found->second, profile.alleles[k], qi);
            else if (profile.alleles[k] > std::get<1>(*prior))
                *prior = std::make_tuple(found->second, profile.alleles[k], qi);
        }
        if (calls.empty()) continue;
        auto read_i = graph_read_by_qname.find(source_read.qname);
        size_t ri;
        if (read_i != graph_read_by_qname.end()) {
            ri = read_i->second;
        } else {
            ReadRecord read;
            read.tid = graph.region.tid;
            read.qname = source_read.qname;
            read.mapq = source_read.mapq;
            read.beg = source_read.beg;
            read.end = source_read.end;
            read.reverse = source_read.reverse;
            ri = graph.reads.size();
            graph.reads.push_back(std::move(read));
            graph.read_var_profile.push_back(ReadVariantProfile{});
            added_obs.emplace_back();
            graph_read_by_qname.try_emplace(graph.reads.back().qname, ri);
        }
        for (auto& call : calls) {
            const auto comp = complement_of_new.find(std::get<0>(call));
            // A hidden complement site takes calls only from reads the graph
            // left silent: the read speaks through the graph row already.
            if (comp != complement_of_new.end() && ri < old_profiles) {
                const ReadVariantProfile& old = graph.read_var_profile[ri];
                bool graph_call = false;
                for (const size_t row : *comp->second) {
                    if (old.start_var_idx < 0 || static_cast<int>(row) < old.start_var_idx ||
                        static_cast<int>(row) > old.end_var_idx) continue;
                    if (old.alleles[row - static_cast<size_t>(old.start_var_idx)] >= 0) graph_call = true;
                }
                if (graph_call) continue;
            }
            added_obs[ri].push_back(call);
        }
    }

    // Rebuild each profile over the merged index space, every channel included.
    for (size_t ri = 0; ri < graph.reads.size(); ++ri) {
        const ReadVariantProfile& old = graph.read_var_profile[ri];
        std::map<size_t, std::array<int, 5>> by_index;  // allele, qi, graph, bam, bam_qi
        const auto slot = [&by_index](size_t i) -> std::array<int, 5>& {
            return by_index.try_emplace(i, std::array<int, 5>{-1, 0, -1, -1, 0}).first->second;
        };
        std::map<size_t, uint8_t> qualities;
        if (old.start_var_idx >= 0) {
            for (size_t k = 0; k < old.alleles.size(); ++k) {
                const size_t ci = static_cast<size_t>(old.start_var_idx) + k;
                if (ci >= graph_to_new.size() || graph_to_new[ci] < 0) continue;
                const size_t ni = static_cast<size_t>(graph_to_new[ci]);
                const bool any = old.alleles[k] >= 0 ||
                    (k < old.graph_alleles.size() && old.graph_alleles[k] >= 0) ||
                    (k < old.bam_alleles.size() && old.bam_alleles[k] != -1);
                if (!any) continue;
                auto& s = slot(ni);
                s[0] = old.alleles[k];
                s[1] = k < old.alt_qi.size() ? old.alt_qi[k] : 0;
                s[2] = k < old.graph_alleles.size() ? old.graph_alleles[k] : -1;
                s[3] = k < old.bam_alleles.size() ? old.bam_alleles[k] : -1;
                s[4] = k < old.bam_qi.size() ? old.bam_qi[k] : 0;
                if (k < old.bam_base_qualities.size() && old.bam_base_qualities[k] > 0)
                    qualities.emplace(ni, old.bam_base_qualities[k]);
            }
        }
        for (const auto& [ni, allele, qi] : added_obs[ri]) {
            auto& s = slot(ni);
            s[0] = allele;
            s[1] = qi;
            s[3] = allele;
            s[4] = qi;
        }
        ReadVariantProfile prof;
        prof.read_id = static_cast<int>(ri);
        prof.bam_mapq = old.bam_mapq;
        if (!by_index.empty()) {
            prof.start_var_idx = static_cast<int>(by_index.begin()->first);
            prof.end_var_idx = static_cast<int>(by_index.rbegin()->first);
            const size_t span = static_cast<size_t>(prof.end_var_idx - prof.start_var_idx + 1);
            prof.alleles.assign(span, -1);
            prof.alt_qi.assign(span, 0);
            const bool graph_channel = !old.graph_alleles.empty();
            const bool bam_channel = !old.bam_alleles.empty() || !added_obs[ri].empty();
            if (graph_channel) prof.graph_alleles.assign(span, -1);
            if (bam_channel) { prof.bam_alleles.assign(span, -1); prof.bam_qi.assign(span, 0); }
            if (!qualities.empty()) prof.bam_base_qualities.assign(span, 0);
            for (const auto& [ni, s] : by_index) {
                const size_t off = ni - static_cast<size_t>(prof.start_var_idx);
                prof.alleles[off] = s[0];
                prof.alt_qi[off] = s[1];
                if (graph_channel) prof.graph_alleles[off] = s[2];
                if (bam_channel) { prof.bam_alleles[off] = s[3]; prof.bam_qi[off] = s[4]; }
            }
            for (const auto& [ni, q] : qualities)
                prof.bam_base_qualities[ni - static_cast<size_t>(prof.start_var_idx)] = q;
        }
        graph.read_var_profile[ri] = std::move(prof);
    }
    graph.candidates = std::move(merged);
    graph_chunk.site_ids = std::move(merged_ids);
    graph_chunk.site_meta = std::move(merged_meta);
    graph_chunk.site_allele_orig_idx = std::move(merged_orig);
    graph_chunk.candidate_allele_contrasts.assign(graph.candidates.size(), std::nullopt);
    graph.haps.resize(graph.reads.size(), 0);
    graph.phase_sets.resize(graph.reads.size(), kUnphasedReadPhaseSet);
    rebuild_read_var_cr(graph);
    return sel.selected.size();
}

} // namespace pgphase_collect
