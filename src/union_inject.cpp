// Union gap phasing, stage 1: the alignment solve's heterozygous sites enter
// the graph chunk (injection, hidden complement sites).

#include "union_internal.hpp"

namespace pgphase_collect {

namespace {

using CandKey = CandidateIdentityKey;

CandKey cand_key_of(const CandidateVariant& cand) {
    return CandKey{cand.key.sort_pos(), static_cast<int>(cand.key.type), cand.key.ref_len,
                   cand.key.alt};
}

}  // namespace

// A site is admitted only if each of its two alleles is seen by at least this
// many reads mapped at kDefaultMinMapq or better.
constexpr int kUnionMinConfidentAlleleReads = 2;

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
    std::vector<std::string> alts;
    if (cand.msa_insertion_alts.size() >= 2) {
        for (const std::string& allele : cand.msa_insertion_alts) {
            const std::string a = anchored(allele);
            if (a.empty()) { alts.clear(); break; }
            alts.push_back(a);
        }
    }
    if (alts.empty() && usable) alts.push_back(vcf_alt);
    GraphSiteMeta meta;
    meta.chrom = contig_name != nullptr ? contig_name : std::string();
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

size_t inject_alignment_sites(GraphChunkBuildResult& graph_chunk, const PhasingChunk& bam,
                              const char* contig_name) {
    PhasingChunk& graph = graph_chunk.chunk;

    // Rows that already vote describe this locus. An alignment call overlapping
    // one is another description of it: it joins as a hidden complement site
    // (no VCF record) carrying the alignment's calls, because the graph often
    // leaves reads at such a locus without any allele.
    std::vector<std::tuple<hts_pos_t, hts_pos_t, size_t>> voting_spans;
    for (size_t ci = 0; ci < graph.candidates.size(); ++ci) {
        const CandidateVariant& row = graph.candidates[ci];
        if ((row.lcd_var_i_to_cate & kCandGermlineVarCate) == 0) continue;
        voting_spans.emplace_back(row.key.pos - 1, row.key.pos + std::max(1, row.key.ref_len), ci);
    }
    std::sort(voting_spans.begin(), voting_spans.end());
    hts_pos_t widest = 0;
    for (const auto& span : voting_spans) widest = std::max(widest, std::get<1>(span) - std::get<0>(span));
    const auto voting_rows_over = [&](hts_pos_t beg, hts_pos_t end) {
        std::vector<size_t> rows;
        auto it = std::lower_bound(voting_spans.begin(), voting_spans.end(),
                                   std::make_tuple(beg - widest, hts_pos_t{0}, size_t{0}));
        for (; it != voting_spans.end() && std::get<0>(*it) < end; ++it)
            if (std::get<1>(*it) > beg) rows.push_back(std::get<2>(*it));
        return rows;
    };
    std::unordered_map<size_t, std::vector<size_t>> complement_rows;  // bam candidate -> graph rows it overlaps

    // Reads admitted at low MAPQ can represent a divergent haplotype but cannot
    // create a site on their own: each of the site's two alleles (a 1/2 site
    // carries no REF reads) needs confidently mapped support.
    std::vector<std::map<int, int>> confident(bam.candidates.size());
    std::vector<std::map<int, int>> observed(bam.candidates.size());
    for (const ReadVariantProfile& prof : bam.read_var_profile) {
        if (prof.read_id < 0 || prof.start_var_idx < 0) continue;
        const ReadRecord& read = bam.reads[static_cast<size_t>(prof.read_id)];
        if (read.is_skipped) continue;
        for (size_t k = 0; k < prof.alleles.size(); ++k) {
            const int a = prof.alleles[k];
            if (a < 0) continue;
            const size_t ci = static_cast<size_t>(prof.start_var_idx) + k;
            ++observed[ci][a];
            if (read.mapq >= kDefaultMinMapq) ++confident[ci][a];
        }
    }
    const auto confidently_supported = [&](size_t bi) {
        std::vector<std::pair<int, int>> top;
        for (const auto& [a, n] : observed[bi]) top.emplace_back(n, a);
        std::sort(top.rbegin(), top.rend());
        if (top.size() < 2) return false;
        for (size_t i = 0; i < 2; ++i) {
            const auto it = confident[bi].find(top[i].second);
            if (it == confident[bi].end() || it->second < kUnionMinConfidentAlleleReads) return false;
        }
        return true;
    };

    // Clean heterozygotes and MSA-verified noisy heterozygotes, keeping the
    // alignment's own category bits.
    std::vector<size_t> selected;
    for (size_t bi = 0; bi < bam.candidates.size(); ++bi) {
        const CandidateVariant& c = bam.candidates[bi];
        const VariantCategory cat = c.counts.category;
        const bool clean = cat == VariantCategory::CleanHetSnp || cat == VariantCategory::CleanHetIndel;
        const bool msa = cat == VariantCategory::NoisyCandHet && c.msa_verified;
        if (!clean && !msa) continue;
        if ((c.lcd_var_i_to_cate & kCandGermlineVarCate) == 0) continue;
        if (c.key.pos < graph.ref_beg || c.key.pos > graph.ref_end) continue;
        if (!confidently_supported(bi)) continue;
        std::vector<size_t> rows = voting_rows_over(c.key.pos - 1, c.key.pos + std::max(1, c.key.ref_len));
        if (!rows.empty()) complement_rows.emplace(bi, std::move(rows));
        selected.push_back(bi);
    }
    if (selected.empty()) return 0;

    // One site per left-aligned allele: other placements become aliases whose
    // read calls fold into the canonical row.
    std::unordered_map<size_t, size_t> alias_of;
    {
        std::map<CandKey, size_t> canonical;
        std::vector<size_t> kept;
        for (const size_t bi : selected) {
            const auto [it, inserted] = canonical.emplace(left_aligned_indel_key(bam.candidates[bi], bam), bi);
            if (inserted) kept.push_back(bi);
            else alias_of.emplace(bi, it->second);
        }
        selected = std::move(kept);
    }

    // Merge in key order; every index-parallel array is rebuilt with the table.
    struct Slot { CandKey key; long graph_index; long bam_index; };
    std::vector<Slot> slots;
    slots.reserve(graph.candidates.size() + selected.size());
    for (size_t ci = 0; ci < graph.candidates.size(); ++ci)
        slots.push_back(Slot{cand_key_of(graph.candidates[ci]), static_cast<long>(ci), -1});
    for (const size_t bi : selected) {
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
            merged_ids.push_back(ci < graph_chunk.site_ids.size() ? graph_chunk.site_ids[ci] : std::string());
            merged_meta.push_back(had_meta ? graph_chunk.site_meta[ci] : GraphSiteMeta{});
            merged_orig.push_back(ci < graph_chunk.site_allele_orig_idx.size()
                                      ? graph_chunk.site_allele_orig_idx[ci] : std::vector<int>{});
            continue;
        }
        const size_t bi = static_cast<size_t>(slot.bam_index);
        bam_to_new.emplace(bi, ni);
        CandidateVariant cand = bam.candidates[bi];
        cand.key.tid = graph.region.tid;
        cand.bam_injected = true;
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
        if (complement_rows.count(bi) != 0) {
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
    for (const auto& [bi, rows] : complement_rows) {
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
    return selected.size();
}
// ── Global EM ───────────────────────────────────────────────────────────────

} // namespace pgphase_collect
