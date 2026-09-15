#include "gap_recovery.hpp"
#include "gap_evidence.hpp"

#include "collect_phase.hpp"
#include "collect_phase_noisy.hpp"

#include <algorithm>
#include <iostream>
#include <map>
#include <set>
#include <tuple>
#include <string_view>
#include <unordered_map>

namespace pgphase_collect {

static constexpr hts_pos_t kGapMsaWindow = 512;
static constexpr hts_pos_t kGapMsaOverlap = 64;

int select_graph_gap_bam_reads(PhasingChunk& proposal, const PhaseGap& gap,
                               const Options& opts) {
    int selected = 0;
    for (size_t ri = 0; ri < proposal.reads.size(); ++ri) {
        auto& read = proposal.reads[ri];
        auto& profile = proposal.read_var_profile[ri];
        hts_pos_t graph_beg = -1, graph_end = -1;
        for (size_t pi = 0; pi < profile.graph_alleles.size(); ++pi) {
            const int vi = profile.start_var_idx + static_cast<int>(pi);
            if (profile.graph_alleles[pi] < 0 || vi < 0 ||
                static_cast<size_t>(vi) >= proposal.candidates.size()) continue;
            const auto pos = proposal.candidates[vi].key.sort_pos();
            if (graph_beg < 0) graph_beg = pos;
            graph_end = pos;
        }
        if (read.is_skipped || !read.alignment || read.mapq < opts.min_mapq ||
            (read.alignment->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) ||
            graph_beg < 0 || graph_beg > gap.right_beg || graph_end < gap.left_end ||
            read.end < gap.left_end || read.beg > gap.right_beg) {
            read.is_skipped = true;
            continue;
        }
        ++selected;
        const bam1_t* bam = read.alignment.get();
        for (size_t pi = 0; pi < profile.alleles.size(); ++pi) {
            const int vi = profile.start_var_idx + static_cast<int>(pi);
            if (vi < 0 || static_cast<size_t>(vi) >= proposal.candidates.size()) continue;
            const auto& candidate = proposal.candidates[vi];
            if (!candidate.msa_verified && pi < profile.bam_alleles.size()) {
                profile.alleles[pi] = profile.bam_alleles[pi];
                if (pi < profile.alt_qi.size())
                    profile.alt_qi[pi] = pi < profile.bam_qi.size() ? profile.bam_qi[pi] : -1;
                continue;
            }
            if (candidate.key.sort_pos() < read.beg || candidate.key.sort_pos() > read.end) {
                profile.alleles[pi] = -1;
                continue;
            }
            if (candidate.key.type != VariantType::Snp || candidate.msa_verified) {
                if (pi < profile.alt_qi.size() &&
                    profile.alt_qi[pi] == kGraphConfirmedAltQi && !candidate.msa_verified)
                    profile.alleles[pi] = -1;
                continue;
            }
            // Read the nucleotide at this reference position from the BAM,
            // including when hybrid injection previously filled a missing allele.
            profile.alleles[pi] = -1;
            hts_pos_t ref = bam->core.pos + 1;
            int query = 0;
            const uint32_t* cigar = bam_get_cigar(bam);
            for (uint32_t ci = 0; ci < bam->core.n_cigar; ++ci) {
                const int op = bam_cigar_op(cigar[ci]);
                const int len = bam_cigar_oplen(cigar[ci]);
                const int consumes = bam_cigar_type(op);
                if ((consumes & 2) && candidate.key.pos >= ref && candidate.key.pos < ref + len) {
                    if (consumes & 1) {
                        const int qi = query + static_cast<int>(candidate.key.pos - ref);
                        if (bam_get_qual(bam)[qi] != 255 && bam_get_qual(bam)[qi] >= opts.min_bq) {
                            const char base = seq_nt16_str[bam_seqi(bam_get_seq(bam), qi)];
                            const char ref_base = candidate.ref_base < 4 ? "ACGT"[candidate.ref_base] : 'N';
                            if (base == ref_base) profile.alleles[pi] = 0;
                            else if (candidate.key.alt.size() == 1 && base == candidate.key.alt[0])
                                profile.alleles[pi] = 1;
                            if (pi < profile.alt_qi.size()) profile.alt_qi[pi] = qi;
                        }
                    }
                    break;
                }
                if (consumes & 1) query += len;
                if (consumes & 2) ref += len;
            }
        }
    }
    return selected;
}

std::vector<PhaseGap> find_phase_gaps(const std::vector<PhasingChunk>& chunks) {
    using BlockKey = std::pair<int, hts_pos_t>;
    std::set<BlockKey> supported;
    for (const auto& chunk : chunks) {
        for (size_t i = 0; i < chunk.haps.size(); ++i) {
            if (chunk.haps[i] != 0 && !chunk.reads[i].is_skipped && chunk.phase_sets[i] >= 0)
                supported.emplace(chunk.region.tid, chunk.phase_sets[i]);
        }
    }
    std::map<BlockKey, std::pair<hts_pos_t, hts_pos_t>> bounds;
    for (const auto& chunk : chunks) {
        for (const auto& v : chunk.candidates) {
            const BlockKey key{chunk.region.tid, v.phase_set};
            if (!supported.count(key) || v.hap_to_cons_alle[1] < 0 ||
                v.hap_to_cons_alle[2] < 0 || v.hap_to_cons_alle[1] == v.hap_to_cons_alle[2]) continue;
            const hts_pos_t pos = v.key.sort_pos();
            auto inserted = bounds.emplace(key, std::make_pair(pos, pos));
            auto& span = inserted.first->second;
            span.first = std::min(span.first, pos);
            span.second = std::max(span.second, pos);
        }
    }
    struct Block { int tid; hts_pos_t ps, beg, end; };
    std::vector<Block> blocks;
    for (const auto& [key, span] : bounds)
        blocks.push_back({key.first, key.second, span.first, span.second});
    std::sort(blocks.begin(), blocks.end(), [](const Block& a, const Block& b) {
        return std::tie(a.tid, a.beg, a.end, a.ps) < std::tie(b.tid, b.beg, b.end, b.ps);
    });
    std::vector<RegionChunk> coverage;
    for (const auto& chunk : chunks) coverage.push_back(chunk.region);
    std::sort(coverage.begin(), coverage.end(), [](const RegionChunk& a, const RegionChunk& b) {
        return std::tie(a.tid, a.beg) < std::tie(b.tid, b.beg);
    });
    std::vector<RegionChunk> components;
    for (const auto& region : coverage) {
        if (components.empty() || components.back().tid != region.tid ||
            region.beg > components.back().end + 1) components.push_back(region);
        else components.back().end = std::max(components.back().end, region.end);
    }
    std::vector<PhaseGap> gaps;
    if (blocks.empty()) return gaps;
    Block left = blocks.front();
    for (size_t i = 1; i < blocks.size(); ++i) {
        const Block& right = blocks[i];
        if (right.tid == left.tid && right.beg > left.end) {
            for (const auto& region : components) {
                if (region.tid == left.tid && region.beg <= left.end && region.end >= right.beg) {
                    gaps.push_back({left.tid, left.ps, right.ps, left.end, right.beg, region.beg, region.end});
                    break;
                }
            }
        }
        // Nested/overlapping blocks do not define an empty genomic gap.
        if (right.tid != left.tid || right.end > left.end) left = right;
    }
    return gaps;
}

static void orient_candidate(CandidateVariant& v, hts_pos_t ps, bool flip) {
    v.phase_set = ps;
    if (!flip) return;
    std::swap(v.hap_to_cons_alle[1], v.hap_to_cons_alle[2]);
    std::swap(v.hap_to_alle_profile[1], v.hap_to_alle_profile[2]);
    if (v.hap_alt == 1 || v.hap_alt == 2) v.hap_alt = 3 - v.hap_alt;
    if (v.hap_ref == 1 || v.hap_ref == 2) v.hap_ref = 3 - v.hap_ref;
}

size_t GapReadIndex::Hash::operator()(const Key& key) const {
    return std::hash<std::string_view>{}(key.second) ^ std::hash<int>{}(key.first);
}

GapReadIndex::GapReadIndex(const std::vector<PhasingChunk>& chunks) {
    for (size_t ci = 0; ci < chunks.size(); ++ci) {
        for (size_t ri = 0; ri < chunks[ci].reads.size(); ++ri) {
            const auto& read = chunks[ci].reads[ri];
            const Key key{read.input_index, read.qname};
            reads[key].emplace_back(ci, ri);
            if (!read.is_skipped && chunks[ci].haps[ri] != 0 &&
                chunks[ci].phase_sets[ri] >= 0) {
                assignments.emplace(
                    key, std::make_pair(chunks[ci].haps[ri],
                                        chunks[ci].phase_sets[ri]));
                established_phase_sets.insert(chunks[ci].phase_sets[ri]);
            }
        }
        // filter_hybrid_reads_by_margin / filter_hybrid_small_phase_sets zero
        // a read's haps/phase_sets when its block loses read support, but
        // never touch CandidateVariant::phase_set -- an orphaned block can
        // leave candidate rows still labelled with its id and the main
        // solve's orientation, with zero committed reads to show for it
        // (exactly why find_phase_gaps has no read support for that id
        // either, and the region becomes a gap). That id is just as real a
        // collision risk for emit_independent_gap_block as a read-supported
        // one: a gap's local re-solve reproducing it would insert new rows
        // under the proposal's own polarity while the orphaned rows keep the
        // main solve's, silently splitting one phase set's orientation.
        for (const auto& candidate : chunks[ci].candidates)
            if (candidate.phase_set > 0) established_phase_sets.insert(candidate.phase_set);
    }
}

GapStitchResult stitch_gap_proposal(std::vector<PhasingChunk>& chunks,
                                    const PhasingChunk& proposal,
                                    const PhaseGap& gap, const Options& opts,
                                    const GapReadIndex* read_index,
                                    bool defer_phase_set_merge,
                                    bool orientation_only,
                                    std::vector<GapLinkEvidence>* evidence) {
    using ReadKey = GapReadIndex::Key;
    const auto local_index = read_index == nullptr ? std::make_unique<GapReadIndex>(chunks) : nullptr;
    const auto& index = read_index == nullptr ? *local_index : *read_index;
    auto first_assignment = [&](const ReadKey& key, const auto& locations) {
        if (read_index != nullptr) {
            const auto frozen = read_index->assignments.find(key);
            if (frozen != read_index->assignments.end()) return frozen->second;
            return std::make_pair(0, static_cast<hts_pos_t>(-1));
        }
        for (const auto& [ci, ri] : locations) {
            const auto& chunk = chunks[ci];
            if (chunk.region.tid != gap.tid || chunk.reads[ri].is_skipped ||
                chunk.haps[ri] == 0 || chunk.phase_sets[ri] < 0) continue;
            return std::make_pair(chunk.haps[ri], chunk.phase_sets[ri]);
        }
        return std::make_pair(0, static_cast<hts_pos_t>(-1));
    };
    // Fallback for a read with no committed pre-recovery hap/PS at all: derive
    // an implied flank side from its own allele agreement with that flank's
    // already-resolved sites (mirrors check_agree_alleles, generalized from a
    // variant-pair edge to a read-vs-block edge). Only consulted when the
    // committed lookup above finds nothing.
    //
    // A read can appear in more than one overlapping tile chunk
    // (GapReadIndex::reads holds every location). CandidateVariant::
    // hap_to_cons_alle for a gap.left_ps/right_ps site is only guaranteed to
    // be in this batch's canonical hap1/hap2 orientation for chunks this
    // batch's orient_candidate calls have actually touched; an untouched
    // duplicate tile copy of the same physical site can still carry its own,
    // independently k-means-derived (and therefore possibly oppositely
    // labelled) sense for the same phase-set id. Trusting whichever location
    // is seen first measurably corrupted phase sets on chr20 (a single PS
    // ending up with contradictory HP calls across its reads -- see
    // CHECKPOINT.md, 2026-09-15). Require every location that resolves a
    // given phase set to agree; any disagreement poisons that phase set for
    // this read rather than picking a side.
    auto implied_assignment = [&](const auto& locations) {
        std::map<hts_pos_t, int> per_ps; // ps -> implied hp (1/2), or -1 once conflicting
        for (const auto& [ci, ri] : locations) {
            const auto& chunk = chunks[ci];
            if (chunk.region.tid != gap.tid || ri >= chunk.reads.size() ||
                chunk.reads[ri].is_skipped || ri >= chunk.read_var_profile.size()) continue;
            const auto& prof = chunk.read_var_profile[ri];
            if (prof.start_var_idx < 0 || prof.end_var_idx < prof.start_var_idx) continue;
            std::map<hts_pos_t, int> local; // this one location's own per-ps implied hp
            for (int vi = prof.start_var_idx;
                 vi <= prof.end_var_idx && static_cast<size_t>(vi) < chunk.candidates.size(); ++vi) {
                const auto& site = chunk.candidates[vi];
                if (site.phase_set != gap.left_ps && site.phase_set != gap.right_ps) continue;
                if (site.hap_to_cons_alle[1] < 0 || site.hap_to_cons_alle[2] < 0 ||
                    site.hap_to_cons_alle[1] == site.hap_to_cons_alle[2]) continue;
                const size_t allele_idx = static_cast<size_t>(vi - prof.start_var_idx);
                if (allele_idx >= prof.alleles.size()) continue;
                const int allele = prof.alleles[allele_idx];
                if (allele < 0) continue;
                int h = 0;
                if (site.hap_to_cons_alle[1] == allele) h = 1;
                else if (site.hap_to_cons_alle[2] == allele) h = 2;
                if (h == 0) continue;
                const auto it = local.find(site.phase_set);
                if (it == local.end()) local.emplace(site.phase_set, h);
                else if (it->second != h) it->second = -1;
            }
            for (const auto& [ps, h] : local) {
                const auto it = per_ps.find(ps);
                if (it == per_ps.end()) per_ps.emplace(ps, h);
                else if (h < 0 || it->second != h) it->second = -1;
            }
        }
        for (const hts_pos_t ps : {gap.left_ps, gap.right_ps}) {
            const auto it = per_ps.find(ps);
            if (it != per_ps.end() && it->second > 0) return std::make_pair(it->second, ps);
        }
        return std::make_pair(0, static_cast<hts_pos_t>(-1));
    };
    std::set<hts_pos_t> supported_phase_sets;
    for (const auto& [key, locations] : index.reads) {
        const auto assignment = first_assignment(key, locations);
        if (assignment.first != 0) supported_phase_sets.insert(assignment.second);
    }
    // Only proposal reads are queried below. Retain the same first valid
    // assignment in chunk order without rebuilding a chromosome-sized name map.
    // Vote tallying (below) and therefore the join/orientation decision itself
    // use ONLY committed pre-recovery assignments -- byte-identical to the
    // decision the default (gap_link_by_alleles off) path makes. Allele-
    // derived agreement fed the same vote directly in an earlier version of
    // this function and measurably corrupted joins on chr20 (see
    // CHECKPOINT.md, 2026-09-15); it is now consulted only after a join is
    // already accepted on committed evidence alone, to additively attach more
    // reads to it (below), the same purely-additive contract
    // emit_independent_gap_block uses.
    std::unordered_map<ReadKey, std::pair<int, hts_pos_t>, GapReadIndex::Hash> original;
    for (const auto& read : proposal.reads) {
        const ReadKey key{read.input_index, read.qname};
        const auto found = index.reads.find(key);
        if (found == index.reads.end()) continue;
        const auto assignment = first_assignment(key, found->second);
        if (assignment.first != 0) original.emplace(key, assignment);
    }
    std::map<ReadKey, size_t> proposal_reads;
    std::map<ReadKey, size_t> observation_reads;
    std::map<hts_pos_t, std::array<std::array<int, 4>, 2>> votes;
    for (size_t i = 0; i < proposal.reads.size(); ++i) {
        const auto& read = proposal.reads[i];
        if (read.is_skipped) continue;
        const ReadKey key{read.input_index, read.qname};
        observation_reads.emplace(key, i);
        if (proposal.haps[i] == 0 || proposal.phase_sets[i] < 0) continue;
        if (!proposal_reads.emplace(key, i).second) continue;
        const auto it = original.find(key);
        if (it == original.end()) continue;
        const auto [hp, ps] = it->second;
        const int side = ps == gap.left_ps ? 0 : ps == gap.right_ps ? 1 : -1;
        if (side < 0) continue;
        ++votes[proposal.phase_sets[i]][side][(hp - 1) * 2 + proposal.haps[i] - 1];
    }
    struct Link { hts_pos_t ps = -1; bool flip = false; int support = -1; };
    std::array<Link, 2> links;
    std::array<Link, 2> bridge;
    int bridge_support = -1;
    if (evidence != nullptr) evidence->clear();
    for (const auto& [ps, sides] : votes) {
        if (evidence != nullptr) evidence->push_back({ps, sides});
        std::array<Link, 2> pair;
        for (int side = 0; side < 2; ++side) {
            bool flip = false;
            if (opts.verbose >= 2) {
                std::cerr << "GapLinkVotes\t" << gap.left_end << '\t' << gap.right_beg
                          << '\t' << ps << '\t' << side;
                for (const int count : sides[side]) std::cerr << '\t' << count;
                std::cerr << '\n';
            }
            if (!select_stitch_orientation(sides[side], &opts, flip)) continue;
            if (opts.gap_hp_link_beg >= 0) {
                // Each original haplotype must independently favor the same
                // orientation before committing a homopolymer rescue.
                const auto& v = sides[side];
                const int first = flip ? v[1] - v[0] : v[0] - v[1];
                const int second = flip ? v[2] - v[3] : v[3] - v[2];
                if (std::min(first, second) < opts.min_block_link_reads) continue;
            }
            const int support = flip ? sides[side][1] + sides[side][2]
                                     : sides[side][0] + sides[side][3];
            pair[side] = {ps, flip, support};
            if (support > links[side].support) links[side] = pair[side];
        }
        const int support = std::min(pair[0].support, pair[1].support);
        if (support > bridge_support) {
            bridge_support = support;
            bridge = pair;
        }
    }
    // Prefer a proposal block that actually connects both flanks over two
    // independently stronger, disconnected one-sided proposal blocks.
    if (bridge_support >= 0) links = bridge;
    GapStitchResult result;
    result.left_linked = links[0].ps >= 0;
    result.right_linked = links[1].ps >= 0;
    result.joined = result.left_linked && result.right_linked && links[0].ps == links[1].ps;
    result.right_flip = result.joined && links[0].flip != links[1].flip;
    if (orientation_only) return result;
    // Failed last-resort trials must not extend or modify either trusted flank.
    if (opts.gap_hp_link_beg >= 0 && !result.joined) return result;
    std::map<hts_pos_t, std::pair<hts_pos_t, bool>> accepted;
    if (result.left_linked) accepted[links[0].ps] = {gap.left_ps, links[0].flip};
    if (result.right_linked && !result.joined)
        accepted[links[1].ps] = {gap.right_ps, links[1].flip};
    if (accepted.empty()) return result;

    const bool right_flip = result.right_flip;
    std::set<ReadKey> added;
    // Flank phase-set id(s) this join actually resolved to (already the
    // decision above; allele agreement below only attaches more reads to it,
    // never influences which side wins or whether it flips).
    std::set<hts_pos_t> accepted_targets;
    if (opts.gap_link_by_alleles)
        for (const auto& [local_ps, tgt] : accepted) accepted_targets.insert(tgt.first);
    for (auto& chunk : chunks) {
        if (chunk.region.tid != gap.tid) continue;
        if (result.joined && !defer_phase_set_merge) {
            for (size_t i = 0; i < chunk.haps.size(); ++i) {
                if (chunk.phase_sets[i] != gap.right_ps) continue;
                chunk.phase_sets[i] = gap.left_ps;
                if (right_flip && chunk.haps[i] != 0) chunk.haps[i] = 3 - chunk.haps[i];
            }
            for (auto& v : chunk.candidates)
                if (v.phase_set == gap.right_ps) orient_candidate(v, gap.left_ps, right_flip);
        }
        for (size_t i = 0; i < chunk.reads.size(); ++i) {
            auto& read = chunk.reads[i];
            if (read.is_skipped || (chunk.haps[i] != 0 && chunk.phase_sets[i] >= 0)) continue;
            const ReadKey key{read.input_index, read.qname};
            // Another overlapping chunk may already own a trusted assignment.
            if (original.count(key)) continue;
            const auto found = proposal_reads.find(key);
            if (found == proposal_reads.end()) continue;
            const size_t pi = found->second;
            const auto link = accepted.find(proposal.phase_sets[pi]);
            if (link == accepted.end()) continue;
            chunk.haps[i] = link->second.second ? 3 - proposal.haps[pi] : proposal.haps[pi];
            chunk.phase_sets[i] = link->second.first;
            added.insert(key);
            const auto& src = proposal.reads[pi];
            read.n_clean_agree_snps = src.n_clean_agree_snps;
            read.n_clean_conflict_snps = src.n_clean_conflict_snps;
            read.n_bridge_agree_snps = src.n_bridge_agree_snps;
            read.n_bridge_conflict_snps = src.n_bridge_conflict_snps;
            read.hap_score_margin = src.hap_score_margin;
            read.n_vars_scored = src.n_vars_scored;
        }
        if (opts.gap_link_by_alleles) {
            // Purely additive: a read with no committed pre-recovery
            // hap/PS at all, and no membership in the gap's own winning
            // local bucket above, can still directly agree by allele with
            // whichever flank this join already resolved to (candidates in
            // this chunk are already in that flank's canonical orientation,
            // reoriented above if this chunk held the merged-away side) --
            // attach it with no further flip. Never touches a read the loop
            // above already claimed, and never affects which side won or its
            // orientation.
            //
            // Deliberately does NOT skip on `added` (unlike a read that is
            // already gated by this chunk's own `chunk.haps[i]!=0` check):
            // `implied_assignment` reads THIS chunk's own live candidates,
            // which only gain the gap's interior sites via this same per-
            // chunk loop's own merge_var_profile call below -- a read's
            // earlier (output-owning) tile copy can still be missing those
            // sites at the point its own attach attempt runs here, while a
            // later, overlapping tile copy already has them and succeeds.
            // Gating on `added` would then leave exactly the output-owning
            // copy unphased while a non-owning copy silently absorbed the
            // assignment. Letting every qualifying copy attempt
            // independently is redundant when they agree (same shared
            // orient_candidate output feeds every copy) but never
            // corrupting, and `added` (a set) still dedupes the count below.
            for (size_t i = 0; i < chunk.reads.size(); ++i) {
                auto& read = chunk.reads[i];
                if (read.is_skipped || (chunk.haps[i] != 0 && chunk.phase_sets[i] >= 0)) continue;
                const ReadKey key{read.input_index, read.qname};
                if (original.count(key)) continue;
                const auto found = index.reads.find(key);
                if (found == index.reads.end()) continue;
                const auto implied = implied_assignment(found->second);
                if (implied.first == 0 || !accepted_targets.count(implied.second)) continue;
                chunk.haps[i] = implied.first;
                chunk.phase_sets[i] = implied.second;
                added.insert(key);
            }
        }
        std::vector<CandidateVariant> sites;
        std::vector<VariantCategory> categories;
        std::vector<size_t> indices;
        for (size_t pi = 0; pi < proposal.candidates.size(); ++pi) {
            const auto& v = proposal.candidates[pi];
            const auto link = accepted.find(v.phase_set);
            if (link == accepted.end() || v.key.sort_pos() < chunk.region.beg ||
                v.key.sort_pos() > chunk.region.end ||
                (v.lcd_var_i_to_cate & kCandGermlineVarCate) == 0) continue;
            sites.push_back(v);
            orient_candidate(sites.back(), link->second.first, link->second.second);
            categories.push_back(v.counts.category);
            indices.push_back(pi);
        }
        if (sites.empty()) continue;
        std::vector<ReadVariantProfile> profiles(chunk.reads.size());
        for (size_t i = 0; i < chunk.reads.size(); ++i) {
            const ReadKey key{chunk.reads[i].input_index, chunk.reads[i].qname};
            const auto found = observation_reads.find(key);
            if (found == observation_reads.end()) continue;
            const auto& source = proposal.read_var_profile[found->second];
            auto& dest = profiles[i];
            dest.start_var_idx = 0;
            dest.end_var_idx = static_cast<int>(sites.size()) - 1;
            dest.alleles.assign(sites.size(), -1);
            dest.alt_qi.assign(sites.size(), -1);
            dest.graph_alleles.assign(sites.size(), -1);
            for (size_t vi = 0; vi < indices.size(); ++vi) {
                const int pi = static_cast<int>(indices[vi]);
                if (pi < source.start_var_idx || pi > source.end_var_idx) continue;
                dest.alleles[vi] = source.alleles[pi - source.start_var_idx];
                if (!source.alt_qi.empty()) dest.alt_qi[vi] = source.alt_qi[pi - source.start_var_idx];
                if (static_cast<size_t>(pi - source.start_var_idx) <
                    source.graph_alleles.size())
                    dest.graph_alleles[vi] = source.graph_alleles[
                        static_cast<size_t>(pi - source.start_var_idx)];
            }
        }
        // Keep trusted core orientations. Previously unsupported exact matches
        // may adopt the proposal's phase; validated repeats also need the MSA
        // category and MSA observations, not their excluded first-pass category.
        VariantKeySet replace_sites;
        std::vector<std::pair<VariantKey, uint32_t>> inserted_flags;
        for (const auto& site : sites) {
            const auto old = std::lower_bound(chunk.candidates.begin(), chunk.candidates.end(), site.key,
                [](const CandidateVariant& v, const VariantKey& key) {
                    return exact_comp_var_site(&v.key, &key) < 0;
                });
            if (old == chunk.candidates.end() || exact_comp_var_site(&old->key, &site.key) != 0) {
                inserted_flags.emplace_back(site.key, site.lcd_var_i_to_cate);
                continue;
            }
            if (supported_phase_sets.count(old->phase_set)) continue;
            replace_sites.insert(site.key);
            inserted_flags.emplace_back(site.key, site.lcd_var_i_to_cate);
        }
        merge_var_profile(chunk, sites, categories, profiles, nullptr, false, false, &replace_sites);
        // The generic MSA merger derives flags from categories. Clean proposal
        // rows may carry a stricter non-anchor mask that must survive insertion.
        for (const auto& [key, flags] : inserted_flags) {
            const auto inserted = std::lower_bound(chunk.candidates.begin(), chunk.candidates.end(), key,
                [](const CandidateVariant& v, const VariantKey& key) {
                    return exact_comp_var_site(&v.key, &key) < 0;
                });
            inserted->lcd_var_i_to_cate = flags;
        }
    }
    result.reads_added = static_cast<int>(added.size());
    return result;
}

int apply_gap_phase_edges(std::vector<PhasingChunk>& chunks,
                          const std::vector<GapPhaseEdge>& edges) {
    std::unordered_map<hts_pos_t, hts_pos_t> parent;
    std::unordered_map<hts_pos_t, bool> parity;
    auto ensure = [&](hts_pos_t ps) {
        parent.emplace(ps, ps);
        parity.emplace(ps, false);
    };
    auto find = [&](hts_pos_t ps) {
        hts_pos_t root = ps;
        bool flip = false;
        while (parent[root] != root) {
            flip ^= parity[root];
            root = parent[root];
        }
        hts_pos_t node = ps;
        bool prefix = false;
        while (parent[node] != node) {
            const hts_pos_t next = parent[node];
            const bool edge_flip = parity[node];
            parent[node] = root;
            parity[node] = flip ^ prefix;
            prefix ^= edge_flip;
            node = next;
        }
        return std::make_pair(root, flip);
    };

    int conflicts = 0;
    for (const GapPhaseEdge& edge : edges) {
        ensure(edge.left_ps);
        ensure(edge.right_ps);
        const auto left = find(edge.left_ps);
        const auto right = find(edge.right_ps);
        if (left.first == right.first) {
            if ((left.second ^ right.second) != edge.right_flip) ++conflicts;
            continue;
        }
        parent[right.first] = left.first;
        parity[right.first] = left.second ^ right.second ^ edge.right_flip;
    }

    for (PhasingChunk& chunk : chunks) {
        for (CandidateVariant& candidate : chunk.candidates) {
            if (candidate.phase_set < 0 || !parent.count(candidate.phase_set)) continue;
            const auto resolved = find(candidate.phase_set);
            orient_candidate(candidate, resolved.first, resolved.second);
        }
        for (size_t i = 0; i < chunk.phase_sets.size(); ++i) {
            const hts_pos_t ps = chunk.phase_sets[i];
            if (ps < 0 || !parent.count(ps)) continue;
            const auto resolved = find(ps);
            chunk.phase_sets[i] = resolved.first;
            if (resolved.second && chunk.haps[i] != 0)
                chunk.haps[i] = 3 - chunk.haps[i];
        }
    }
    return conflicts;
}

int emit_independent_gap_block(std::vector<PhasingChunk>& chunks,
                               const PhasingChunk& proposal,
                               const PhaseGap& gap, const Options& opts,
                               const GapReadIndex& read_index,
                               int min_reads,
                               hts_pos_t* emitted_ps,
                               std::set<hts_pos_t>* emitted_this_round) {
    // Individual read confidence (min-read-margin) was already applied
    // upstream by filter_hybrid_reads_by_margin before this is called, so
    // proposal.haps[pi] != 0 here already reflects that gate.
    if (emitted_ps != nullptr) *emitted_ps = -1;
    if (min_reads <= 0) return 0;
    using ReadKey = GapReadIndex::Key;

    // A gap-only read is one with no pre-recovery haplotype/phase-set
    // assignment at all -- neither flank, nor any other original block. Only
    // such reads are eligible: anything with an existing assignment is left
    // strictly to stitch_gap_proposal's flank-vote logic, never to this path.
    //
    // `read_index.assignments` is a frozen pre-recovery snapshot, so by the
    // time this gap's turn comes up in the strictly sequential batch loop
    // that calls this function, a read with no assignment in the snapshot
    // may already have been claimed by an earlier-processed bridge join or
    // independent block this same round (an adjacent gap's window can
    // overlap this one). Tried additionally excluding those via a live
    // `chunks` lookup so the chosen group reflects only what is still
    // genuinely available; measured on chr20 it changed which of two
    // competing local groups wins often enough to make 116 reads that were
    // fine before (114 concordant) no longer independent-block members at
    // all, for a wash in net reads gained -- reverted. The existing apply
    // loop below already never overwrites a read with any live assignment,
    // which is the actual safety property that matters; which group wins
    // the largest-frozen-count tie only needs to be deterministic, not
    // live-accurate.
    std::map<hts_pos_t, std::vector<size_t>> candidates_by_ps;
    int locally_phased = 0, has_original = 0;
    for (size_t pi = 0; pi < proposal.reads.size(); ++pi) {
        const auto& read = proposal.reads[pi];
        if (read.is_skipped || proposal.haps[pi] == 0 || proposal.phase_sets[pi] < 0) continue;
        ++locally_phased;
        const ReadKey key{read.input_index, read.qname};
        if (read_index.assignments.count(key)) { ++has_original; continue; }  // has an original assignment; not gap-only
        // The local re-solve's phase-set id is a genome position, same as a
        // real established block's -- a wide-window local solve can, by
        // coincidence, anchor on the same first het site as an adjacent
        // established block and reproduce its exact id. Never treat such a
        // group as "new": emitting into it would silently pollute a real
        // block with reads whose orientation was never vote-checked against
        // it (see established_phase_sets).
        if (read_index.established_phase_sets.count(proposal.phase_sets[pi])) { ++has_original; continue; }
        // `read_index` is a whole-batch snapshot frozen before this
        // sequential loop starts, so it cannot know about an id another gap
        // already emitted into earlier in this same loop. Adjacent gaps'
        // windows overlap by construction (kGapRecoveryFlank on both sides),
        // so two gaps can independently reconstruct the same leftmost-het id
        // from overlapping gap-only read pools; without this check they could
        // each emit into "chosen_ps" with an independently-derived, possibly
        // opposite orientation.
        if (emitted_this_round != nullptr && emitted_this_round->count(proposal.phase_sets[pi])) {
            ++has_original;
            continue;
        }
        candidates_by_ps[proposal.phase_sets[pi]].push_back(pi);
    }

    // Pick the largest qualifying group. A gap rarely has more than one real
    // internal split; taking the single best avoids emitting several tiny,
    // low-confidence fragments from noise in the same window.
    hts_pos_t chosen_ps = -1;
    size_t chosen_size = 0;
    for (const auto& [ps, members] : candidates_by_ps) {
        if (members.size() > chosen_size) {
            chosen_ps = ps;
            chosen_size = members.size();
        }
    }
    if (opts.verbose >= 2) {
        std::cerr << "GapIndependentBlock\t" << gap.left_end << '\t' << gap.right_beg
                  << "\tlocally_phased=" << locally_phased << "\thas_original=" << has_original
                  << "\tgap_only_ps_groups=" << candidates_by_ps.size()
                  << "\tchosen_ps=" << chosen_ps << "\tchosen_size=" << chosen_size << '\n';
    }
    if (chosen_ps < 0 || static_cast<int>(chosen_size) < min_reads) return 0;
    if (emitted_ps != nullptr) *emitted_ps = chosen_ps;
    if (emitted_this_round != nullptr) emitted_this_round->insert(chosen_ps);

    std::map<ReadKey, size_t> proposal_reads;
    std::map<ReadKey, size_t> observation_reads;
    for (size_t pi = 0; pi < proposal.reads.size(); ++pi) {
        const auto& read = proposal.reads[pi];
        if (read.is_skipped) continue;
        const ReadKey key{read.input_index, read.qname};
        observation_reads.emplace(key, pi);
        if (proposal.haps[pi] != 0 && proposal.phase_sets[pi] == chosen_ps)
            proposal_reads.emplace(key, pi);
    }

    int applied = 0;
    for (auto& chunk : chunks) {
        if (chunk.region.tid != gap.tid) continue;
        for (size_t i = 0; i < chunk.reads.size(); ++i) {
            auto& read = chunk.reads[i];
            // Never overwrite a read that already carries any assignment,
            // from this gap's own bridge attempt or any other source.
            if (read.is_skipped || (chunk.haps[i] != 0 && chunk.phase_sets[i] >= 0)) continue;
            const ReadKey key{read.input_index, read.qname};
            const auto found = proposal_reads.find(key);
            if (found == proposal_reads.end()) continue;
            const size_t pi = found->second;
            chunk.haps[i] = proposal.haps[pi];
            chunk.phase_sets[i] = proposal.phase_sets[pi];
            ++applied;
            const auto& src = proposal.reads[pi];
            read.n_clean_agree_snps = src.n_clean_agree_snps;
            read.n_clean_conflict_snps = src.n_clean_conflict_snps;
            read.n_bridge_agree_snps = src.n_bridge_agree_snps;
            read.n_bridge_conflict_snps = src.n_bridge_conflict_snps;
            read.hap_score_margin = src.hap_score_margin;
            read.n_vars_scored = src.n_vars_scored;
        }

        std::vector<CandidateVariant> sites;
        std::vector<VariantCategory> categories;
        std::vector<size_t> indices;
        for (size_t pi = 0; pi < proposal.candidates.size(); ++pi) {
            const auto& v = proposal.candidates[pi];
            if (v.phase_set != chosen_ps || v.key.sort_pos() < chunk.region.beg ||
                v.key.sort_pos() > chunk.region.end ||
                (v.lcd_var_i_to_cate & kCandGermlineVarCate) == 0) continue;
            sites.push_back(v);  // kept as-is: no flank to orient against.
            categories.push_back(v.counts.category);
            indices.push_back(pi);
        }
        if (sites.empty()) continue;

        std::vector<ReadVariantProfile> profiles(chunk.reads.size());
        for (size_t i = 0; i < chunk.reads.size(); ++i) {
            const ReadKey key{chunk.reads[i].input_index, chunk.reads[i].qname};
            const auto found = observation_reads.find(key);
            if (found == observation_reads.end()) continue;
            const auto& source = proposal.read_var_profile[found->second];
            auto& dest = profiles[i];
            dest.start_var_idx = 0;
            dest.end_var_idx = static_cast<int>(sites.size()) - 1;
            dest.alleles.assign(sites.size(), -1);
            dest.alt_qi.assign(sites.size(), -1);
            dest.graph_alleles.assign(sites.size(), -1);
            for (size_t vi = 0; vi < indices.size(); ++vi) {
                const int pi = static_cast<int>(indices[vi]);
                if (pi < source.start_var_idx || pi > source.end_var_idx) continue;
                dest.alleles[vi] = source.alleles[pi - source.start_var_idx];
                if (!source.alt_qi.empty()) dest.alt_qi[vi] = source.alt_qi[pi - source.start_var_idx];
                if (static_cast<size_t>(pi - source.start_var_idx) < source.graph_alleles.size())
                    dest.graph_alleles[vi] = source.graph_alleles[
                        static_cast<size_t>(pi - source.start_var_idx)];
            }
        }

        VariantKeySet replace_sites;
        std::vector<std::pair<VariantKey, uint32_t>> inserted_flags;
        for (const auto& site : sites) {
            const auto old = std::lower_bound(chunk.candidates.begin(), chunk.candidates.end(), site.key,
                [](const CandidateVariant& v, const VariantKey& key) {
                    return exact_comp_var_site(&v.key, &key) < 0;
                });
            inserted_flags.emplace_back(site.key, site.lcd_var_i_to_cate);
            // A brand-new independent block never has standing to replace an
            // existing row -- unlike a confirmed flank bridge, nothing here
            // has been validated against a trusted anchor.
            if (old != chunk.candidates.end() && exact_comp_var_site(&old->key, &site.key) == 0)
                inserted_flags.pop_back();
        }
        merge_var_profile(chunk, sites, categories, profiles, nullptr, false, false, &replace_sites);
        for (const auto& [key, flags] : inserted_flags) {
            const auto inserted = std::lower_bound(chunk.candidates.begin(), chunk.candidates.end(), key,
                [](const CandidateVariant& v, const VariantKey& key) {
                    return exact_comp_var_site(&v.key, &key) < 0;
                });
            if (inserted != chunk.candidates.end() && exact_comp_var_site(&inserted->key, &key) == 0)
                inserted->lcd_var_i_to_cate = flags;
        }
    }
    if (opts.verbose >= 2) {
        std::cerr << "GapIndependentBlockApplied\t" << gap.left_end << '\t' << gap.right_beg
                  << "\tchosen_size=" << chosen_size << "\tapplied=" << applied << '\n';
    }
    return applied;
}

int anchor_orphan_msa_sites(std::vector<PhasingChunk>& chunks, const Options& opts,
                            const std::vector<PhaseGap>* gaps) {
    std::set<std::pair<int, hts_pos_t>> supported;
    for (const auto& chunk : chunks)
        for (size_t ri = 0; ri < chunk.reads.size(); ++ri)
            if (!chunk.reads[ri].is_skipped && chunk.haps[ri] > 0 && chunk.phase_sets[ri] >= 0)
                supported.emplace(chunk.region.tid, chunk.phase_sets[ri]);
    int anchored = 0;
    int64_t* overlaps = nullptr;
    int64_t capacity = 0;
    for (auto& chunk : chunks) {
        if (!chunk.read_var_cr) continue;
        for (size_t vi = 0; vi < chunk.candidates.size(); ++vi) {
            auto& var = chunk.candidates[vi];
            if (gaps && std::none_of(gaps->begin(), gaps->end(), [&](const PhaseGap& gap) {
                    return gap_owns_variant(gap, var.key);
                })) continue;
            if (!var.msa_verified || var.phase_set < 0 ||
                var.hap_to_cons_alle[1] < 0 || var.hap_to_cons_alle[2] < 0 ||
                var.hap_to_cons_alle[1] == var.hap_to_cons_alle[2] ||
                supported.count({chunk.region.tid, var.phase_set})) continue;
            std::map<hts_pos_t, std::array<int, 4>> votes;
            std::set<std::pair<int, std::string_view>> seen;
            const auto n = cr_overlap(chunk.read_var_cr.get(), "cr", vi, vi + 1, &overlaps, &capacity);
            for (int64_t oi = 0; oi < n; ++oi) {
                const auto ri = static_cast<size_t>(cr_label(chunk.read_var_cr.get(), overlaps[oi]));
                const auto& read = chunk.reads[ri];
                const int hap = chunk.haps[ri];
                if (read.is_skipped || hap < 1 || hap > 2 || chunk.phase_sets[ri] < 0) continue;
                const auto& profile = chunk.read_var_profile[ri];
                if (profile.start_var_idx < 0 || vi < static_cast<size_t>(profile.start_var_idx) ||
                    vi > static_cast<size_t>(profile.end_var_idx)) continue;
                const int allele = profile.alleles[vi - profile.start_var_idx];
                const int side = allele == var.hap_to_cons_alle[1] ? 0 :
                                 allele == var.hap_to_cons_alle[2] ? 1 : -1;
                if (side < 0 || !seen.emplace(read.input_index, read.qname).second) continue;
                ++votes[chunk.phase_sets[ri]][2 * (hap - 1) + side];
            }
            hts_pos_t best_ps = -1;
            int best_support = -1;
            bool best_flip = false;
            for (const auto& [ps, v] : votes) {
                bool flip = false;
                if (!select_stitch_orientation(v, &opts, flip)) continue;
                const int first = flip ? v[1] - v[0] : v[0] - v[1];
                const int second = flip ? v[2] - v[3] : v[3] - v[2];
                const int support = std::min(first, second);
                if (support < opts.min_block_link_reads || support <= best_support) continue;
                best_ps = ps;
                best_support = support;
                best_flip = flip;
            }
            if (best_ps < 0) continue;
            orient_candidate(var, best_ps, best_flip);
            ++anchored;
        }
    }
    free(overlaps);
    return anchored;
}

void prepare_gap_msa_regions(PhasingChunk& chunk, hts_pos_t beg, hts_pos_t end) {
    beg = std::max(beg, chunk.ref_beg);
    end = std::min(end, chunk.ref_end);
    std::vector<Interval> regions;
    for (const auto& region : chunk.noisy_regions)
        if (region.end >= beg && region.beg <= end) regions.push_back(region);
    std::sort(regions.begin(), regions.end(), [](const Interval& a, const Interval& b) {
        return a.beg < b.beg;
    });
    // Tile clean stretches too: a phasing break need not trigger noise detection.
    auto tile = [&](hts_pos_t start, hts_pos_t stop) {
        for (hts_pos_t pos = start; pos <= stop; pos += kGapMsaWindow - kGapMsaOverlap) {
            const hts_pos_t last = std::min(stop, pos + kGapMsaWindow - 1);
            chunk.noisy_regions.push_back({pos, last, 0});
            if (last == stop) break;
        }
    };
    chunk.noisy_regions = regions;
    hts_pos_t next = beg;
    for (const auto& region : regions) {
        if (region.beg > next) tile(next, region.beg - 1);
        next = std::max(next, region.end + 1);
    }
    if (next <= end) tile(next, end);
}

void run_gap_msa_tier(PhasingChunk& chunk, const Options& opts,
                      hts_pos_t beg, hts_pos_t end, bool snp_only) {
    Options msa_opts = opts;
    msa_opts.private_msa = true;
    msa_opts.private_msa_admit_all_in_region = true;
    msa_opts.skip_noisy_kmeans = false;
    const auto sorted = sort_noisy_regs(chunk);
    std::vector<bool> done(chunk.noisy_regions.size(), false);
    while (true) {
        bool new_sites = false;
        for (const int index : sorted) {
            const auto& reg = chunk.noisy_regions[index];
            if (done[index] || reg.end < beg || reg.beg > end) continue;
            const int admitted = collect_noisy_vars1(chunk, msa_opts, index, nullptr, snp_only);
            if (admitted > 0) {
                const hts_pos_t backfill_beg = msa_opts.gap_recovery_beg >= 0
                    ? std::max(reg.beg, msa_opts.gap_recovery_beg) : reg.beg;
                const hts_pos_t backfill_end = msa_opts.gap_recovery_end >= 0
                    ? std::min(reg.end, msa_opts.gap_recovery_end) : reg.end;
                if (backfill_beg <= backfill_end)
                    backfill_msa_observations(chunk, msa_opts, backfill_beg, backfill_end);
            }
            if (admitted >= 0) done[index] = true;
            new_sites |= admitted > 0;
        }
        if (!new_sites) break;
        assign_hap_based_on_germline_het_vars_kmeans(chunk, msa_opts, kCandGermlineVarCate);
    }
}

} // namespace pgphase_collect
