#include "gap_recovery.hpp"

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
    for (size_t ci = 0; ci < chunks.size(); ++ci)
        for (size_t ri = 0; ri < chunks[ci].reads.size(); ++ri) {
            const auto& read = chunks[ci].reads[ri];
            const Key key{read.input_index, read.qname};
            reads[key].emplace_back(ci, ri);
            if (!read.is_skipped && chunks[ci].haps[ri] != 0 &&
                chunks[ci].phase_sets[ri] >= 0) {
                assignments.emplace(
                    key, std::make_pair(chunks[ci].haps[ri],
                                        chunks[ci].phase_sets[ri]));
            }
        }
}

GapStitchResult stitch_gap_proposal(std::vector<PhasingChunk>& chunks,
                                    const PhasingChunk& proposal,
                                    const PhaseGap& gap, const Options& opts,
                                    const GapReadIndex* read_index,
                                    bool defer_phase_set_merge,
                                    bool orientation_only) {
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
    std::set<hts_pos_t> supported_phase_sets;
    for (const auto& [key, locations] : index.reads) {
        const auto assignment = first_assignment(key, locations);
        if (assignment.first != 0) supported_phase_sets.insert(assignment.second);
    }
    // Only proposal reads are queried below. Retain the same first valid
    // assignment in chunk order without rebuilding a chromosome-sized name map.
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
    for (const auto& [ps, sides] : votes) {
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
