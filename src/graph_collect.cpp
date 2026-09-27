#include "graph_collect.hpp"
#include "collect_pipeline.hpp"

#include "arg_parse.hpp"
#include "collect_output.hpp"
#include "collect_phase.hpp"
#include "collect_phase_pgbam.hpp"
#include "collect_types.hpp"
#include "collect_var.hpp"
#include "gbz_ffi.h"
#include "graph_bam_adapter.hpp"
#include "graph_query.hpp"
#include "graph_sites.hpp"
#include "noise_filter.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <atomic>
#include <cstdio>
#include <cstdint>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <memory>
#include <map>
#include <set>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <unistd.h>

#include <htslib/faidx.h>
#include <htslib/sam.h>

namespace pgphase_collect {

namespace {

// Graph alignments encode path support directly, and the chr20 gap audit
// found useful linking reads between MAPQ 5 and 29. Keep this separate
// from longcallD/BAM, whose established default remains MAPQ 30.
constexpr int kDefaultGraphMinMapq = 5;
// One megabase gives the graph and its independent BAM solve enough shared
// context to stabilize read assignments while retaining bounded worker memory.
constexpr hts_pos_t kDefaultGraphRecoveryChunkSize = 1000000;

static bam_hdr_t* build_synthetic_header(faidx_t* fai) {
    bam_hdr_t* hdr = sam_hdr_init();
    if (!hdr) throw std::runtime_error("failed to allocate synthetic BAM header");
    const int n_seq = faidx_nseq(fai);
    for (int i = 0; i < n_seq; ++i) {
        const char* name = faidx_iseq(fai, i);
        const hts_pos_t len = faidx_seq_len(fai, name);
        const std::string len_str = std::to_string(len);
        if (sam_hdr_add_line(hdr, "SQ", "SN", name, "LN", len_str.c_str(), nullptr) < 0) {
            sam_hdr_destroy(hdr);
            throw std::runtime_error(std::string("failed to add SQ line for ") + name);
        }
    }
    return hdr;
}

static std::optional<bool> local_run_boundary_flip(
        const PhasingChunk& chunk, size_t source_i, size_t target_i);

// Before relabeling an established graph block, require molecules to link
// its first and last shared sites without reversing their allele orientation.
static bool local_graph_block_supported(
        const PhasingChunk& chunk,
        const std::vector<const RecoverySourceSite*>& sites,
        hts_pos_t graph_ps) {
    std::optional<size_t> first;
    std::optional<size_t> last;
    for (const RecoverySourceSite* site : sites) {
        if (!site->can_adopt || site->graph_phase_set != graph_ps ||
            site->candidate_index >= chunk.candidates.size() ||
            chunk.candidates[site->candidate_index].phase_set != graph_ps)
            continue;
        if (!first ||
            chunk.candidates[site->candidate_index].key.sort_pos() <
                chunk.candidates[*first].key.sort_pos())
            first = site->candidate_index;
        if (!last ||
            chunk.candidates[site->candidate_index].key.sort_pos() >
                chunk.candidates[*last].key.sort_pos())
            last = site->candidate_index;
    }
    if (!first || !last) return false;
    if (*first == *last) return true;
    if (chunk.candidates[*first].key.sort_pos() >=
        chunk.candidates[*last].key.sort_pos())
        return false;
    const std::optional<bool> flip =
        local_run_boundary_flip(chunk, *first, *last);
    return flip && !*flip;
}

// Salvage only a locally supported run from a source block that has a weak
// cut elsewhere. Keep the original source PS on all other sites and reads.
// A run can join two adjacent graph blocks only when both are fully covered
// inside that run and each block has direct support between its end sites.
static void adopt_local_bam_source_runs(
        PhasingChunk& chunk, hts_pos_t source_ps,
        const std::vector<const RecoverySourceSite*>& sites,
        const std::vector<const RecoverySourceRead*>& reads,
        const std::vector<hts_pos_t>& weak_cuts,
        const std::set<hts_pos_t>& approved_graph_ps,
        const std::set<size_t>& eligible_indices,
        const std::map<hts_pos_t, size_t>& local_graph_component,
        const std::vector<RecoverySeam>& seams) {
    std::vector<const RecoverySourceSite*> usable;
    for (const RecoverySourceSite* site : sites)
        if (site->can_adopt && site->candidate_index < chunk.candidates.size())
            usable.push_back(site);
    constexpr size_t kMinSupportedSourceSites = 2;
    if (usable.size() < kMinSupportedSourceSites || weak_cuts.empty()) return;
    std::sort(usable.begin(), usable.end(),
              [&chunk](const RecoverySourceSite* a,
                       const RecoverySourceSite* b) {
                  const hts_pos_t a_pos =
                      chunk.candidates[a->candidate_index].key.sort_pos();
                  const hts_pos_t b_pos =
                      chunk.candidates[b->candidate_index].key.sort_pos();
                  return a_pos == b_pos ? a->candidate_index < b->candidate_index
                                        : a_pos < b_pos;
              });
    std::vector<int> local_index(chunk.candidates.size(), -1);
    for (size_t i = 0; i < usable.size(); ++i)
        local_index[usable[i]->candidate_index] = static_cast<int>(i);
    const size_t component_count = weak_cuts.size() + 1;
    std::vector<size_t> component_of(usable.size());
    for (size_t i = 0; i < usable.size(); ++i) {
        const hts_pos_t pos =
            chunk.candidates[usable[i]->candidate_index].key.sort_pos();
        component_of[i] = static_cast<size_t>(
            std::lower_bound(weak_cuts.begin(), weak_cuts.end(), pos) -
            weak_cuts.begin());
    }

    const auto flip_hap = [](int hap) {
        return hap == 1 ? 2 : (hap == 2 ? 1 : hap);
    };
    for (size_t component = 0; component < component_count; ++component) {
        std::map<hts_pos_t, int> root_parity;
        size_t component_sites = 0;
        for (size_t i = 0; i < usable.size(); ++i) {
            if (component_of[i] != component) continue;
            ++component_sites;
            const RecoverySourceSite& site = *usable[i];
            if (site.graph_phase_set <= 0 ||
                approved_graph_ps.count(site.graph_phase_set) == 0)
                continue;
            const CandidateVariant& candidate =
                chunk.candidates[site.candidate_index];
            if (candidate.phase_set <= 0 || candidate.phase_set == source_ps ||
                candidate.bam_injected)
                continue;
            const int parity =
                candidate.hap_to_cons_alle[1] == site.hap1_allele ? 0 :
                candidate.hap_to_cons_alle[2] == site.hap1_allele ? 1 : -1;
            const auto [it, inserted] =
                root_parity.emplace(candidate.phase_set, parity);
            if (!inserted && it->second != parity) it->second = -1;
        }
        if (component_sites < kMinSupportedSourceSites ||
            root_parity.empty() || root_parity.size() > 2 ||
            std::any_of(root_parity.begin(), root_parity.end(),
                        [](const auto& root) { return root.second < 0; }))
            continue;
        if (root_parity.size() == 2) {
            const auto left = root_parity.begin();
            const auto right = std::next(left);
            const auto left_component =
                local_graph_component.find(left->first);
            const auto right_component =
                local_graph_component.find(right->first);
            if (left_component == local_graph_component.end() ||
                right_component == local_graph_component.end() ||
                left_component->second != component ||
                right_component->second != component)
                continue;
            const bool adjacent_seam = std::any_of(
                seams.begin(), seams.end(), [&](const RecoverySeam& seam) {
                    return seam.left_phase_set == left->first &&
                           seam.right_phase_set == right->first;
                });
            if (!adjacent_seam) continue;
            if (!local_graph_block_supported(chunk, usable, left->first) ||
                !local_graph_block_supported(chunk, usable, right->first))
                continue;
            const bool flip_right = left->second != right->second;
            const hts_pos_t left_ps = left->first;
            const hts_pos_t right_ps = right->first;
            for (CandidateVariant& candidate : chunk.candidates) {
                if (candidate.phase_set != right_ps) continue;
                if (flip_right) {
                    std::swap(candidate.hap_to_cons_alle[1],
                              candidate.hap_to_cons_alle[2]);
                    std::swap(candidate.hap_to_alle_profile[1],
                              candidate.hap_to_alle_profile[2]);
                    candidate.hap_alt = flip_hap(candidate.hap_alt);
                    candidate.hap_ref = flip_hap(candidate.hap_ref);
                }
                candidate.phase_set = left_ps;
            }
            for (size_t ri = 0; ri < chunk.reads.size(); ++ri) {
                if (chunk.phase_sets[ri] != right_ps) continue;
                if (flip_right)
                    chunk.haps[ri] = flip_hap(chunk.haps[ri]);
                chunk.phase_sets[ri] = left_ps;
            }
            root_parity.erase(right_ps);
        }
        if (root_parity.size() != 1) continue;
        const hts_pos_t root_ps = root_parity.begin()->first;
        const bool flip = root_parity.begin()->second == 1;
        std::set<size_t> adopted;
        bool has_new_site = false;
        for (size_t i = 0; i < usable.size(); ++i) {
            if (component_of[i] != component) continue;
            const RecoverySourceSite& site = *usable[i];
            if (eligible_indices.count(site.candidate_index) == 0 ||
                (site.graph_phase_set > 0 &&
                 approved_graph_ps.count(site.graph_phase_set) == 0))
                continue;
            const CandidateVariant& candidate =
                chunk.candidates[site.candidate_index];
            if (candidate.phase_set > 0 && candidate.phase_set != source_ps &&
                candidate.phase_set != root_ps)
                continue;
            adopted.insert(site.candidate_index);
            has_new_site |= candidate.phase_set != root_ps;
        }
        if (!has_new_site) continue;
        for (size_t i = 0; i < usable.size(); ++i) {
            const RecoverySourceSite& site = *usable[i];
            if (adopted.count(site.candidate_index) == 0) continue;
            CandidateVariant& candidate =
                chunk.candidates[site.candidate_index];
            if (candidate.phase_set == root_ps) continue;
            candidate.phase_set = root_ps;
            candidate.hap_to_cons_alle[1] =
                flip ? site.hap2_allele : site.hap1_allele;
            candidate.hap_to_cons_alle[2] =
                flip ? site.hap1_allele : site.hap2_allele;
        }
        for (const RecoverySourceRead* read : reads) {
            const size_t ri = read->read_index;
            if (ri >= chunk.phase_sets.size() ||
                ri >= chunk.read_var_profile.size() ||
                (chunk.phase_sets[ri] > 0 &&
                 chunk.phase_sets[ri] != source_ps))
                continue;
            const ReadVariantProfile& profile = chunk.read_var_profile[ri];
            if (profile.start_var_idx < 0) continue;
            bool observes_adopted = false;
            bool observes_other_component = false;
            for (size_t offset = 0; offset < profile.alleles.size(); ++offset) {
                const size_t ci =
                    static_cast<size_t>(profile.start_var_idx) + offset;
                if (ci >= local_index.size()) break;
                if (profile.alleles[offset] < 0 || local_index[ci] < 0)
                    continue;
                const size_t local = static_cast<size_t>(local_index[ci]);
                if (component_of[local] != component)
                    observes_other_component = true;
                else if (adopted.count(ci) != 0)
                    observes_adopted = true;
            }
            if (!observes_adopted || observes_other_component) continue;
            chunk.haps[ri] = flip ? flip_hap(read->hap) : read->hap;
            chunk.phase_sets[ri] = root_ps;
        }
    }
}

// A shared clean site establishes the allele gauge, while spanning reads test
// whether both source haplotypes agree with that gauge despite occasional errors.
static bool source_graph_vote_supported(int same, int cross, int parity) {
    constexpr double kSourceGraphMaxP = 0.01;
    if (same == cross || (parity == 0 ? same < cross : cross < same))
        return false;
    const int winner = std::max(same, cross);
    const int total = same + cross;
    double term = std::exp(std::lgamma(static_cast<double>(total + 1)) -
                           std::lgamma(static_cast<double>(winner + 1)) -
                           std::lgamma(static_cast<double>(total - winner + 1)) -
                           static_cast<double>(total) * std::log(2.0));
    double tail = term;
    for (int k = winner; k < total; ++k) {
        term *= static_cast<double>(total - k) /
                static_cast<double>(k + 1);
        tail += term;
    }
    return tail <= kSourceGraphMaxP;
}

// A graph block may be reused at a second seam when one complete BAM source
// path corroborates its two endpoints and every shared site's polarity. This
// supplies the same end-to-end certificate as a single spanning read without
// requiring one molecule to cover the entire established graph block.
static std::set<hts_pos_t> complete_graph_source_paths(
        const GraphChunkBuildResult& gc) {
    struct GraphExtent {
        size_t first = 0;
        size_t last = 0;
        size_t site_count = 0;
    };
    std::map<hts_pos_t, GraphExtent> extents;
    for (size_t ci = 0; ci < gc.chunk.candidates.size(); ++ci) {
        const CandidateVariant& candidate = gc.chunk.candidates[ci];
        if (candidate.phase_set <= 0 || candidate.bam_injected ||
            candidate.hap_to_cons_alle[1] < 0 ||
            candidate.hap_to_cons_alle[2] < 0 ||
            candidate.hap_to_cons_alle[1] == candidate.hap_to_cons_alle[2])
            continue;
        const auto [it, inserted] = extents.try_emplace(
            candidate.phase_set, GraphExtent{ci, ci, 0});
        (void)inserted;
        it->second.last = ci;
        ++it->second.site_count;
    }
    struct SharedPath {
        int parity = -1;
        bool conflict = false;
        bool first = false;
        bool last = false;
        size_t shared_sites = 0;
    };
    std::map<std::pair<hts_pos_t, hts_pos_t>, SharedPath> paths;
    for (const RecoverySourceSite& site : gc.recovery_source_sites) {
        const auto extent = extents.find(site.graph_phase_set);
        const auto complete = gc.recovery_source_path_supported.find(site.phase_set);
        if (extent == extents.end() || complete == gc.recovery_source_path_supported.end() ||
            !complete->second || !site.can_adopt ||
            site.candidate_index >= gc.chunk.candidates.size() ||
            site.graph_hap1_allele < 0 || site.graph_hap1_allele > 1 ||
            site.hap1_allele < 0 || site.hap1_allele > 1)
            continue;
        SharedPath& path = paths[{site.graph_phase_set, site.phase_set}];
        const int parity = site.graph_hap1_allele != site.hap1_allele;
        if (path.parity >= 0 && path.parity != parity) path.conflict = true;
        path.parity = parity;
        path.first |= site.candidate_index == extent->second.first;
        path.last |= site.candidate_index == extent->second.last;
        ++path.shared_sites;
    }
    std::set<hts_pos_t> certified;
    for (const auto& [key, path] : paths) {
        const auto extent = extents.find(key.first);
        if (extent == extents.end() || extent->second.site_count < 2 ||
            gc.chunk.candidates[extent->second.first].key.sort_pos() >=
                gc.chunk.candidates[extent->second.last].key.sort_pos() ||
            path.conflict || !path.first || !path.last ||
            path.shared_sites < 2)
            continue;
        for (const RecoveryPhaseGauge& gauge : gc.recovery_phase_gauges) {
            const auto vote = std::find_if(
                gauge.block_votes.begin(), gauge.block_votes.end(),
                [&](const RecoveryBlockGaugeVote& candidate) {
                    return candidate.graph_phase_set == key.first &&
                           candidate.bam_phase_set == key.second;
                });
            if (vote == gauge.block_votes.end()) continue;
            const bool both_haps =
                vote->counts[0][0] + vote->counts[0][1] > 0 &&
                vote->counts[1][0] + vote->counts[1][1] > 0;
            const int same = vote->counts[0][0] + vote->counts[1][1];
            const int cross = vote->counts[0][1] + vote->counts[1][0];
            if (both_haps && source_graph_vote_supported(same, cross, path.parity))
                certified.insert(key.first);
        }
    }
    return certified;
}

// Attach complete BAM source blocks without splitting established graph blocks.
// A source with weak cuts can transfer only its independently supported runs.
// One run may bridge adjacent graph blocks when both flank votes and local
// graph-block support agree; no label crosses the source's weak cut.
static void attach_bam_source_phase_sets(
        GraphChunkBuildResult& graph_chunk,
        const std::set<hts_pos_t>& locally_bridged_sources) {
    PhasingChunk& chunk = graph_chunk.chunk;
    std::map<hts_pos_t, std::vector<const RecoverySourceSite*>> sites_by_ps;
    std::map<hts_pos_t, std::vector<const RecoverySourceRead*>> reads_by_ps;
    for (const RecoverySourceSite& site : graph_chunk.recovery_source_sites)
        if (site.candidate_index < chunk.candidates.size() && site.phase_set > 0)
            sites_by_ps[site.phase_set].push_back(&site);
    for (const RecoverySourceRead& read : graph_chunk.recovery_source_reads)
        if (read.read_index < chunk.reads.size() && read.phase_set > 0)
            reads_by_ps[read.phase_set].push_back(&read);

    for (const auto& [source_ps, sites] : sites_by_ps) {
        if (locally_bridged_sources.count(source_ps) != 0) continue;
        if (sites.size() < 2) continue;
        std::map<hts_pos_t, int> graph_parity;
        for (const RecoverySourceSite* site : sites) {
            if (!site->can_adopt || site->graph_phase_set <= 0) continue;
            const int parity =
                site->graph_hap1_allele == site->hap1_allele ? 0 : 1;
            const auto [it, inserted] = graph_parity.emplace(
                site->graph_phase_set, parity);
            if (!inserted && it->second != parity) it->second = -1;
        }
        const auto source_path =
            graph_chunk.recovery_source_path_supported.find(source_ps);
        const bool complete_source_path =
            source_path != graph_chunk.recovery_source_path_supported.end() &&
            source_path->second;
        const auto weak = graph_chunk.recovery_source_weak_cuts.find(source_ps);
        std::map<hts_pos_t, size_t> local_graph_component;
        if (!complete_source_path &&
            weak != graph_chunk.recovery_source_weak_cuts.end() &&
            !weak->second.empty()) {
            std::map<hts_pos_t, size_t> graph_site_counts;
            for (const CandidateVariant& candidate : chunk.candidates) {
                if (graph_parity.count(candidate.phase_set) != 0 &&
                    !candidate.bam_injected &&
                    candidate.hap_to_cons_alle[1] >= 0 &&
                    candidate.hap_to_cons_alle[2] >= 0 &&
                    candidate.hap_to_cons_alle[1] !=
                        candidate.hap_to_cons_alle[2])
                    ++graph_site_counts[candidate.phase_set];
            }
            std::map<hts_pos_t, std::set<size_t>> matched_graph_sites;
            std::set<hts_pos_t> component_conflicts;
            for (const RecoverySourceSite* site : sites) {
                if (!site->can_adopt || site->graph_phase_set <= 0 ||
                    site->candidate_index >= chunk.candidates.size())
                    continue;
                const CandidateVariant& candidate =
                    chunk.candidates[site->candidate_index];
                if (candidate.phase_set != site->graph_phase_set ||
                    candidate.bam_injected)
                    continue;
                const size_t component = static_cast<size_t>(
                    std::lower_bound(
                        weak->second.begin(), weak->second.end(),
                        candidate.key.sort_pos()) - weak->second.begin());
                const auto [it, inserted] = local_graph_component.try_emplace(
                    site->graph_phase_set, component);
                if (!inserted && it->second != component)
                    component_conflicts.insert(site->graph_phase_set);
                matched_graph_sites[site->graph_phase_set].insert(
                    site->candidate_index);
            }
            for (auto it = local_graph_component.begin();
                 it != local_graph_component.end();) {
                if (component_conflicts.count(it->first) != 0 ||
                    matched_graph_sites[it->first].size() !=
                        graph_site_counts[it->first])
                    it = local_graph_component.erase(it);
                else ++it;
            }
        }
        std::set<hts_pos_t> locally_bridgeable_graph_ps;
        for (const RecoverySeam& seam : graph_chunk.recovery_windows) {
            const auto left = local_graph_component.find(seam.left_phase_set);
            const auto right = local_graph_component.find(seam.right_phase_set);
            if (left == local_graph_component.end() ||
                right == local_graph_component.end() ||
                left->second != right->second ||
                !local_graph_block_supported(
                    chunk, sites, seam.left_phase_set) ||
                !local_graph_block_supported(
                    chunk, sites, seam.right_phase_set))
                continue;
            locally_bridgeable_graph_ps.insert(seam.left_phase_set);
            locally_bridgeable_graph_ps.insert(seam.right_phase_set);
        }
        std::set<hts_pos_t> approved_graph_ps;
        for (const auto& [graph_ps, parity] : graph_parity) {
            if (parity < 0) continue;
            for (const RecoveryPhaseGauge& gauge : graph_chunk.recovery_phase_gauges) {
                const auto vote = std::find_if(
                    gauge.block_votes.begin(), gauge.block_votes.end(),
                    [graph_ps, source_ps](const RecoveryBlockGaugeVote& v) {
                        return v.graph_phase_set == graph_ps &&
                               v.bam_phase_set == source_ps;
                    });
                if (vote == gauge.block_votes.end()) continue;
                const int same = vote->counts[0][0] + vote->counts[1][1];
                const int cross = vote->counts[0][1] + vote->counts[1][0];
                // Both source haplotypes must independently identify the same
                // local orientation. A conflicting read leaves the graph block
                // in its original phase set for this first transfer pass.
                const bool both_haps =
                    vote->counts[0][0] + vote->counts[0][1] > 0 &&
                    vote->counts[1][0] + vote->counts[1][1] > 0;
                // A complete source path keeps the same orientation across
                // all its sites, so one read error can be judged statistically.
                // Only a fully covered pair of adjacent graph blocks can use
                // a statistical vote from a weak-cut source. Other runs keep
                // the conflict-free rule for a one-flank attachment.
                const bool supported_local_component =
                    locally_bridgeable_graph_ps.count(graph_ps) != 0;
                const bool approved =
                    complete_source_path || supported_local_component
                        ? source_graph_vote_supported(same, cross, parity)
                        : ((parity == 0 && same > 0 && cross == 0) ||
                           (parity == 1 && cross > 0 && same == 0));
                if (both_haps && approved)
                    approved_graph_ps.insert(graph_ps);
                break;
            }
        }
        // One well-supported flank is enough to orient this source block.
        // The opposite graph flank stays independent unless it passes its own
        // exact-site and two-haplotype read checks.
        std::vector<const RecoverySeam*> attachable_seams;
        for (const RecoverySeam& seam : graph_chunk.recovery_windows)
            if (approved_graph_ps.count(seam.left_phase_set) != 0 ||
                approved_graph_ps.count(seam.right_phase_set) != 0)
                attachable_seams.push_back(&seam);
        if (attachable_seams.empty()) continue;

        const auto source_reads = reads_by_ps.find(source_ps);
        if (source_reads == reads_by_ps.end()) continue;

        std::set<size_t> eligible_indices;
        for (const RecoverySourceSite* site : sites) {
            if (!site->can_adopt) continue;
            if (site->graph_phase_set > 0 &&
                approved_graph_ps.count(site->graph_phase_set) == 0)
                continue;
            const hts_pos_t pos =
                chunk.candidates[site->candidate_index].key.sort_pos();
            const bool inside_attachable_seam = std::any_of(
                attachable_seams.begin(), attachable_seams.end(),
                [pos](const RecoverySeam* seam) {
                    return seam->beg <= pos && pos <= seam->end;
                });
            if (inside_attachable_seam)
                eligible_indices.insert(site->candidate_index);
        }

        if (source_path == graph_chunk.recovery_source_path_supported.end() ||
            !source_path->second) {
            const auto weak_cuts =
                graph_chunk.recovery_source_weak_cuts.find(source_ps);
            if (weak_cuts != graph_chunk.recovery_source_weak_cuts.end())
                adopt_local_bam_source_runs(chunk, source_ps, sites,
                                            source_reads->second,
                                            weak_cuts->second,
                                            approved_graph_ps, eligible_indices,
                                            local_graph_component,
                                            graph_chunk.recovery_windows);
            continue;
        }

        // Earlier seams may have relabeled and flipped the graph block. Use
        // its current PS and allele orientation, not the pre-stitch metadata,
        // then move the entire current block into the supported BAM gauge.
        std::map<hts_pos_t, int> current_parity;
        for (const RecoverySourceSite* site : sites) {
            if (!site->can_adopt || site->graph_phase_set <= 0 ||
                approved_graph_ps.count(site->graph_phase_set) == 0)
                continue;
            const CandidateVariant& candidate =
                chunk.candidates[site->candidate_index];
            if (candidate.phase_set <= 0 || candidate.phase_set == source_ps)
                continue;
            const int parity =
                candidate.hap_to_cons_alle[1] == site->hap1_allele ? 0 :
                candidate.hap_to_cons_alle[2] == site->hap1_allele ? 1 : -1;
            const auto [it, inserted] =
                current_parity.emplace(candidate.phase_set, parity);
            if (!inserted && it->second != parity) it->second = -1;
        }
        std::map<hts_pos_t, int> transferred_graph_ps;
        for (const RecoverySourceSite* site : sites) {
            if (site->graph_phase_set <= 0 ||
                eligible_indices.count(site->candidate_index) == 0)
                continue;
            const hts_pos_t current_ps =
                chunk.candidates[site->candidate_index].phase_set;
            const auto parity = current_parity.find(current_ps);
            if (parity != current_parity.end() && parity->second >= 0)
                transferred_graph_ps.emplace(current_ps, parity->second);
        }
        std::set<size_t> adopted_indices;
        for (const RecoverySourceSite* site : sites) {
            if (eligible_indices.count(site->candidate_index) == 0) continue;
            const hts_pos_t current_ps =
                chunk.candidates[site->candidate_index].phase_set;
            if (current_ps > 0 && current_ps != source_ps &&
                transferred_graph_ps.count(current_ps) == 0)
                continue;
            adopted_indices.insert(site->candidate_index);
        }
        const auto flip_hap = [](int hap) {
            return hap == 1 ? 2 : (hap == 2 ? 1 : hap);
        };
        for (const auto& [graph_ps, parity] : transferred_graph_ps) {
            const bool flip = parity == 1;
            for (CandidateVariant& candidate : chunk.candidates) {
                if (candidate.phase_set != graph_ps) continue;
                if (flip) {
                    std::swap(candidate.hap_to_cons_alle[1],
                              candidate.hap_to_cons_alle[2]);
                    std::swap(candidate.hap_to_alle_profile[1],
                              candidate.hap_to_alle_profile[2]);
                    candidate.hap_alt = flip_hap(candidate.hap_alt);
                    candidate.hap_ref = flip_hap(candidate.hap_ref);
                }
                candidate.phase_set = source_ps;
            }
            for (size_t ri = 0; ri < chunk.reads.size(); ++ri) {
                if (chunk.phase_sets[ri] != graph_ps) continue;
                if (flip) chunk.haps[ri] = flip_hap(chunk.haps[ri]);
                chunk.phase_sets[ri] = source_ps;
            }
        }
        for (const RecoverySourceSite* site : sites) {
            if (adopted_indices.count(site->candidate_index) == 0) continue;
            CandidateVariant& candidate =
                chunk.candidates[site->candidate_index];
            candidate.phase_set = source_ps;
            candidate.hap_to_cons_alle[1] = site->hap1_allele;
            candidate.hap_to_cons_alle[2] = site->hap2_allele;
        }
        if (adopted_indices.empty()) continue;
        for (const RecoverySourceRead* read : source_reads->second) {
            const size_t ri = read->read_index;
            const hts_pos_t current_ps = chunk.phase_sets[ri];
            if (current_ps > 0 && current_ps != source_ps) continue;
            const ReadVariantProfile& profile = chunk.read_var_profile[ri];
            if (profile.start_var_idx < 0) continue;
            bool observes_adopted_site = false;
            for (size_t offset = 0; offset < profile.alleles.size(); ++offset) {
                if (profile.alleles[offset] >= 0 &&
                    adopted_indices.count(
                        static_cast<size_t>(profile.start_var_idx) + offset) != 0) {
                    observes_adopted_site = true;
                    break;
                }
            }
            if (!observes_adopted_site) continue;
            chunk.haps[ri] = read->hap;
            chunk.phase_sets[ri] = source_ps;
        }
    }
}

// Test one preselected boundary pair, counting each read once. The candidate
// pair is chosen by coordinate before observing its alleles; the two disjoint
// read halves guard against a bridge driven by one read cohort.
static std::optional<bool> local_run_boundary_flip(
        const PhasingChunk& chunk, size_t source_i, size_t target_i) {
    constexpr double kLocalRunMaxP = 0.05;
    std::array<std::pair<int, int>, 2> halves{};
    std::array<int, 2> source_alleles{};
    const CandidateVariant& source = chunk.candidates[source_i];
    const CandidateVariant& target = chunk.candidates[target_i];
    for (size_t ri = 0; ri < chunk.read_var_profile.size(); ++ri) {
        if (ri >= chunk.reads.size() || chunk.reads[ri].is_skipped) continue;
        const ReadVariantProfile& profile = chunk.read_var_profile[ri];
        if (profile.start_var_idx < 0 ||
            source_i < static_cast<size_t>(profile.start_var_idx) ||
            target_i < static_cast<size_t>(profile.start_var_idx) ||
            source_i > static_cast<size_t>(profile.end_var_idx) ||
            target_i > static_cast<size_t>(profile.end_var_idx))
            continue;
        const size_t source_offset =
            source_i - static_cast<size_t>(profile.start_var_idx);
        const size_t target_offset =
            target_i - static_cast<size_t>(profile.start_var_idx);
        if (source_offset >= profile.alleles.size() ||
            target_offset >= profile.alleles.size())
            continue;
        const int source_allele = profile.alleles[source_offset];
        const int target_allele = profile.alleles[target_offset];
        if ((source_allele != 0 && source_allele != 1) ||
            (target_allele != 0 && target_allele != 1))
            continue;
        ++source_alleles[static_cast<size_t>(source_allele)];
        uint64_t hash = 14695981039346656037ULL;
        for (const unsigned char byte : chunk.reads[ri].qname) {
            hash ^= byte;
            hash *= 1099511628211ULL;
        }
        auto& votes = halves[static_cast<size_t>(hash & 1ULL)];
        if ((source_allele == source.hap_to_cons_alle[1]) ==
            (target_allele == target.hap_to_cons_alle[1]))
            ++votes.first;
        else
            ++votes.second;
    }
    if (source_alleles[0] == 0 || source_alleles[1] == 0) return std::nullopt;
    const int same = halves[0].first + halves[1].first;
    const int cross = halves[0].second + halves[1].second;
    const int winner = std::max(same, cross);
    const int total = same + cross;
    if (winner == 0 || same == cross ||
        halves[0].first == halves[0].second ||
        halves[1].first == halves[1].second)
        return std::nullopt;
    const bool flip = cross > same;
    for (const auto& half : halves)
        if ((half.second > half.first) != flip) return std::nullopt;
    double term = std::exp(
        std::lgamma(static_cast<double>(total + 1)) -
        std::lgamma(static_cast<double>(winner + 1)) -
        std::lgamma(static_cast<double>(total - winner + 1)) -
        static_cast<double>(total) * std::log(2.0));
    double tail = term;
    for (int k = winner; k < total; ++k) {
        term *= static_cast<double>(total - k) /
                static_cast<double>(k + 1);
        tail += term;
    }
    return tail <= kLocalRunMaxP ? std::optional<bool>(flip) : std::nullopt;
}

// A one-site BAM block has no exact shared allele to orient a neighboring
// graph block. The closest graph SNP can itself have a few bad observations;
// use the next preselected clean SNP only when molecules independently
// confirm its orientation to both the boundary SNP and the block endpoint.
static void attach_singleton_bam_sources(GraphChunkBuildResult& gc) {
    PhasingChunk& chunk = gc.chunk;
    std::map<hts_pos_t, std::vector<const RecoverySourceSite*>> sources;
    std::map<hts_pos_t, std::vector<size_t>> graph_snps;
    for (const RecoverySourceSite& site : gc.recovery_source_sites)
        if (site.phase_set > 0) sources[site.phase_set].push_back(&site);
    for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
        const CandidateVariant& candidate = chunk.candidates[ci];
        if (candidate.phase_set > 0 && !candidate.bam_injected &&
            candidate.counts.category == VariantCategory::CleanHetSnp &&
            candidate.hap_to_cons_alle[1] >= 0 &&
            candidate.hap_to_cons_alle[2] >= 0 &&
            candidate.hap_to_cons_alle[1] != candidate.hap_to_cons_alle[2])
            graph_snps[candidate.phase_set].push_back(ci);
    }
    for (const auto& [source_ps, sites] : sources) {
        if (sites.size() != 1 || !sites.front()->can_adopt ||
            sites.front()->candidate_index >= chunk.candidates.size()) continue;
        const size_t source_i = sites.front()->candidate_index;
        CandidateVariant& source = chunk.candidates[source_i];
        if (source.phase_set != source_ps || !source.bam_injected ||
            !source.msa_verified || !source.alignment_verified ||
            source.hap_to_cons_alle[1] < 0 ||
            source.hap_to_cons_alle[2] < 0 ||
            source.hap_to_cons_alle[1] == source.hap_to_cons_alle[2])
            continue;

        const hts_pos_t source_pos = source.key.sort_pos();
        std::optional<hts_pos_t> nearest_distance;
        hts_pos_t target_ps = 0;
        bool source_is_left = false;
        for (const auto& [graph_ps, indices] : graph_snps) {
            if (indices.size() < 3) continue;
            const hts_pos_t first =
                chunk.candidates[indices.front()].key.sort_pos();
            const hts_pos_t last =
                chunk.candidates[indices.back()].key.sort_pos();
            if (source_pos >= first && source_pos <= last) continue;
            const hts_pos_t distance =
                source_pos < first ? first - source_pos : source_pos - last;
            if (!nearest_distance || distance < *nearest_distance) {
                nearest_distance = distance;
                target_ps = graph_ps;
                source_is_left = source_pos < first;
            } else if (distance == *nearest_distance) {
                target_ps = 0;
            }
        }
        if (target_ps <= 0 || target_ps == source_ps) continue;
        const auto& indices = graph_snps.at(target_ps);
        const size_t boundary = source_is_left
            ? indices.front() : indices.back();
        const size_t next = source_is_left
            ? indices[1] : indices[indices.size() - 2];
        const size_t far = source_is_left
            ? indices.back() : indices.front();
        const hts_pos_t boundary_pos =
            chunk.candidates[boundary].key.sort_pos();
        const hts_pos_t next_pos = chunk.candidates[next].key.sort_pos();
        const hts_pos_t far_pos = chunk.candidates[far].key.sort_pos();
        if (boundary_pos == next_pos || next_pos == far_pos) continue;
        const std::optional<bool> first_link =
            local_run_boundary_flip(chunk, source_i, boundary);
        const std::optional<bool> second_link =
            local_run_boundary_flip(chunk, source_i, next);
        const std::optional<bool> local_graph_link =
            local_run_boundary_flip(chunk, boundary, next);
        const std::optional<bool> whole_graph_link =
            local_run_boundary_flip(chunk, boundary, far);
        if (!second_link || !local_graph_link || *local_graph_link ||
            !whole_graph_link || *whole_graph_link ||
            (first_link && *first_link != *second_link)) continue;

        source.phase_set = target_ps;
        if (*second_link) {
            std::swap(source.hap_to_cons_alle[1], source.hap_to_cons_alle[2]);
            std::swap(source.hap_to_alle_profile[1],
                      source.hap_to_alle_profile[2]);
            std::swap(source.hap_alt, source.hap_ref);
        }
        for (const RecoverySourceRead& read : gc.recovery_source_reads) {
            if (read.phase_set != source_ps ||
                (read.hap != 1 && read.hap != 2) ||
                read.read_index >= chunk.phase_sets.size() ||
                read.read_index >= chunk.read_var_profile.size() ||
                chunk.phase_sets[read.read_index] != source_ps) continue;
            const ReadVariantProfile& profile =
                chunk.read_var_profile[read.read_index];
            if (profile.start_var_idx < 0 ||
                source_i < static_cast<size_t>(profile.start_var_idx) ||
                source_i > static_cast<size_t>(profile.end_var_idx)) continue;
            const size_t offset =
                source_i - static_cast<size_t>(profile.start_var_idx);
            if (offset >= profile.alleles.size() ||
                (profile.alleles[offset] != 0 &&
                 profile.alleles[offset] != 1)) continue;
            chunk.phase_sets[read.read_index] = target_ps;
            chunk.haps[read.read_index] = *second_link
                ? (read.hap == 1 ? 2 : read.hap == 2 ? 1 : 0)
                : read.hap;
        }
    }
}

// A weak BAM source block may still have a supported local run at one end.
// Attach such a run only when exactly one neighboring phase set has a decisive
// molecule link. This never relabels the distant side of the weak source cut.
static void attach_unanchored_bam_source_runs(
        GraphChunkBuildResult& gc,
        const std::set<hts_pos_t>& locally_bridged_sources,
        const std::set<hts_pos_t>* only_sources = nullptr) {
    PhasingChunk& chunk = gc.chunk;
    std::map<hts_pos_t, std::vector<const RecoverySourceSite*>> by_source;
    for (const RecoverySourceSite& site : gc.recovery_source_sites)
        if (site.phase_set > 0 && site.candidate_index < chunk.candidates.size())
            by_source[site.phase_set].push_back(&site);
    const auto verified = [](const CandidateVariant& candidate) {
        return candidate.counts.category == VariantCategory::CleanHetSnp ||
               candidate.counts.category == VariantCategory::CleanHetIndel ||
               (candidate.bam_injected && candidate.msa_verified &&
                candidate.alignment_verified);
    };
    const auto oriented = [](const CandidateVariant& candidate) {
        return candidate.phase_set > 0 &&
               candidate.hap_to_cons_alle[1] >= 0 &&
               candidate.hap_to_cons_alle[1] <= 1 &&
               candidate.hap_to_cons_alle[2] >= 0 &&
               candidate.hap_to_cons_alle[2] <= 1 &&
               candidate.hap_to_cons_alle[1] != candidate.hap_to_cons_alle[2];
    };
    for (const auto& [source_ps, sites] : by_source) {
        if (locally_bridged_sources.count(source_ps) != 0 ||
            (only_sources != nullptr && only_sources->count(source_ps) == 0))
            continue;
        const auto cuts = gc.recovery_source_weak_cuts.find(source_ps);
        if (cuts == gc.recovery_source_weak_cuts.end() ||
            cuts->second.empty()) continue;
        const size_t run_count = cuts->second.size() + 1;
        std::vector<std::vector<const RecoverySourceSite*>> runs(run_count);
        for (const RecoverySourceSite* site : sites) {
            const hts_pos_t pos =
                chunk.candidates[site->candidate_index].key.sort_pos();
            const size_t run = static_cast<size_t>(std::lower_bound(
                cuts->second.begin(), cuts->second.end(), pos) -
                cuts->second.begin());
            runs[run].push_back(site);
        }
        for (size_t run = 0; run < run_count; ++run) {
            const auto& sites_in_run = runs[run];
            if (sites_in_run.empty()) continue;
            // Shared graph rows or an earlier exact-key attachment already
            // determine this run's owner; a boundary vote cannot override it.
            bool independent = true;
            for (const RecoverySourceSite* site : sites_in_run) {
                const CandidateVariant& candidate =
                    chunk.candidates[site->candidate_index];
                if (!site->can_adopt || candidate.phase_set != source_ps)
                    independent = false;
            }
            if (!independent) continue;
            size_t left_source = sites_in_run.front()->candidate_index;
            size_t right_source = left_source;
            for (const RecoverySourceSite* site : sites_in_run) {
                const size_t ci = site->candidate_index;
                const hts_pos_t pos = chunk.candidates[ci].key.sort_pos();
                if (pos < chunk.candidates[left_source].key.sort_pos() ||
                    (pos == chunk.candidates[left_source].key.sort_pos() &&
                     ci < left_source)) left_source = ci;
                if (pos > chunk.candidates[right_source].key.sort_pos() ||
                    (pos == chunk.candidates[right_source].key.sort_pos() &&
                     ci < right_source)) right_source = ci;
            }
            const hts_pos_t left_pos =
                chunk.candidates[left_source].key.sort_pos();
            const hts_pos_t right_pos =
                chunk.candidates[right_source].key.sort_pos();
            std::optional<size_t> left_target;
            std::optional<size_t> right_target;
            for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
                const CandidateVariant& target = chunk.candidates[ci];
                if (target.phase_set == source_ps || !oriented(target) ||
                    !verified(target)) continue;
                const hts_pos_t pos = target.key.sort_pos();
                if (pos < left_pos &&
                    (!left_target || pos >
                     chunk.candidates[*left_target].key.sort_pos()))
                    left_target = ci;
                if (pos > right_pos &&
                    (!right_target || pos <
                     chunk.candidates[*right_target].key.sort_pos()))
                    right_target = ci;
            }
            std::optional<std::pair<hts_pos_t, bool>> attachment;
            bool ambiguous = false;
            for (const auto [source_i, target_i] :
                 {std::make_pair(left_source, left_target),
                  std::make_pair(right_source, right_target)}) {
                if (!target_i || !verified(chunk.candidates[source_i]))
                    continue;
                const std::optional<bool> flip =
                    local_run_boundary_flip(chunk, source_i, *target_i);
                if (!flip) continue;
                const hts_pos_t target_ps =
                    chunk.candidates[*target_i].phase_set;
                if (attachment && attachment->first != target_ps) {
                    ambiguous = true;
                    break;
                }
                if (attachment && attachment->second != *flip) {
                    ambiguous = true;
                    break;
                }
                attachment = std::make_pair(target_ps, *flip);
            }
            if (ambiguous || !attachment) continue;
            // A different run from the same original source may already have
            // attached to this root. Reusing it would bypass the weak cut.
            const bool root_used_elsewhere = std::any_of(
                sites.begin(), sites.end(), [&](const RecoverySourceSite* site) {
                    const hts_pos_t pos = chunk.candidates[
                        site->candidate_index].key.sort_pos();
                    const size_t other_run = static_cast<size_t>(
                        std::lower_bound(cuts->second.begin(),
                                         cuts->second.end(), pos) -
                        cuts->second.begin());
                    return other_run != run &&
                        chunk.candidates[site->candidate_index].phase_set ==
                            attachment->first;
                });
            if (root_used_elsewhere) continue;
            std::set<size_t> adopted;
            for (const RecoverySourceSite* site : sites_in_run) {
                CandidateVariant& candidate =
                    chunk.candidates[site->candidate_index];
                adopted.insert(site->candidate_index);
                candidate.phase_set = attachment->first;
                if (attachment->second) {
                    std::swap(candidate.hap_to_cons_alle[1],
                              candidate.hap_to_cons_alle[2]);
                    std::swap(candidate.hap_to_alle_profile[1],
                              candidate.hap_to_alle_profile[2]);
                    std::swap(candidate.hap_alt, candidate.hap_ref);
                }
            }
            for (const RecoverySourceRead& read : gc.recovery_source_reads) {
                if (read.phase_set != source_ps ||
                    read.read_index >= chunk.phase_sets.size() ||
                    read.read_index >= chunk.read_var_profile.size() ||
                    chunk.phase_sets[read.read_index] != source_ps)
                    continue;
                const ReadVariantProfile& profile =
                    chunk.read_var_profile[read.read_index];
                if (profile.start_var_idx < 0) continue;
                bool observes_run = false;
                bool observes_other_run = false;
                for (size_t offset = 0; offset < profile.alleles.size(); ++offset) {
                    if (profile.alleles[offset] < 0) continue;
                    const size_t ci =
                        static_cast<size_t>(profile.start_var_idx) + offset;
                    if (adopted.count(ci) != 0) observes_run = true;
                    else if (ci < chunk.candidates.size() &&
                             chunk.candidates[ci].phase_set == source_ps)
                        observes_other_run = true;
                }
                if (!observes_run || observes_other_run) continue;
                chunk.phase_sets[read.read_index] = attachment->first;
                chunk.haps[read.read_index] = attachment->second
                    ? (read.hap == 1 ? 2 : (read.hap == 2 ? 1 : 0))
                    : read.hap;
            }
        }
    }
}

/// In-chunk recovery for one chunk.
///
/// The BAM sub-solve phases injected sites and reads in an independent local
/// gauge. Once merged, decisive allele edges stitch left graph PS -> local PS
/// -> right graph PS. Existing graph blocks are never globally re-clustered.
static void run_in_chunk_recovery(GraphChunkBuildResult& gc,
                                  const Options& opts,
                                  WorkerContext& ctx,
                                  const char* contig) {
    if (!recover_phase_set_seams_in_place(gc, opts, ctx, contig)) return;

    Options stitch_opts = opts;
    constexpr int kRecoveryLinkWindow = 128;
    stitch_opts.block_link_window = kRecoveryLinkWindow;
    // Match longcallD's two-read link floor. A single molecule may nominate an
    // edge for diagnostics, but cannot merge two production phase blocks.
    stitch_opts.min_block_link_reads = 1;
    stitch_opts.link_by_alleles = true;

    // This is the complete replay state: established graph blocks, imported
    // local blocks, read assignments and every allele observation at the seam.
    if (!opts.phase_matrix_dump_prefix.empty()) {
        for (const RecoveryPhaseGauge& gauge : gc.recovery_phase_gauges) {
            for (const PhaseSetGaugeVote& vote : gauge.graph_votes)
                std::fprintf(stderr, "[recovery-gauge] %" PRId64 "-%" PRId64 " ps=%" PRId64 " same=%d cross=%d\n",
                             static_cast<int64_t>(gauge.beg), static_cast<int64_t>(gauge.end),
                             static_cast<int64_t>(vote.phase_set), vote.same, vote.cross);
            for (const RecoveryBlockGaugeVote& vote : gauge.block_votes)
                std::fprintf(stderr, "[recovery-block-gauge] %" PRId64 "-%" PRId64
                                     " graph_ps=%" PRId64 " bam_ps=%" PRId64
                                     " counts=%d,%d,%d,%d shared_candidates=%d,%d\n",
                             static_cast<int64_t>(gauge.beg),
                             static_cast<int64_t>(gauge.end),
                             static_cast<int64_t>(vote.graph_phase_set),
                             static_cast<int64_t>(vote.bam_phase_set),
                             vote.counts[0][0], vote.counts[0][1],
                             vote.counts[1][0], vote.counts[1][1],
                             vote.shared_candidate_same,
                             vote.shared_candidate_cross);
        }
    }
    std::set<hts_pos_t> pre_attach_sources;
    for (const RecoveryPhaseGauge& gauge : gc.recovery_phase_gauges)
        for (const RecoveryPhysicalSnpBridge& bridge :
             gauge.physical_snp_bridges)
            if (bridge.pre_attach_source_phase_set > 0)
                pre_attach_sources.insert(bridge.pre_attach_source_phase_set);
    if (!pre_attach_sources.empty())
        attach_unanchored_bam_source_runs(gc, {}, &pre_attach_sources);
    const std::set<hts_pos_t> reusable_graph_paths =
        complete_graph_source_paths(gc);
    dump_recovery_phase_state(gc.chunk, stitch_opts, "recovery-input");
    std::set<hts_pos_t> locally_bridged_sources;
    stitch_recovery_phase_sets_left_to_right(
        gc.chunk, gc.recovery_windows, gc.recovery_phase_gauges, stitch_opts,
        &gc.recovery_source_path_supported,
        &gc.recovery_source_weak_cuts,
        &gc.recovery_source_quality_cuts,
        &locally_bridged_sources, &reusable_graph_paths);
    // These source blocks already oriented both graph flanks through a local
    // validated path. A later one-sided run attachment would split that join.
    attach_bam_source_phase_sets(gc, locally_bridged_sources);
    attach_unanchored_bam_source_runs(gc, locally_bridged_sources);
    attach_singleton_bam_sources(gc);
    dump_recovery_phase_state(gc.chunk, stitch_opts, "recovery-final");
}

static CandidateTable graph_chunks_to_candidate_table(
    const std::vector<GraphChunkBuildResult>& graph_chunks,
    const std::unordered_map<std::string, int>& contig_to_tid,
    const Options& opts)
{
    CandidateTable result;

    for (const GraphChunkBuildResult& graph_chunk : graph_chunks) {
        const PhasingChunk& chunk = graph_chunk.chunk;
        for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
            const CandidateVariant& mcand = chunk.candidates[ci];

            if (ci >= graph_chunk.site_meta.size()) continue;
            const GraphSiteMeta& meta = graph_chunk.site_meta[ci];

            // A site merged in by the in-chunk recovery (no catalog id) is
            // written only when it carries phase information. The alignment's
            // in-gap discovery also calls homozygous variants, and emitting
            // those would add ~12,850 hom records to chr20 -- a 23% larger VCF
            // whose extra content is alignment-discovered calls appearing ONLY
            // inside recovery windows, which is a biased subset of the genome.
            // The arm's output is the catalog's sites plus what recovery
            // phased, not a variant call set.
            if (ci < graph_chunk.site_ids.size() && graph_chunk.site_ids[ci].empty()) {
                const auto& h = mcand.hap_to_cons_alle;
                const int n_alleles = static_cast<int>(mcand.counts.alle_covs.size());
                // Both haplotype consensus alleles must be VALID INDICES for this
                // site's allele set as merged, and they must differ. A merged site
                // is collapsed to biallelic (alle_covs = {ref_cov, alt_cov}), so a
                // consensus index of 2 survives from a wider alignment
                // representation and has no allele here; such sites were emitted
                // 1|1 -- 186 of them on chr20, carrying hap_to_cons_alle (2,1) or
                // (1,2).
                if (h.size() < 3 || h[1] < 0 || h[2] < 0 || h[1] >= n_alleles ||
                    h[2] >= n_alleles || h[1] == h[2])
                    continue;
            }

            auto tid_it = contig_to_tid.find(meta.chrom);
            if (tid_it == contig_to_tid.end()) continue;
            const int fai_tid = tid_it->second;

            const std::vector<int>* orig_idx =
                (ci < graph_chunk.site_allele_orig_idx.size())
                    ? &graph_chunk.site_allele_orig_idx[ci]
                    : nullptr;

            const std::vector<int>& alle_covs = mcand.counts.alle_covs;
            const int n_new_alleles = static_cast<int>(alle_covs.size());

            for (int new_a = 1; new_a < n_new_alleles; ++new_a) {
                // Map surviving allele index back to original walk index in GraphSite.
                int orig_walk_idx = new_a;
                if (orig_idx != nullptr && new_a < static_cast<int>(orig_idx->size())) {
                    orig_walk_idx = (*orig_idx)[new_a];
                }

                // orig_walk_idx: 0 = ref walk, 1 = first alt walk → meta.alts[0], etc.
                const int alt_idx = orig_walk_idx - 1;
                if (alt_idx < 0 || alt_idx >= static_cast<int>(meta.alts.size())) continue;
                const std::string& raw_alt = meta.alts[static_cast<size_t>(alt_idx)];
                if (raw_alt.empty() || raw_alt == "*") continue;

                if (meta.ref.empty()) continue;

                // Normalize to minimal VCF form before deriving the key. Catalog
                // alleles carry flanking repeat context on both sides; an equal-
                // length SNP padded with repeat bases (e.g. the AGGG array at
                // chr20:49031440) otherwise falls through to the MNP branch with a
                // misleading multi-bp ref/alt and the wrong POS, and the same SNP
                // is emitted once per overlapping snarl. Trimming the shared suffix
                // then prefix yields the BAM pipeline's canonical representation so
                // keys match and duplicates collapse. (apply_graph_noise_filter
                // already trims for the noise check; this trims the key itself.)
                std::string ref_seq = meta.ref;
                std::string alt_seq = raw_alt;
                hts_pos_t var_pos = meta.pos;
                trim_to_minimal_vcf(var_pos, ref_seq, alt_seq);

                CandidateVariant cand;
                cand.key.tid = fai_tid;

                // Derive VariantKey from VCF-anchored alleles matching BAM-path convention:
                //   SNP : key.pos = site.pos, key.alt = alt base, ref_len = 1
                //   INS : key.pos = site.pos+1 (after anchor), key.alt = inserted bases
                //   DEL : key.pos = site.pos+1 (first deleted base), key.alt = deleted bases
                if (ref_seq.size() == 1 && alt_seq.size() == 1) {
                    // SNP
                    cand.key.type = VariantType::Snp;
                    cand.key.pos = var_pos;
                    cand.key.alt = alt_seq;
                    cand.key.ref_len = 1;
                } else if (alt_seq.size() > ref_seq.size() && ref_seq[0] == alt_seq[0]) {
                    // Left-anchored insertion: strip shared prefix; pos after anchor.
                    cand.key.type = VariantType::Insertion;
                    cand.key.pos = var_pos + 1;
                    cand.key.alt = alt_seq.substr(ref_seq.size());
                    cand.key.ref_len = 0;
                } else if (ref_seq.size() > alt_seq.size() && ref_seq[0] == alt_seq[0]) {
                    // Left-anchored deletion: empty alt matches BAM-path convention.
                    cand.key.type = VariantType::Deletion;
                    cand.key.pos = var_pos + static_cast<hts_pos_t>(alt_seq.size());
                    cand.key.alt = "";
                    cand.key.ref_len = static_cast<int>(ref_seq.size() - alt_seq.size());
                } else {
                    // Complex / MNP: no shared anchor base — classify by net length change.
                    cand.key.pos = var_pos;
                    cand.key.alt = alt_seq;
                    cand.key.ref_len = static_cast<int>(ref_seq.size());
                    if (alt_seq.size() > ref_seq.size()) {
                        cand.key.type = VariantType::Insertion;
                    } else if (ref_seq.size() > alt_seq.size()) {
                        cand.key.type = VariantType::Deletion;
                    } else {
                        cand.key.type = VariantType::Snp;  // MNP (equal length)
                    }
                }

                // Heterozygous between two ALTERNATE alleles: the reference is not
                // one of this site's haplotypes, so every quantity defined against
                // it is degenerate -- ref_cov is 0, and the allele fraction
                // alt/(ref+alt) is 1.0 for BOTH alleles, which reads as a
                // homozygous alt and then as LOW_AF. The site's own coverage is the
                // denominator that means something: at 4,785,719 that turns two
                // fractions of 1.00 into 0.34 and 0.66.
                const auto& hcons = mcand.hap_to_cons_alle;
                const bool het_by_consensus = hcons.size() > 2 && hcons[1] >= 0 &&
                                              hcons[2] >= 0 && hcons[1] != hcons[2];
                const int ref_cov = alle_covs.empty() ? 0 : alle_covs[0];
                const int alt_cov = alle_covs[static_cast<size_t>(new_a)];
                int total_cov = ref_cov + alt_cov;
                if (het_by_consensus) {
                    total_cov = 0;
                    for (const int c : alle_covs) total_cov += c;
                }
                cand.counts.ref_cov = ref_cov;
                cand.counts.alt_cov = alt_cov;
                cand.counts.total_cov = total_cov;
                // Preserve the strand split computed during biallelic decomposition
                // (graph_bam_adapter.cpp). Candidates are biallelic here, so the
                // chunk candidate's forward/reverse fields already correspond to
                // this ref/alt pair; copying them keeps REVERSE counts non-zero.
                cand.counts.forward_ref = mcand.counts.forward_ref;
                cand.counts.reverse_ref = mcand.counts.reverse_ref;
                cand.counts.forward_alt = mcand.counts.forward_alt;
                cand.counts.reverse_alt = mcand.counts.reverse_alt;
                cand.counts.allele_fraction =
                    total_cov > 0 ? static_cast<double>(alt_cov) / total_cov : 0.0;
                cand.counts.n_uniq_alles = 2;
                cand.counts.alle_covs = {ref_cov, alt_cov};

                // A site whose two haplotype consensus alleles DIFFER is
                // heterozygous even when no read carries the reference: a 1|2
                // site has ref_cov == 0 by construction. Deciding hom from depth
                // alone discards exactly the sites that bridge a gap -- measured
                // in chr20:4,766,928-4,792,960, where 4,785,719 ('ATTTT' at 22
                // reads against a pure 25 bp deletion at 43) and 4,791,668 (16 T
                // at 34 against 17 T at 24) are the two heterozygotes the gap
                // needs and both were called CleanHom here.
                const bool is_hom_alt = (ref_cov == 0 && alt_cov >= opts.min_alt_depth) &&
                                        !het_by_consensus;
                if (alt_cov < opts.min_alt_depth || total_cov < opts.min_depth) {
                    cand.counts.category = VariantCategory::LowCoverage;
                    cand.counts.candvarcate_initial = VariantCategory::LowCoverage;
                    cand.lcd_var_i_to_cate = kLongcalldLowCovVar;
                } else if (is_hom_alt) {
                    cand.counts.category = VariantCategory::CleanHom;
                    cand.counts.candvarcate_initial = VariantCategory::CleanHom;
                    cand.lcd_var_i_to_cate = kCandCleanHom;
                } else if (cand.counts.allele_fraction < opts.min_af ||
                           cand.counts.allele_fraction > opts.max_af) {
                    cand.counts.category = VariantCategory::LowAlleleFraction;
                    cand.counts.candvarcate_initial = VariantCategory::LowAlleleFraction;
                    cand.lcd_var_i_to_cate = kLongcalldLowAfVar;
                } else if (cand.key.type == VariantType::Snp) {
                    cand.counts.category = VariantCategory::CleanHetSnp;
                    cand.counts.candvarcate_initial = VariantCategory::CleanHetSnp;
                    cand.lcd_var_i_to_cate = kCandCleanHetSnp;
                } else if (mcand.counts.category == VariantCategory::RepeatHetIndel) {
                    // Preserve the noise-filter demotion (apply_graph_noise_filter).
                    // Re-classifying from scratch would re-promote homopolymer/STR
                    // indels to CleanHetIndel, letting them pass the germline output
                    // gate as phased het indels even though they were excluded from
                    // k-means (hence never properly phased).
                    cand.counts.category = VariantCategory::RepeatHetIndel;
                    cand.counts.candvarcate_initial = VariantCategory::RepeatHetIndel;
                    cand.lcd_var_i_to_cate = kLongcalldRepHetVar;
                } else {
                    cand.counts.category = VariantCategory::CleanHetIndel;
                    cand.counts.candvarcate_initial = VariantCategory::CleanHetIndel;
                    cand.lcd_var_i_to_cate = kCandCleanHetIndel;
                }

                // Translate multi-allelic hap_to_cons_alle to biallelic for this alt.
                // Homozygous alt: both haplotypes carry the alt allele, no phase set.
                const int hap1 = mcand.hap_to_cons_alle[1];
                const int hap2 = mcand.hap_to_cons_alle[2];
                cand.hap_to_cons_alle[0] = -1;
                if (is_hom_alt) {
                    cand.hap_to_cons_alle[1] = 1;
                    cand.hap_to_cons_alle[2] = 1;
                    cand.hap_alt = 1;
                    cand.hap_ref = 1;
                    // Homozygous candidates have no phase block. Keep the BAM
                    // candidate convention (0) rather than the read sentinel
                    // (-1), so both pipelines serialize the same state.
                    cand.phase_set = kUnsetCandidatePhaseSet;
                } else {
                    cand.hap_to_cons_alle[1] = (hap1 == new_a) ? 1 : 0;
                    cand.hap_to_cons_alle[2] = (hap2 == new_a) ? 1 : 0;
                    cand.hap_alt = cand.hap_to_cons_alle[1];
                    cand.hap_ref = cand.hap_to_cons_alle[2];
                    cand.phase_set = mcand.phase_set;
                }
                cand.alt_ref_base = 4;  // use FASTA anchor (BAM-path default)
                cand.lcd_make_variants_region_pass = true;

                // A site merged in by the in-chunk recovery is written only when
                // this writer's own classification calls it a het. Two reasons.
                // The reclassification above sets CleanHom whenever ref_cov == 0
                // (is_hom_alt, line ~196), and an alignment candidate merged from
                // a gap often carries ref_cov = 0 -- measured, ref_cov 0 against
                // alt_cov 57 -- so the depth synthesis manufactures hom calls:
                // 408 extra 1|1 records on chr20 against 125 in the default. And
                // a hom site carries no phase information in any case, so a VCF
                // whose contract is the catalog's sites plus what recovery phased
                // has no reason to gain hom calls that appear only inside
                // recovery windows.
                if (ci < graph_chunk.site_ids.size() && graph_chunk.site_ids[ci].empty()) {
                    const VariantCategory c = cand.counts.category;
                    if (c != VariantCategory::CleanHetSnp &&
                        c != VariantCategory::CleanHetIndel &&
                        c != VariantCategory::NoisyCandHet)
                        continue;
                }

                result.push_back(std::move(cand));
            }
        }
    }

    std::stable_sort(result.begin(), result.end(),
                     [](const CandidateVariant& a, const CandidateVariant& b) {
                         return exact_comp_cand_var(&a, &b) < 0;
                     });

    // Collapse duplicate variants. Overlapping/nested snarls in the catalog can
    // emit the same physical variant from several sites; after minimal-VCF
    // normalization these share an identity and sort adjacently. Keep the
    // highest-coverage copy so the output VCF/TSV carries each variant once.
    // (Phasing already ran upstream per snarl site; this only de-duplicates the
    // emitted table.) If two copies disagree on the phased haplotype, the
    // higher-coverage call wins and the conflict is reported under --verbose.
    if (!result.empty()) {
        size_t write = 0;
        int collapsed = 0;
        int hap_conflicts = 0;
        for (size_t read = 1; read < result.size(); ++read) {
            if (exact_comp_cand_var(&result[write], &result[read]) == 0) {
                ++collapsed;
                if (result[write].phase_set == result[read].phase_set &&
                    result[write].hap_alt != result[read].hap_alt) {
                    ++hap_conflicts;
                }
                if (result[read].counts.total_cov > result[write].counts.total_cov) {
                    result[write] = std::move(result[read]);
                }
            } else {
                ++write;
                if (write != read) result[write] = std::move(result[read]);
            }
        }
        result.resize(write + 1);
        if (opts.verbose && collapsed > 0) {
            std::cerr << "graph: collapsed " << collapsed
                      << " duplicate variant record(s) from overlapping snarls ("
                      << hap_conflicts << " with conflicting haplotype calls)\n";
        }
    }

    make_colocated_alleles_complementary(result, opts.min_alt_depth);
    drop_conflicting_haplotype_alleles(result);

    return result;
}

// A graph SNP can look heterozygous when the physical molecules carry its ALT
// base on one haplotype and a deletion over its REF base on the other. The
// graph's REF/ALT labels then give the wrong read orientation. Check only
// graph-clean, biallelic SNPs and require decisive physical evidence before
// excluding one from the graph solve. This uses the existing per-thread BAM
// handles and scans each chunk once rather than seeking for every site.
static bool exclude_ref_absent_graph_snps(
        GraphSiteCatalog& catalog, const GraphChunkBuildResult& built,
        const RegionChunk& region, const std::string& contig,
        WorkerContext& context, const Options& opts) {
    struct SiteEvidence {
        hts_pos_t pos;
        size_t catalog_index;
        char ref;
        char alt;
        int ref_count = 0;
        int alt_count = 0;
        int deletion_count = 0;
        int other_count = 0;
    };
    std::unordered_map<std::string, size_t> catalog_index;
    catalog_index.reserve(catalog.sites.size());
    for (size_t i = 0; i < catalog.sites.size(); ++i)
        catalog_index.emplace(graph_site_key_str(catalog.sites[i]), i);

    std::vector<SiteEvidence> sites;
    for (size_t i = 0; i < built.chunk.candidates.size(); ++i) {
        if (built.chunk.candidates[i].counts.category != VariantCategory::CleanHetSnp)
            continue;
        const auto found = catalog_index.find(built.site_ids[i]);
        if (found == catalog_index.end()) continue;
        const GraphSite& site = catalog.sites[found->second];
        if (site.ref.size() != 1 || site.alts.size() != 1 ||
            site.alts[0].size() != 1 || site.ref == site.alts[0])
            continue;
        sites.push_back({site.pos, found->second, site.ref[0], site.alts[0][0]});
    }
    if (sites.empty()) return false;
    std::sort(sites.begin(), sites.end(),
              [](const SiteEvidence& a, const SiteEvidence& b) {
                  return a.pos < b.pos;
              });

    // Graph GAF phasing admits MAPQ 5, but rejecting a graph allele needs
    // the BAM pipeline's high-confidence alignment floor.
    constexpr int kMinPhysicalValidationMapq = 30;
    const int min_mapq = std::max(opts.min_mapq, kMinPhysicalValidationMapq);
    for (size_t input = 0; input < context.bams.size(); ++input) {
        const int tid = sam_hdr_name2tid(context.headers[input].get(), contig.c_str());
        if (tid < 0) continue;
        std::unique_ptr<hts_itr_t, decltype(&hts_itr_destroy)> itr(
            sam_itr_queryi(context.indexes[input].get(), tid,
                           region.beg - 1, region.end), hts_itr_destroy);
        if (!itr) throw std::runtime_error("failed to query BAM for graph SNP validation: " + contig);
        std::unique_ptr<bam1_t, decltype(&bam_destroy1)> read(bam_init1(), bam_destroy1);
        if (!read) throw std::runtime_error("failed to allocate BAM record for graph SNP validation");
        int status = 0;
        while ((status = sam_itr_next(context.bams[input]->get(), itr.get(), read.get())) >= 0) {
            const bam1_core_t& core = read->core;
            if ((core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) ||
                (!opts.include_filtered && (core.flag & (BAM_FQCFAIL | BAM_FDUP))) ||
                core.qual < min_mapq || core.qual == 255)
                continue;
            hts_pos_t ref_pos = core.pos + 1;
            int query_pos = 0;
            const uint32_t* cigar = bam_get_cigar(read.get());
            const uint8_t* sequence = bam_get_seq(read.get());
            const uint8_t* quality = bam_get_qual(read.get());
            for (uint32_t op_i = 0; op_i < core.n_cigar; ++op_i) {
                const int op = bam_cigar_op(cigar[op_i]);
                const int length = bam_cigar_oplen(cigar[op_i]);
                const int consumed = bam_cigar_type(op);
                if (consumed & 2) {
                    auto it = std::lower_bound(
                        sites.begin(), sites.end(), ref_pos,
                        [](const SiteEvidence& site, hts_pos_t pos) {
                            return site.pos < pos;
                        });
                    for (; it != sites.end() && it->pos < ref_pos + length; ++it) {
                        if (op == BAM_CDEL) {
                            ++it->deletion_count;
                        } else if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
                            const int qi = query_pos + static_cast<int>(it->pos - ref_pos);
                            if (qi < 0 || qi >= core.l_qseq || quality[qi] < opts.min_bq)
                                continue;
                            const char base = seq_nt16_str[bam_seqi(sequence, qi)];
                            if (base == it->ref) ++it->ref_count;
                            else if (base == it->alt) ++it->alt_count;
                            else ++it->other_count;
                        }
                    }
                    ref_pos += length;
                }
                if (consumed & 1) query_pos += length;
            }
        }
        if (status < -1)
            throw std::runtime_error("failed to read BAM for graph SNP validation: " + contig);
    }

    // A true REF/ALT heterozygote has probability 2^-n of yielding zero REF
    // bases in n callable observations. Correct for all tested graph SNPs.
    // Also require a substantial deletion allele; otherwise a homozygous ALT
    // site with incidental indel errors could be excluded here.
    constexpr int kMinDeletionObservations = 10;
    constexpr double kMinDeletionFraction = 0.2;
    constexpr double kFamilywiseError = 0.01;
    bool excluded = false;
    for (const SiteEvidence& site : sites) {
        const int callable = site.ref_count + site.alt_count;
        const int total = callable + site.deletion_count;
        if (site.ref_count != 0 || site.other_count != 0 ||
            site.deletion_count < kMinDeletionObservations || total == 0 ||
            static_cast<double>(site.deletion_count) / total < kMinDeletionFraction ||
            std::ldexp(1.0, -callable) * sites.size() > kFamilywiseError)
            continue;
        GraphSite& graph_site = catalog.sites[site.catalog_index];
        graph_site.eligible = false;
        graph_site.skip_reason = "bam_alt_deletion_no_ref";
        excluded = true;
        if (opts.verbose)
            std::cerr << "graph: excluded REF-absent SNP " << contig << ':' << site.pos
                      << " (REF " << site.ref_count << ", ALT " << site.alt_count
                      << ", deletion " << site.deletion_count << ")\n";
    }
    return excluded;
}

// A graph chunk boundary can cut through a BAM-supported seam without leaving
// any phased read in both graph chunks. Re-solve only short, unjoined boundary
// gaps and transfer the result when exact clean SNPs on each original block
// agree with one physically validated bridge in the local replay.
template <class SolveReplay>
static void bridge_graph_chunk_boundaries(
        std::vector<PhasingChunk>& chunks, SolveReplay&& solve_replay) {
    constexpr hts_pos_t kReplayFlank = 50000;
    constexpr hts_pos_t kMinReplayGap = 10000;
    constexpr hts_pos_t kMaxReplayGap = 50000;
    constexpr size_t kMinSharedSnpsPerSide = 2;
    const auto anchor = [](const CandidateVariant& site) {
        return !site.bam_injected && site.phase_set > 0 &&
            site.key.type == VariantType::Snp &&
            site.counts.category == VariantCategory::CleanHetSnp &&
            site.hap_to_cons_alle[1] >= 0 &&
            site.hap_to_cons_alle[1] <= 1 &&
            site.hap_to_cons_alle[2] >= 0 &&
            site.hap_to_cons_alle[2] <= 1 &&
            site.hap_to_cons_alle[1] != site.hap_to_cons_alle[2];
    };
    const auto same_key = [](const VariantKey& a, const VariantKey& b) {
        return a.pos == b.pos && a.type == b.type &&
            a.ref_len == b.ref_len && a.alt == b.alt;
    };
    for (size_t i = 1; i < chunks.size(); ++i) {
        PhasingChunk& left = chunks[i - 1];
        PhasingChunk& right = chunks[i];
        if (left.region.tid != right.region.tid ||
            left.region.end + 1 != right.region.beg)
            continue;
        const CandidateVariant* left_anchor = nullptr;
        const CandidateVariant* right_anchor = nullptr;
        for (const CandidateVariant& site : left.candidates)
            if (anchor(site) &&
                (left_anchor == nullptr ||
                 site.key.pos > left_anchor->key.pos))
                left_anchor = &site;
        for (const CandidateVariant& site : right.candidates)
            if (anchor(site) &&
                (right_anchor == nullptr ||
                 site.key.pos < right_anchor->key.pos))
                right_anchor = &site;
        if (left_anchor == nullptr || right_anchor == nullptr ||
            left_anchor->phase_set == right_anchor->phase_set)
            continue;
        const hts_pos_t gap = right_anchor->key.pos - left_anchor->key.pos;
        if (gap < kMinReplayGap || gap > kMaxReplayGap)
            continue;
        RegionChunk region;
        region.tid = left.region.tid;
        region.beg = std::max(left.region.beg,
                              left.region.end - kReplayFlank + 1);
        region.end = std::min(right.region.end,
                              left.region.end + kReplayFlank);
        region.chunk_id = left.region.chunk_id;
        region.reg_chunk_i = left.region.reg_chunk_i;
        if (left_anchor->key.pos < region.beg ||
            right_anchor->key.pos > region.end)
            continue;
        GraphChunkBuildResult replay = solve_replay(region);
        const auto find_replay = [&](const VariantKey& key)
                -> std::optional<size_t> {
            for (size_t ci = 0; ci < replay.chunk.candidates.size(); ++ci)
                if (anchor(replay.chunk.candidates[ci]) &&
                    same_key(replay.chunk.candidates[ci].key, key))
                    return ci;
            return std::nullopt;
        };
        const auto replay_left = find_replay(left_anchor->key);
        const auto replay_right = find_replay(right_anchor->key);
        if (!replay_left || !replay_right)
            continue;
        const hts_pos_t replay_ps =
            replay.chunk.candidates[*replay_left].phase_set;
        if (replay_ps != replay.chunk.candidates[*replay_right].phase_set)
            continue;
        const auto original_graph_ps = [&](size_t candidate_i)
                -> std::optional<hts_pos_t> {
            std::optional<hts_pos_t> source;
            for (const RecoverySourceSite& site : replay.recovery_source_sites) {
                if (site.candidate_index != candidate_i ||
                    site.graph_phase_set <= 0)
                    continue;
                if (source && *source != site.graph_phase_set)
                    return std::nullopt;
                source = site.graph_phase_set;
            }
            return source;
        };
        const auto source_left = original_graph_ps(*replay_left);
        const auto source_right = original_graph_ps(*replay_right);
        // A boundary SNP can be absent from source_sites when recovery creates
        // its phase assignment. The physical bridge still records the original
        // graph phase-set ID, so use it for that side when the exact SNP match
        // and the block-wide parity check below corroborate the replay.
        const hts_pos_t left_ps = left_anchor->phase_set;
        const hts_pos_t right_ps = right_anchor->phase_set;
        bool physical_bridge = false;
        for (const RecoveryPhaseGauge& gauge : replay.recovery_phase_gauges) {
            for (const RecoveryPhysicalSnpBridge& bridge :
                 gauge.physical_snp_bridges) {
                const bool left_matches = source_left
                    ? bridge.left_phase_set == *source_left
                    : bridge.left_phase_set == left_ps;
                const bool right_matches = source_right
                    ? bridge.right_phase_set == *source_right
                    : bridge.right_phase_set == right_ps;
                physical_bridge |= left_matches && right_matches &&
                    bridge.left_phase_set != bridge.right_phase_set;
            }
        }
        if (!physical_bridge)
            continue;
        const auto block_parity = [&](const PhasingChunk& original,
                                      hts_pos_t phase_set)
                -> std::optional<bool> {
            std::optional<bool> parity;
            size_t shared = 0;
            for (const CandidateVariant& site : original.candidates) {
                if (!anchor(site) || site.phase_set != phase_set ||
                    site.key.pos < region.beg || site.key.pos > region.end)
                    continue;
                const auto found = find_replay(site.key);
                if (!found) continue;
                const CandidateVariant& matched =
                    replay.chunk.candidates[*found];
                if (matched.phase_set != replay_ps)
                    return std::nullopt;
                const bool current =
                    site.hap_to_cons_alle[1] != matched.hap_to_cons_alle[1];
                if (parity && *parity != current)
                    return std::nullopt;
                parity = current;
                ++shared;
            }
            return shared >= kMinSharedSnpsPerSide ? parity : std::nullopt;
        };
        const auto left_parity = block_parity(left, left_ps);
        const auto right_parity = block_parity(right, right_ps);
        if (!left_parity || !right_parity)
            continue;
        const bool flip = *left_parity != *right_parity;
        for (size_t j = i; j < chunks.size(); ++j) {
            PhasingChunk& chunk = chunks[j];
            if (chunk.region.tid != right.region.tid) break;
            for (CandidateVariant& site : chunk.candidates) {
                if (site.phase_set != right_ps) continue;
                if (flip)
                    std::swap(site.hap_to_cons_alle[1],
                              site.hap_to_cons_alle[2]);
                site.phase_set = left_ps;
            }
            for (size_t ri = 0; ri < chunk.reads.size(); ++ri) {
                if (ri >= chunk.phase_sets.size() ||
                    chunk.phase_sets[ri] != right_ps)
                    continue;
                if (flip && ri < chunk.haps.size() &&
                    (chunk.haps[ri] == 1 || chunk.haps[ri] == 2))
                    chunk.haps[ri] = 3 - chunk.haps[ri];
                chunk.phase_sets[ri] = left_ps;
            }
        }
    }
}

// Processes one batch of graph chunks in parallel (one thread pool per reg_chunk_i batch,
// mirroring collect_chunk_batch_parallel in collect_pipeline.cpp).
// Each worker queries overlapping reads via the gbz-base FFI (one SQLite connection
// per thread), then build_graph_chunk + assign_hap k-means.
// Peak memory = threads × (reads_per_chunk + sites_per_chunk).
// After all workers join: populate_graph_chunk_overlaps + stitch_chunk_haps.
static std::vector<GraphChunkBuildResult> process_graph_chunk_batch(
    const std::string& sites_vcf,
    const std::vector<RegionChunk>& chunks,
    size_t batch_begin,
    size_t batch_end,
    const bam_hdr_t* header,
    const GraphQueryConfig& qconfig,
    const std::string& ref_sample,
    const std::unordered_map<std::string, std::string>& fai_full_to_suffix,
    const std::unordered_map<std::string, std::string>& chrom_remap,
    const Options& opts,
    const PgbamSidecarData* pgbam_sidecar)
{
    const size_t batch_size = batch_end - batch_begin;
    std::vector<GraphChunkBuildResult> graph_chunks(batch_size);

    const std::string batch_contig = header->target_name[chunks[batch_begin].tid];

    // Pre-compute the query contig (pangenome suffix) once for the batch.
    const std::string batch_query_contig = [&]() -> std::string {
        auto it = fai_full_to_suffix.find(batch_contig);
        return (it != fai_full_to_suffix.end()) ? it->second : batch_contig;
    }();

    // Query the genome-coordinate range of the reference path so we can
    // skip or clamp chunks that fall outside the GBZ subgraph.  This is
    // done once on the main thread before spawning workers.
    PathRange path_range;
    {
        char* err = nullptr;
        void* gbz_tmp = pgphase_gbz_open(qconfig.gbz_db.c_str(), &err);
        if (gbz_tmp) {
            path_range = query_gbz_path_range(gbz_tmp, ref_sample, batch_query_contig);
            pgphase_gbz_close(gbz_tmp);
            if (opts.verbose >= 2 && path_range.valid)
                std::cerr << "GBZ path range for " << ref_sample << "#" << batch_query_contig
                          << ": [" << path_range.start << ", " << path_range.end << ")\n";
        } else {
            if (err) pgphase_gbz_free_string(err);
        }
    }

    const size_t worker_count = std::min<size_t>(static_cast<size_t>(opts.threads), batch_size);
    std::atomic<size_t> next_offset{0};
    std::exception_ptr first_error;
    std::mutex error_mutex;
    std::vector<std::thread> workers;
    workers.reserve(worker_count);

    for (size_t w = 0; w < worker_count; ++w) {
        workers.emplace_back([&]() {
            char* err = nullptr;
            void* gbz_h = pgphase_gbz_open(qconfig.gbz_db.c_str(), &err);
            if (!gbz_h) {
                std::string msg = err ? std::string(err) : "unknown error";
                if (err) pgphase_gbz_free_string(err);
                std::lock_guard<std::mutex> lock(error_mutex);
                if (!first_error)
                    first_error = std::make_exception_ptr(
                        std::runtime_error("failed to open GBZ: " + msg));
                return;
            }
            void* gaf_h = pgphase_gaf_open(qconfig.gaf_db.c_str(), &err);
            if (!gaf_h) {
                std::string msg = err ? std::string(err) : "unknown error";
                if (err) pgphase_gbz_free_string(err);
                pgphase_gbz_close(gbz_h);
                std::lock_guard<std::mutex> lock(error_mutex);
                if (!first_error)
                    first_error = std::make_exception_ptr(
                        std::runtime_error("failed to open GAF-base: " + msg));
                return;
            }
            if (pgphase_gbz_gaf_validate(gbz_h, gaf_h, &err) != 0) {
                std::string msg = err ? std::string(err) : "unknown error";
                if (err) pgphase_gbz_free_string(err);
                pgphase_gaf_close(gaf_h);
                pgphase_gbz_close(gbz_h);
                std::lock_guard<std::mutex> lock(error_mutex);
                if (!first_error)
                    first_error = std::make_exception_ptr(
                        std::runtime_error("GBZ/GAF-base incompatible: " + msg));
                return;
            }
            struct HandleCleanup {
                void* gbz; void* gaf;
                ~HandleCleanup() {
                    if (gaf) pgphase_gaf_close(gaf);
                    if (gbz) pgphase_gbz_close(gbz);
                }
            } cleanup{gbz_h, gaf_h};

            try {
                // One sites VCF handle per thread.
                SitesVcfHandle sites_handle(sites_vcf);

                // Per-thread reference index for noise detection.
                std::unique_ptr<faidx_t, FaiDeleter> thread_fai(
                    load_reference_index(opts.ref_fasta));

                // Per-thread recovery context. htslib handles are not shareable,
                // which is why the parent pipeline opens a set per worker, and it
                // is the only reason recovery could not already live here. Built
                // once per thread rather than per chunk.
                std::unique_ptr<WorkerContext> thread_recovery_ctx;
                if (!opts.bam_files.empty())
                    thread_recovery_ctx = std::make_unique<WorkerContext>(opts);

                while (true) {
                    const size_t offset = next_offset.fetch_add(1);
                    if (offset >= batch_size) break;
                    const RegionChunk& region = chunks[batch_begin + offset];

                    // Load sites for this chunk's region via tabix.
                    GraphSiteCatalog chunk_catalog = load_sites_for_region(
                        sites_handle, batch_contig, region.beg, region.end);
                    // Normalize contig names to match FAI convention.
                    for (GraphSite& s : chunk_catalog.sites) {
                        auto it = chrom_remap.find(s.chrom);
                        if (it != chrom_remap.end()) s.chrom = it->second;
                        if (!s.ref_contig.empty()) {
                            auto it2 = chrom_remap.find(s.ref_contig);
                            if (it2 != chrom_remap.end()) s.ref_contig = it2->second;
                        }
                    }

                    GraphSiteCatalogView chunk_view = chunk_catalog.view_all();

                    std::vector<GraphReadAllele> chunk_rows;
                    if (!chunk_view.empty()) {
                        hts_pos_t q_beg = region.beg - 1;
                        hts_pos_t q_end = region.end;
                        if (path_range.valid) {
                            q_beg = std::max(q_beg, path_range.start);
                            q_end = std::min(q_end, path_range.end);
                        }
                        if (q_beg < q_end) {
                            chunk_rows = query_gbz_interval_gaf_ffi(
                                gbz_h, gaf_h, ref_sample, batch_query_contig,
                                q_beg, q_end,
                                chunk_view, qconfig.min_mapq);
                        }
                    }

                    graph_chunks[offset] = build_graph_chunk(
                        chunk_view,
                        chunk_rows,
                        batch_contig,
                        region.beg - 1,
                        region.end,
                        region.chunk_id,
                        opts);
                    if (thread_recovery_ctx != nullptr &&
                        exclude_ref_absent_graph_snps(chunk_catalog, graph_chunks[offset],
                                                      region, batch_contig,
                                                      *thread_recovery_ctx, opts)) {
                        graph_chunks[offset] = build_graph_chunk(
                            chunk_view, chunk_rows, batch_contig, region.beg - 1,
                            region.end, region.chunk_id, opts);
                    }

                    // Noise filter: fetch reference slice and reclassify
                    // indels in homopolymer/repeat/low-complexity contexts.
                    {
                        hts_pos_t ref_len = 0;
                        char* ref_raw = faidx_fetch_seq64(
                            thread_fai.get(), batch_contig.c_str(),
                            region.beg - 1, region.end - 1, &ref_len);
                        if (ref_raw && ref_len > 0) {
                            std::string ref_slice(ref_raw, ref_raw + ref_len);
                            std::free(ref_raw);
                            apply_graph_noise_filter(
                                graph_chunks[offset], ref_slice,
                                region.beg, region.beg + ref_len - 1,
                                opts.noisy_reg_max_xgaps);
                        } else {
                            std::free(ref_raw);
                        }
                    }

                    // A repeat-context indel may earn its way back in before the solve.
                    promote_link_supported_repeat_indels(graph_chunks[offset], opts);

                    assign_hap_based_on_germline_het_vars_kmeans(
                        graph_chunks[offset].chunk, opts, kCandGermlineClean);

                    // Recover and stitch local BAM blocks before cross-chunk
                    // overlaps are computed. Each worker owns one chunk, so the
                    // merge and its index rebuild require no synchronization.
                    if (thread_recovery_ctx != nullptr)
                        run_in_chunk_recovery(graph_chunks[offset], opts, *thread_recovery_ctx,
                                              batch_contig.c_str());
                    if (thread_recovery_ctx != nullptr) {
                        recover_independent_bam_read_blocks_in_place(
                            graph_chunks[offset], opts, *thread_recovery_ctx,
                            batch_contig.c_str());
                    }
                }
            } catch (...) {
                std::lock_guard<std::mutex> lock(error_mutex);
                if (!first_error) first_error = std::current_exception();
            }
        });
    }
    for (std::thread& w : workers) w.join();
    if (first_error) std::rethrow_exception(first_error);

    populate_graph_chunk_overlaps(graph_chunks);

    std::vector<PhasingChunk> phasing_chunks;
    phasing_chunks.reserve(batch_size);
    for (GraphChunkBuildResult& gc : graph_chunks)
        phasing_chunks.push_back(std::move(gc.chunk));
    stitch_chunk_haps(phasing_chunks, &opts, pgbam_sidecar);
    if (batch_size > 1 && !opts.bam_files.empty() &&
        pgbam_sidecar == nullptr) {
        Options replay_opts = opts;
        replay_opts.threads = 1;
        replay_opts.phase_matrix_dump_prefix.clear();
        bridge_graph_chunk_boundaries(phasing_chunks,
            [&](const RegionChunk& region) {
                const std::vector<RegionChunk> replay_region{region};
                auto result = process_graph_chunk_batch(
                    sites_vcf, replay_region, 0, 1, header, qconfig,
                    ref_sample, fai_full_to_suffix, chrom_remap,
                    replay_opts, nullptr);
                return std::move(result.front());
            });
    }
    for (size_t i = 0; i < batch_size; ++i) {
        graph_chunks[i].chunk = std::move(phasing_chunks[i]);
        rescue_unphased_graph_reads(
            graph_chunks[i].chunk, graph_chunks[i].recovery_windows);
        apply_independent_bam_read_blocks(graph_chunks[i].chunk);
    }

    return graph_chunks;
}

// Per-chunk tabix queries on an indexed GAF file.  Each worker seeks directly
// to the overlapping region — only the relevant reads are decompressed and
// parsed, making this efficient even for very large (100+ GB) GAF files.
static std::vector<GraphChunkBuildResult> process_graph_chunk_batch_indexed_gaf(
    const std::string& sites_vcf,
    const std::vector<RegionChunk>& chunks,
    size_t batch_begin,
    size_t batch_end,
    const bam_hdr_t* header,
    const std::string& gaf_file,
    int min_mapq,
    const std::unordered_map<std::string, std::string>& fai_full_to_suffix,
    const std::unordered_map<std::string, std::string>& chrom_remap,
    const Options& opts,
    const PgbamSidecarData* pgbam_sidecar)
{
    const size_t batch_size = batch_end - batch_begin;
    std::vector<GraphChunkBuildResult> graph_chunks(batch_size);

    const std::string batch_contig_gaf = header->target_name[chunks[batch_begin].tid];
    const std::string batch_query_contig_gaf = [&]() -> std::string {
        auto it = fai_full_to_suffix.find(batch_contig_gaf);
        return (it != fai_full_to_suffix.end()) ? it->second : batch_contig_gaf;
    }();

    const size_t worker_count = std::min<size_t>(static_cast<size_t>(opts.threads), batch_size);
    std::atomic<size_t> next_offset{0};
    std::exception_ptr first_error;
    std::mutex error_mutex;
    std::vector<std::thread> workers;
    workers.reserve(worker_count);

    for (size_t w = 0; w < worker_count; ++w) {
        workers.emplace_back([&]() {
            try {
                IndexedGafHandle gaf_handle(gaf_file);
                SitesVcfHandle sites_handle(sites_vcf);

                // Per-thread reference index for noise detection.
                std::unique_ptr<faidx_t, FaiDeleter> thread_fai(
                    load_reference_index(opts.ref_fasta));

                // Per-thread recovery context. htslib handles are not shareable,
                // which is why the parent pipeline opens a set per worker, and it
                // is the only reason recovery could not already live here. Built
                // once per thread rather than per chunk.
                std::unique_ptr<WorkerContext> thread_recovery_ctx;
                if (!opts.bam_files.empty())
                    thread_recovery_ctx = std::make_unique<WorkerContext>(opts);

                while (true) {
                    const size_t offset = next_offset.fetch_add(1);
                    if (offset >= batch_size) break;
                    const RegionChunk& region = chunks[batch_begin + offset];

                    GraphSiteCatalog chunk_catalog = load_sites_for_region(
                        sites_handle, batch_contig_gaf, region.beg, region.end);
                    for (GraphSite& s : chunk_catalog.sites) {
                        auto it = chrom_remap.find(s.chrom);
                        if (it != chrom_remap.end()) s.chrom = it->second;
                        if (!s.ref_contig.empty()) {
                            auto it2 = chrom_remap.find(s.ref_contig);
                            if (it2 != chrom_remap.end()) s.ref_contig = it2->second;
                        }
                    }

                    GraphSiteCatalogView chunk_view = chunk_catalog.view_all();

                    std::vector<GraphReadAllele> chunk_rows;
                    if (!chunk_view.empty()) {
                        const hts_pos_t pad = static_cast<hts_pos_t>(opts.gaf_pad);
                        chunk_rows = scan_indexed_gaf_chunk(
                            gaf_handle, batch_query_contig_gaf,
                            std::max<hts_pos_t>(0, region.beg - 1 - pad), region.end + pad,
                            chunk_view, min_mapq);
                    }

                    graph_chunks[offset] = build_graph_chunk(
                        chunk_view,
                        chunk_rows,
                        batch_contig_gaf,
                        region.beg - 1,
                        region.end,
                        region.chunk_id,
                        opts);
                    if (thread_recovery_ctx != nullptr &&
                        exclude_ref_absent_graph_snps(chunk_catalog, graph_chunks[offset],
                                                      region, batch_contig_gaf,
                                                      *thread_recovery_ctx, opts)) {
                        graph_chunks[offset] = build_graph_chunk(
                            chunk_view, chunk_rows, batch_contig_gaf, region.beg - 1,
                            region.end, region.chunk_id, opts);
                    }

                    // Noise filter: fetch reference slice and reclassify
                    // indels in homopolymer/repeat/low-complexity contexts.
                    {
                        hts_pos_t ref_len = 0;
                        char* ref_raw = faidx_fetch_seq64(
                            thread_fai.get(), batch_contig_gaf.c_str(),
                            region.beg - 1, region.end - 1, &ref_len);
                        if (ref_raw && ref_len > 0) {
                            std::string ref_slice(ref_raw, ref_raw + ref_len);
                            std::free(ref_raw);
                            apply_graph_noise_filter(
                                graph_chunks[offset], ref_slice,
                                region.beg, region.beg + ref_len - 1,
                                opts.noisy_reg_max_xgaps);
                        } else {
                            std::free(ref_raw);
                        }
                    }

                    // A repeat-context indel may earn its way back in before the solve.
                    promote_link_supported_repeat_indels(graph_chunks[offset], opts);

                    assign_hap_based_on_germline_het_vars_kmeans(
                        graph_chunks[offset].chunk, opts, kCandGermlineClean);

                    // Recover and stitch local BAM blocks before cross-chunk
                    // overlaps are computed. Each worker owns one chunk, so the
                    // merge and its index rebuild require no synchronization.
                    if (thread_recovery_ctx != nullptr)
                        run_in_chunk_recovery(graph_chunks[offset], opts, *thread_recovery_ctx,
                                              batch_contig_gaf.c_str());
                    if (thread_recovery_ctx != nullptr) {
                        recover_independent_bam_read_blocks_in_place(
                            graph_chunks[offset], opts, *thread_recovery_ctx,
                            batch_contig_gaf.c_str());
                    }
                }
            } catch (...) {
                std::lock_guard<std::mutex> lock(error_mutex);
                if (!first_error) first_error = std::current_exception();
            }
        });
    }
    for (std::thread& w : workers) w.join();
    if (first_error) std::rethrow_exception(first_error);

    if (opts.verbose >= 1) graph_query_report_match_stats();
    populate_graph_chunk_overlaps(graph_chunks);

    std::vector<PhasingChunk> phasing_chunks;
    phasing_chunks.reserve(batch_size);
    for (GraphChunkBuildResult& gc : graph_chunks)
        phasing_chunks.push_back(std::move(gc.chunk));
    stitch_chunk_haps(phasing_chunks, &opts, pgbam_sidecar);
    if (batch_size > 1 && !opts.bam_files.empty() &&
        pgbam_sidecar == nullptr) {
        Options replay_opts = opts;
        replay_opts.threads = 1;
        replay_opts.phase_matrix_dump_prefix.clear();
        bridge_graph_chunk_boundaries(phasing_chunks,
            [&](const RegionChunk& region) {
                const std::vector<RegionChunk> replay_region{region};
                auto result = process_graph_chunk_batch_indexed_gaf(
                    sites_vcf, replay_region, 0, 1, header, gaf_file,
                    min_mapq, fai_full_to_suffix, chrom_remap,
                    replay_opts, nullptr);
                return std::move(result.front());
            });
    }
    for (size_t i = 0; i < batch_size; ++i) {
        graph_chunks[i].chunk = std::move(phasing_chunks[i]);
        rescue_unphased_graph_reads(
            graph_chunks[i].chunk, graph_chunks[i].recovery_windows);
        apply_independent_bam_read_blocks(graph_chunks[i].chunk);
    }

    return graph_chunks;
}

void run_collect_graph_variation(const Options& opts) {
    const bool use_indexed_gaf = !opts.gaf_file.empty();
    if (!use_indexed_gaf) {
        if (opts.gbz_db.empty())
            throw std::runtime_error("--gbz-db is required for collect-graph-variation without --gaf");
        if (opts.gaf_db.empty())
            throw std::runtime_error("--gaf-db is required for collect-graph-variation without --gaf");
    } else {
        // The --gaf path requires a tabix-indexed, bgzip-compressed GAF with
        // annotated coordinate columns so per-chunk region queries are efficient.
        require_indexed_gaf(opts.gaf_file);
        if (!opts.gbz_db.empty() || !opts.gaf_db.empty()) {
            std::cerr << "Warning: --gaf provided; ignoring --gbz-db/--gaf-db\n";
        }
    }

    // Load optional .pgbam sidecar for fallback chunk stitching.
    std::unique_ptr<PgbamSidecarData> pgbam_sidecar;
    if (!opts.pgbam_file.empty()) {
        pgbam_sidecar = std::make_unique<PgbamSidecarData>(load_pgbam_sidecar(opts.pgbam_file));
        if (opts.verbose >= 1)
            std::cerr << "Loaded pgbam sidecar with "
                      << pgbam_sidecar->set_to_threads.size() << " sets from "
                      << opts.pgbam_file << "\n";
    }

    // 1. Reference FASTA index first — needed to resolve autosome contig names for
    //    region filters, which are then passed to the VCF loader so only sites in the
    //    requested regions are parsed (tabix-assisted for bgzipped + indexed VCFs).
    std::unique_ptr<faidx_t, FaiDeleter> fai(load_reference_index(opts.ref_fasta));
    std::unique_ptr<bam_hdr_t, HeaderDeleter> header(build_synthetic_header(fai.get()));

    // 2. Build a chrom alias map between the FASTA and VCF naming conventions.
    //    Pangenome FASTAs use "SAMPLE#HAP#CHROM" (e.g. "CHM13#0#chr20") while graph VCFs
    //    use the plain reference name ("chr20"), or vice versa.  We inspect the FAI names
    //    and build a bidirectional suffix map so both directions resolve automatically.
    //
    //    fai_suffix_to_full : "chr20" → "CHM13#0#chr20"  (used when VCF is short, FAI is full)
    //    fai_full_to_suffix : "CHM13#0#chr20" → "chr20"  (used when VCF is full, FAI is short)
    std::unordered_map<std::string, std::string> fai_suffix_to_full;
    std::unordered_map<std::string, std::string> fai_full_to_suffix;
    {
        const int nseq = faidx_nseq(fai.get());
        fai_suffix_to_full.reserve(static_cast<size_t>(nseq));
        for (int i = 0; i < nseq; ++i) {
            const std::string full(faidx_iseq(fai.get(), i));
            const size_t h = full.rfind('#');
            if (h != std::string::npos) {
                fai_suffix_to_full.emplace(full.substr(h + 1), full);
                fai_full_to_suffix.emplace(full, full.substr(h + 1));
            }
        }
    }

    // Resolves a contig name to its canonical FAI form (the name that exists in the header).
    // Handles both "chr20" → "CHM13#0#chr20" and "CHM13#0#chr20" → "chr20".
    auto resolve_contig = [&](const std::string& name) -> std::string {
        if (faidx_has_seq(fai.get(), name.c_str())) return name;
        // Short name → full pangenome name ("chr20" → "CHM13#0#chr20")
        auto it = fai_suffix_to_full.find(name);
        if (it != fai_suffix_to_full.end()) return it->second;
        // Full pangenome name → short name ("CHM13#0#chr20" → "chr20")
        auto it2 = fai_full_to_suffix.find(name);
        if (it2 != fai_full_to_suffix.end() && faidx_has_seq(fai.get(), it2->second.c_str()))
            return it2->second;
        return name; // unchanged; add_filter_chunks will throw a clear error
    };

    std::vector<RegionFilter> region_filters;
    for (const std::string& r : opts.regions) {
        RegionFilter f = parse_region(r);
        f.chrom = resolve_contig(f.chrom);
        region_filters.push_back(std::move(f));
    }
    if (!opts.region_file.empty()) {
        auto bed = load_bed_regions(opts.region_file);
        for (RegionFilter& f : bed) f.chrom = resolve_contig(f.chrom);
        region_filters.insert(region_filters.end(), bed.begin(), bed.end());
    }
    if (opts.autosome) {
        // With pangenome FASTAs, "chr1"–"chr22" won't be found directly; resolve_contig
        // maps them to the full name (e.g. "CHM13#0#chr1").
        for (int i = 1; i <= 22; ++i) {
            for (const std::string& candidate :
                     {"chr" + std::to_string(i), std::to_string(i)}) {
                const std::string resolved = resolve_contig(candidate);
                if (faidx_has_seq(fai.get(), resolved.c_str())) {
                    region_filters.push_back(RegionFilter{true, resolved, 1, -1});
                    break;
                }
            }
        }
    }

    // When no region filters were specified, build whole-chromosome filters for every
    // FASTA contig so the VCF loader can use tabix instead of streaming the entire file.
    if (region_filters.empty()) {
        const int nseq = faidx_nseq(fai.get());
        region_filters.reserve(static_cast<size_t>(nseq));
        for (int i = 0; i < nseq; ++i)
            region_filters.push_back(
                RegionFilter{true, std::string(faidx_iseq(fai.get(), i)), 1, -1});
    }

    // 3. Build VCF-name → FAI-name chrom remap for per-chunk site normalization.
    //    Sites are loaded per-chunk via tabix, so no global catalog is needed.
    std::unordered_map<std::string, std::string> chrom_remap;
    {
        // VCF uses short name, FAI uses full ("chr20" → "CHM13#0#chr20")
        for (const auto& kv : fai_suffix_to_full) chrom_remap.emplace(kv.first, kv.second);
        // VCF uses full name, FAI uses short ("CHM13#0#chr20" → "chr20")
        for (const auto& kv : fai_full_to_suffix) {
            if (faidx_has_seq(fai.get(), kv.second.c_str()))
                chrom_remap.emplace(kv.first, kv.second);
        }
    }

    // 5. contig name → FAI-order tid (used when remapping candidate tids for output).
    std::unordered_map<std::string, int> contig_to_tid;
    contig_to_tid.reserve(static_cast<size_t>(header->n_targets));
    for (int32_t tid = 0; tid < header->n_targets; ++tid) {
        contig_to_tid[header->target_name[tid]] = tid;
    }

    // 6. Tile genome into RegionChunks using the resolved region_filters (which may have
    //    contig names like "CHM13#0#chr20" resolved from a user-supplied "chr20").
    // Alignment recovery, when a BAM is given. One context for the whole run:
    // it owns the BAM handles the targeted sub-solves reopen per window.
    std::unique_ptr<WorkerContext> bam_recovery_ctx;
    if (!opts.bam_files.empty()) {
        bam_recovery_ctx = std::make_unique<WorkerContext>(opts);
        std::cerr << "Alignment recovery enabled from " << opts.primary_bam_file() << "\n";
    }
    const std::vector<RegionChunk> chunks =
        build_region_chunks(opts, header.get(), fai.get(), region_filters);
    if (chunks.empty()) {
        std::cerr << "No region chunks to process\n";
        return;
    }
    if (opts.verbose >= 1) {
        std::cerr << "Tiled genome into " << chunks.size() << " chunks ("
                  << opts.chunk_size << " bp, " << opts.threads << " thread(s))\n";
    }

    // 7. Build per-chunk query config for the GBZ/GAF-base FFI path.
    GraphQueryConfig qconfig;
    qconfig.gbz_db   = opts.gbz_db;
    qconfig.gaf_db   = opts.gaf_db;
    qconfig.min_mapq  = opts.min_mapq;

    // Derive the reference sample name for GBZ interval queries.
    // For pangenome FASTAs ("CHM13#0#chr20") extract "CHM13".
    // For plain FASTAs ("chr20") leave empty so query uses its GENERIC_SAMPLE default.
    std::string ref_sample = opts.graph_sample;
    if (ref_sample.empty()) {
        const int nseq = faidx_nseq(fai.get());
        for (int i = 0; i < nseq && ref_sample.empty(); ++i) {
            const std::string name(faidx_iseq(fai.get(), i));
            const size_t h = name.find('#');
            if (h != std::string::npos && h > 0)
                ref_sample = name.substr(0, h);
        }
    }
    if (!use_indexed_gaf && opts.verbose >= 1 && !ref_sample.empty())
        std::cerr << "Using reference sample \"" << ref_sample << "\" for GBZ interval queries\n";

    // 8. Open output streams.
    std::ofstream variant_out(opts.output_tsv);
    if (!variant_out) throw std::runtime_error("failed to open output: " + opts.output_tsv);
    write_variants_tsv_header(variant_out);

    std::ofstream vcf_out;
    if (!opts.output_vcf.empty()) {
        vcf_out.open(opts.output_vcf);
        if (!vcf_out) throw std::runtime_error("failed to open VCF output: " + opts.output_vcf);
        write_variants_vcf_header(vcf_out, opts, header.get());
    }
    std::ofstream phased_vcf_out;
    if (!opts.output_phased_vcf.empty()) {
        phased_vcf_out.open(opts.output_phased_vcf);
        if (!phased_vcf_out)
            throw std::runtime_error("failed to open phased VCF output: " + opts.output_phased_vcf);
        write_phased_variants_vcf_header(phased_vcf_out, opts, header.get());
    }

    ReferenceCache ref(fai.get());

    // Phased BAM: open output file and write header.
    struct SamFileCloser { void operator()(samFile* fp) const { if (fp) hts_close(fp); } };
    std::unique_ptr<samFile, SamFileCloser> phased_bam_fp;
    std::unique_ptr<sam_hdr_t, decltype(&sam_hdr_destroy)> phased_bam_hdr(nullptr, &sam_hdr_destroy);
    std::unordered_map<std::string, PhaseReadOutputRow> phased_bam_rows;
    std::unordered_set<std::string> phased_bam_emitted;
    const bool emit_phased_bam = !opts.output_phased_bam.empty();
    if (emit_phased_bam) {
        samFile* raw_fp = hts_open(opts.output_phased_bam.c_str(), "wb");
        if (!raw_fp)
            throw std::runtime_error("failed to open phased BAM: " + opts.output_phased_bam);
        phased_bam_fp.reset(raw_fp);
        sam_hdr_t* hdr = sam_hdr_init();
        if (!hdr) throw std::runtime_error("failed to allocate phased BAM header");
        phased_bam_hdr.reset(hdr);
        if (sam_hdr_write(phased_bam_fp.get(), phased_bam_hdr.get()) < 0)
            throw std::runtime_error("failed to write phased BAM header");
    }

    // 10. Process chunks in reg_chunk_i batches (one contig per batch), streaming output.
    //     Mirrors run_collect_bam_variation's batch loop exactly.
    size_t n_variants = 0;
    size_t n_filtered = 0;
    // Diagnostic: why catalog sites never became candidates. Streamed alongside
    // the batch loop so a whole-chromosome run does not buffer millions of rows.
    std::unique_ptr<std::FILE, int (*)(std::FILE*)> filtered_out(nullptr, std::fclose);
    if (!opts.output_filtered_sites.empty()) {
        std::FILE* raw = std::fopen(opts.output_filtered_sites.c_str(), "w");
        if (raw == nullptr)
            throw std::runtime_error("failed to open filtered sites file: " +
                                     opts.output_filtered_sites);
        filtered_out.reset(raw);
        std::fprintf(filtered_out.get(),
                     "CHROM\tPOS\tSITE_ID\tREF_COV\tALT_COV\tTOTAL_COV\tAF\tREASON\n");
    }
    std::ofstream phase_sites_out;
    if (!opts.output_phase_sites.empty()) {
        phase_sites_out.open(opts.output_phase_sites);
        if (!phase_sites_out)
            throw std::runtime_error("failed to open phase sites file: " +
                                     opts.output_phase_sites);
        write_graph_phase_sites_tsv_header(phase_sites_out);
    }
    // Per-read phasing evidence, accumulated across chunks. Used to diagnose
    // which reads get a haplotype on thin or contradictory evidence: a read is
    // assigned by init_assign_read_hap_based_on_cons_alle with no minimum-observation or margin
    // requirement, so a single informative site is enough to commit it.
    struct PhaseReadDiag {
        int hap = 0;
        hts_pos_t phase_set = kUnphasedReadPhaseSet;
        int n_obs = 0;
        int agree = 0;
        int conflict = 0;
        int score_margin = 0;
        int n_scored = 0;
    };
    std::unordered_map<std::string, PhaseReadDiag> phase_read_diag;
    const bool emit_phase_reads = !opts.output_phase_reads.empty();

    size_t batch_begin = 0;
    while (batch_begin < chunks.size()) {
        size_t batch_end = batch_begin + 1;
        while (batch_end < chunks.size() &&
               chunks[batch_end].reg_chunk_i == chunks[batch_begin].reg_chunk_i) {
            ++batch_end;
        }

        std::vector<GraphChunkBuildResult> graph_chunks =
            use_indexed_gaf
                ? process_graph_chunk_batch_indexed_gaf(
                      opts.graph_sites_vcf, chunks, batch_begin, batch_end,
                      header.get(), opts.gaf_file, opts.min_mapq,
                      fai_full_to_suffix, chrom_remap, opts,
                      pgbam_sidecar.get())
                : process_graph_chunk_batch(
                      opts.graph_sites_vcf, chunks, batch_begin, batch_end,
                      header.get(), qconfig, ref_sample, fai_full_to_suffix,
                      chrom_remap, opts, pgbam_sidecar.get());


        if (emit_phase_reads) {
            for (const GraphChunkBuildResult& gc : graph_chunks) {
                const PhasingChunk& pc = gc.chunk;
                for (const ReadVariantProfile& profile : pc.read_var_profile) {
                    const size_t read_i = static_cast<size_t>(profile.read_id);
                    if (read_i >= pc.reads.size()) continue;
                    const ReadRecord& rr = pc.reads[read_i];
                    PhaseReadDiag& d = phase_read_diag[rr.qname];
                    int obs = 0;
                    for (int allele : profile.alleles)
                        if (allele >= 0) ++obs;
                    d.n_obs += obs;
                    d.agree += rr.n_clean_agree_snps;
                    d.conflict += rr.n_clean_conflict_snps;
                    if (rr.hap_score_margin > d.score_margin)
                        d.score_margin = rr.hap_score_margin;
                    if (rr.n_vars_scored > d.n_scored) d.n_scored = rr.n_vars_scored;
                    const int hap = read_i < pc.haps.size() ? pc.haps[read_i] : 0;
                    if (hap != 0) {
                        d.hap = hap;
                        d.phase_set = read_i < pc.phase_sets.size()
                                          ? pc.phase_sets[read_i]
                                          : kUnphasedReadPhaseSet;
                    }
                }
            }
        }

        for (const GraphChunkBuildResult& gc : graph_chunks) {
            n_filtered += gc.filtered_sites.size();
            if (filtered_out) {
                for (const FilteredGraphSite& fs : gc.filtered_sites) {
                    std::fprintf(filtered_out.get(), "%s\t%lld\t%s\t%d\t%d\t%d\t%.4f\t%s\n",
                                 fs.chrom.c_str(), static_cast<long long>(fs.pos),
                                 fs.site_id.c_str(), fs.ref_cov, fs.alt_cov,
                                 fs.total_cov, fs.allele_fraction,
                                 fs.filter_reason.c_str());
                }
            }
            if (phase_sites_out) {
                write_graph_phase_sites_tsv_rows(phase_sites_out, gc);
            }
        }

        CandidateTable variants =
            graph_chunks_to_candidate_table(graph_chunks, contig_to_tid, opts);

        n_variants += variants.size();
        write_variants_tsv_records(variant_out, header.get(), ref, variants);
        if (!opts.output_vcf.empty())
            write_variants_vcf_records(vcf_out, opts, header.get(), ref, variants);
        if (!opts.output_phased_vcf.empty())
            write_phased_variants_vcf_records(phased_vcf_out, opts, header.get(), ref, variants);

        // Phased BAM: accumulate post-stitch read assignments, then flush.
        // Batches are per-contig so cross-batch read overlap is impossible;
        // flush everything after each batch.
        if (emit_phased_bam) {
            for (const GraphChunkBuildResult& gc : graph_chunks)
                merge_graph_chunk_into_read_rows(phased_bam_rows, gc,
                                                 opts.min_read_hap_margin);
            flush_graph_phase_bam_after_merge(
                phased_bam_fp.get(), phased_bam_hdr.get(),
                phased_bam_rows, nullptr, phased_bam_emitted);
        }

        batch_begin = batch_end;
    }

    if (emit_phase_reads) {
        std::FILE* fp = std::fopen(opts.output_phase_reads.c_str(), "w");
        if (fp == nullptr)
            throw std::runtime_error("failed to open phase reads file: " +
                                     opts.output_phase_reads);
        std::fprintf(fp, "READ\tHAP\tPHASE_SET\tN_OBS\tCLEAN_AGREE\tCLEAN_CONFLICT\tSCORE_MARGIN\tN_SCORED\n");
        for (const auto& [qname, d] : phase_read_diag) {
            std::fprintf(fp, "%s\t%d\t%lld\t%d\t%d\t%d\t%d\t%d\n", qname.c_str(), d.hap,
                         static_cast<long long>(d.phase_set), d.n_obs, d.agree, d.conflict,
                         d.score_margin, d.n_scored);
        }
        std::fclose(fp);
        std::cerr << "Wrote per-read phasing evidence to " << opts.output_phase_reads << "\n";
    }

    std::cerr << "Processed " << chunks.size() << " region chunks with " << opts.threads
              << " worker thread(s)\n";
    std::cerr << "Collected " << n_variants << " candidate variant sites ("
              << n_filtered << " filtered) into " << opts.output_tsv << "\n";
    if (!opts.output_vcf.empty())
        std::cerr << "Wrote candidate VCF to " << opts.output_vcf << "\n";
    if (!opts.output_phased_vcf.empty())
        std::cerr << "Wrote phased candidate VCF to " << opts.output_phased_vcf << "\n";
    if (!opts.output_phase_sites.empty())
        std::cerr << "Wrote retained graph sites to " << opts.output_phase_sites << "\n";
    if (emit_phased_bam)
        std::cerr << "Wrote phased BAM to " << opts.output_phased_bam << "\n";
}

static void print_graph_collect_help() {
    std::cout
        << "Usage: pgphase collect-graph-variation [options]\n"
        << "\n"
        << "Required:\n"
        << "      --ref FILE                Reference FASTA (indexed)\n"
        << "      --sites FILE              Sites VCF from build-snarl-catalog (bgzipped + tabix-indexed)\n"
        << "\n"
        << "Read input:\n"
        << "      --gaf FILE                Coordinate-indexed GAF from pggaf (bgzipped + tabix-indexed)\n"
        << "      --gbz-db FILE             GBZ graph database (legacy --gaf-db path)\n"
        << "      --gaf-db FILE             GAF-base read alignment database (legacy path)\n"
        << "\n"
        << "Options:\n"
        << "  -o, --output FILE             Output TSV [output.tsv]\n"
        << "  -v, --vcf-output FILE         Candidate VCF output\n"
        << "      --phased-vcf-out FILE     Phased VCF with GT:DP:AD:VAF:GQ:PS\n"
        << "      --phased-bam-out FILE     Unaligned BAM with HP/PS tags per read\n"
        << "      --recovery-audit-out FILE One row per candidate the recovery sub-solve\n"
        << "                                found, and what the merge did with it\n"
        << "      --link-earned-repeat-indels  Re-admit a repeat-context het indel when it\n"
        << "                                agrees with a nearby clean het SNP on >= 15 reads\n"
        << "      --bam FILE                Indexed BAM used to recover seams between\n"
        << "                                neighboring graph phase sets; imports sites\n"
        << "                                as independent local phase blocks, then stitches\n"
        << "                                left to right on decisive allele evidence\n"
        << "      --filtered-sites-out FILE Diagnostic TSV of dropped catalog sites and why\n"
        << "      --phase-sites-out FILE    Diagnostic TSV of retained graph sites with SITE_ID\n"
        << "      --phase-reads-out FILE    Diagnostic TSV of per-read phasing evidence\n"
        << "      --phase-matrix-dump PATH  Dump phasing inputs and incoming assignments\n"
        << "      --graph-indel-af-margin F  Max |AF-0.5| for a het-indel k-means anchor [0.11]\n"
        << "      --graph-indel-min-alt INT  Min alt support for a het-indel k-means anchor [0]\n"
        << "      --min-read-margin INT     Min clean-SNP (agree-conflict) to phase a read [0=off]\n"
        << "      --stitch-min-margin INT   Abstain on chunk seams below this vote margin [0]\n"
        << "      --stitch-rule INT         0=net-margin 1=both-strands 2=literal 3=both+margin [0]\n"
        << "      --anchor-af-margin F      Max |AF-0.5| for a site to vote in k-means [0.5=off]\n"
        << "      --min-block-link-reads INT Spanning reads needed to carry a phase block [2]\n"
        << "      --block-link-window INT   Preceding het variants searched for that link [1]\n"
        << "      --link-by-alleles         Let untagged reads link variants by allele pattern\n"
        << "      --emit-nonanchor-hets     Emit/phase hets outside --anchor-af-margin (never anchor)\n"
        << "      --gaf-pad INT             Widen the per-chunk GAF read query by INT bp [0]\n"
        << "      --snarl-allele-phasing    Score multi-allelic snarls as alt-vs-other, not alt-vs-ref\n"
        << "      --snarl-keep-whole        Keep multi-allelic snarls as single n-allelic anchors\n"
        << "      --snarl-top2-frac FLOAT   Min read share on a snarl's top 2 alleles to anchor [0.9]\n"
        << "      --af-vs-site-depth        Score allele fraction against total site depth,\n"
        << "                                recovering hets between two non-reference alleles\n"
        << "  -t, --threads INT             Worker threads [1]\n"
        << "  -q, --min-mapq INT            Minimum read mapping quality [5]\n"
        << "  -D, --min-depth INT           Minimum total depth [5]\n"
        << "      --min-alt-depth INT       Minimum alt depth [2]\n"
        << "      --min-af FLOAT            Minimum allele fraction [0.20]\n"
        << "      --max-af FLOAT            Maximum allele fraction [0.80]\n"
        << "      --min-sv-len INT          Min SV length for SVTYPE/SVLEN tags [30]\n"
        << "      --chunk-size INT          Region chunk size [500000; 1000000 with --bam]\n"
        << "  -r, --region STR              Restrict to region (may be repeated)\n"
        << "      --region-file FILE        BED file of regions\n"
        << "      --autosome                Process chr1-22 / 1-22 only\n"
        << "      --sample NAME             Reference sample name for GBZ interval queries\n"
        << "                                (auto-derived from FASTA if not provided)\n"

        << "      --hifi                    HiFi read mode [default]\n"
        << "      --ont                     ONT read mode (enables strand-bias filter)\n"
        << "      --strand-bias-pval FLOAT  Max p-value for ONT strand-bias filter [0.01]\n"
        << "\n"
        << "Pgbam stitching:\n"
        << "      --pgbam-file FILE         Optional .pgbam sidecar for fallback chunk stitching\n"
        << "      --pgbam-primary-margin INT         Thread polarity margin for primary stitching [2]\n"
        << "      --pgbam-primary-min-winning INT    Winning shared polarized threads for primary stitching [2]\n"
        << "      --no-pgbam-cleanup-pass            Disable final .pgbam cleanup pass\n"
        << "      --pgbam-cleanup-margin INT         Thread polarity margin for cleanup pass [2]\n"
        << "      --pgbam-cleanup-min-winning INT    Winning shared polarized threads for cleanup pass [1]\n"
        << "      --no-pgbam-relaxed-cleanup-pass    Disable relaxed .pgbam cleanup pass\n"
        << "      --pgbam-relaxed-cleanup-margin INT Thread polarity margin for relaxed cleanup [1]\n"
        << "      --pgbam-relaxed-cleanup-min-winning INT Winning threads for relaxed cleanup [1]\n"
        << "\n"
        << "  -V, --verbose INT             Verbosity level [0]\n"
        << "  -h, --help                    Print this help\n"
        << "\n"
        << "Examples:\n"
        << "  pgphase collect-graph-variation \\\n"
        << "      --ref ref.fa \\\n"
        << "      --sites sites.vcf.gz \\\n"
        << "      --gaf reads.gaf \\\n"
        << "      --phased-vcf-out phased.vcf \\\n"
        << "      --phased-bam-out phased.bam \\\n"
        << "      -t 8\n"
        << "\n"
        << "  pgphase collect-graph-variation \\\n"
        << "      --ref ref.fa \\\n"
        << "      --sites sites.vcf.gz \\\n"
        << "      --gbz-db reads.gaf.db \\\n"
        << "      --ont \\\n"
        << "      --phased-vcf-out phased.vcf \\\n"
        << "      -t 16\n";
}

enum GraphCollectOption {
    kGcRecoveryBam = 2000,
    kGcMinAltDepth = 1000,
    kGcMinAf,
    kGcMaxAf,
    kGcMinSvLen,
    kGcChunkSize,
    kGcPhasedVcf,
    kGcGbzDb,
    kGcGafFile,
    kGcGafDb,
    kGcRegionFile,
    kGcAutosome,
    kGcSample,
    kGcOnt,
    kGcHifi,
    kGcStrandBiasPval,
    kGcPhasedBam,
    kGcLinkEarnedRepeatIndels,
    kGcRecoveryAuditOut,
    kGcRef,
    kGcSites,
    kGcPgbamFile,
    kGcPgbamPrimaryMargin,
    kGcPgbamPrimaryMinWinning,
    kGcNoPgbamCleanupPass,
    kGcPgbamCleanupMargin,
    kGcPgbamCleanupMinWinning,
    kGcNoPgbamRelaxedCleanupPass,
    kGcPgbamRelaxedCleanupMargin,
    kGcPgbamRelaxedCleanupMinWinning,
    kGcFilteredSitesOut,
    kGcPhaseSitesOut,
    kGcPhaseReadsOut,
    kGcGraphIndelAfMargin,
    kGcGraphIndelMinAlt,
    kGcMinReadHapMargin,
    kGcStitchMinMargin,
    kGcStitchRule,
    kGcAnchorAfMargin,
    kGcAfVsSiteDepth,
    kGcBlockLink,
    kGcBlockLinkWindow,
    kGcLinkByAlleles,
    kGcEmitNonAnchorHets,
    kGcGafPad,
    kGcSnarlAllelePhasing,
    kGcSnarlKeepWhole,
    kGcSnarlTop2Frac,
    kGcPhaseMatrixDump,
};

} // namespace

} // namespace pgphase_collect

int collect_graph_variation(int argc, char* argv[]) {
    using namespace pgphase_collect;
    Options opts;
    opts.min_mapq = kDefaultGraphMinMapq;
    bool chunk_size_explicit = false;

    {
        std::ostringstream cmd;
        cmd << "pgphase collect-graph-variation";
        for (int i = 1; i < argc; ++i) cmd << ' ' << argv[i];
        opts.command_line = cmd.str();
    }

    optind = 1;
    const struct option long_options[] = {
        {"output",            required_argument, nullptr, 'o'},
        {"vcf-output",        required_argument, nullptr, 'v'},
        {"phased-vcf-out",    required_argument, nullptr, kGcPhasedVcf},
        {"phased-bam-out",   required_argument, nullptr, kGcPhasedBam},
        {"link-earned-repeat-indels", no_argument, nullptr, kGcLinkEarnedRepeatIndels},
        {"recovery-audit-out", required_argument, nullptr, kGcRecoveryAuditOut},
        {"bam",              required_argument, nullptr, kGcRecoveryBam},
        {"filtered-sites-out", required_argument, nullptr, kGcFilteredSitesOut},
        {"phase-sites-out",   required_argument, nullptr, kGcPhaseSitesOut},
        {"phase-reads-out",   required_argument, nullptr, kGcPhaseReadsOut},
        {"phase-matrix-dump", required_argument, nullptr, kGcPhaseMatrixDump},
        {"graph-indel-af-margin", required_argument, nullptr, kGcGraphIndelAfMargin},
        {"graph-indel-min-alt",   required_argument, nullptr, kGcGraphIndelMinAlt},
        {"min-read-margin",   required_argument, nullptr, kGcMinReadHapMargin},
        {"stitch-min-margin", required_argument, nullptr, kGcStitchMinMargin},
        {"stitch-rule",       required_argument, nullptr, kGcStitchRule},
        {"anchor-af-margin",  required_argument, nullptr, kGcAnchorAfMargin},
        {"af-vs-site-depth",  no_argument,       nullptr, kGcAfVsSiteDepth},
        {"min-block-link-reads", required_argument, nullptr, kGcBlockLink},
        {"block-link-window",    required_argument, nullptr, kGcBlockLinkWindow},
        {"link-by-alleles",      no_argument,       nullptr, kGcLinkByAlleles},
        {"emit-nonanchor-hets",  no_argument,       nullptr, kGcEmitNonAnchorHets},
        {"gaf-pad",              required_argument, nullptr, kGcGafPad},
        {"snarl-allele-phasing", no_argument,       nullptr, kGcSnarlAllelePhasing},
        {"snarl-keep-whole",     no_argument,       nullptr, kGcSnarlKeepWhole},
        {"snarl-top2-frac",      required_argument, nullptr, kGcSnarlTop2Frac},
        {"threads",           required_argument, nullptr, 't'},
        {"min-mapq",          required_argument, nullptr, 'q'},
        {"min-depth",         required_argument, nullptr, 'D'},
        {"min-alt-depth",     required_argument, nullptr, kGcMinAltDepth},
        {"min-af",            required_argument, nullptr, kGcMinAf},
        {"max-af",            required_argument, nullptr, kGcMaxAf},
        {"min-sv-len",        required_argument, nullptr, kGcMinSvLen},
        {"chunk-size",        required_argument, nullptr, kGcChunkSize},
        {"region",            required_argument, nullptr, 'r'},
        {"region-file",       required_argument, nullptr, kGcRegionFile},
        {"autosome",          no_argument,       nullptr, kGcAutosome},
        {"gaf",               required_argument, nullptr, kGcGafFile},
        {"gaf-file",          required_argument, nullptr, kGcGafFile},
        {"gbz-db",            required_argument, nullptr, kGcGbzDb},
        {"gaf-db",            required_argument, nullptr, kGcGafDb},
        {"sample",            required_argument, nullptr, kGcSample},
        {"hifi",              no_argument,       nullptr, kGcHifi},
        {"ont",               no_argument,       nullptr, kGcOnt},
        {"strand-bias-pval",  required_argument, nullptr, kGcStrandBiasPval},
        {"ref",               required_argument, nullptr, kGcRef},
        {"sites",             required_argument, nullptr, kGcSites},
        {"pgbam-file",                required_argument, nullptr, kGcPgbamFile},
        {"pgbam-primary-margin",      required_argument, nullptr, kGcPgbamPrimaryMargin},
        {"pgbam-primary-min-winning", required_argument, nullptr, kGcPgbamPrimaryMinWinning},
        {"no-pgbam-cleanup-pass",     no_argument,       nullptr, kGcNoPgbamCleanupPass},
        {"pgbam-cleanup-margin",      required_argument, nullptr, kGcPgbamCleanupMargin},
        {"pgbam-cleanup-min-winning", required_argument, nullptr, kGcPgbamCleanupMinWinning},
        {"no-pgbam-relaxed-cleanup-pass", no_argument,   nullptr, kGcNoPgbamRelaxedCleanupPass},
        {"pgbam-relaxed-cleanup-margin", required_argument, nullptr, kGcPgbamRelaxedCleanupMargin},
        {"pgbam-relaxed-cleanup-min-winning", required_argument, nullptr, kGcPgbamRelaxedCleanupMinWinning},
        {"verbose",           required_argument, nullptr, 'V'},
        {"help",              no_argument,       nullptr, 'h'},
        {nullptr, 0, nullptr, 0}
    };

    int opt = 0;
    int long_index = 0;
    while ((opt = getopt_long(argc, argv, "o:v:t:q:D:r:V:h", long_options, &long_index)) != -1) {
        switch (opt) {
            case 'o': opts.output_tsv = optarg; break;
            case 'v': opts.output_vcf = optarg; break;
            case kGcPhasedVcf:    opts.output_phased_vcf = optarg; break;
            case kGcPhasedBam:    opts.output_phased_bam = optarg; break;
            case kGcLinkEarnedRepeatIndels: opts.link_earned_repeat_indels = true; break;
            case kGcRecoveryAuditOut: opts.recovery_audit_out = optarg; break;
            // Recovery only. Targeted BAM solves supply private sites and
            // local phase blocks for seams the graph catalog could not join.
            case kGcRecoveryBam:  opts.bam_files.push_back(optarg); break;

            case kGcFilteredSitesOut: opts.output_filtered_sites = optarg; break;
            case kGcPhaseSitesOut: opts.output_phase_sites = optarg; break;
            case kGcPhaseReadsOut: opts.output_phase_reads = optarg; break;
            case kGcPhaseMatrixDump: opts.phase_matrix_dump_prefix = optarg; break;
            case kGcGraphIndelAfMargin:
                opts.graph_indel_af_margin = parse_double_arg(optarg, "--graph-indel-af-margin");
                break;
            case kGcGraphIndelMinAlt:
                opts.graph_indel_min_alt = parse_int_arg(optarg, "--graph-indel-min-alt");
                break;
            case kGcMinReadHapMargin:
                opts.min_read_hap_margin = parse_int_arg(optarg, "--min-read-margin");
                break;
            case kGcStitchMinMargin:
                opts.stitch_min_margin = parse_int_arg(optarg, "--stitch-min-margin");
                break;
            case kGcStitchRule:
                opts.stitch_rule = parse_int_arg(optarg, "--stitch-rule");
                break;
            case kGcAnchorAfMargin:
                opts.anchor_af_margin = parse_double_arg(optarg, "--anchor-af-margin");
                break;
            case kGcAfVsSiteDepth: opts.af_vs_site_depth = true; break;
            case kGcBlockLink:
                opts.min_block_link_reads = parse_int_arg(optarg, "--min-block-link-reads");
                break;
            case kGcBlockLinkWindow:
                opts.block_link_window = parse_int_arg(optarg, "--block-link-window");
                break;
            case kGcLinkByAlleles: opts.link_by_alleles = true; break;
            case kGcEmitNonAnchorHets: opts.emit_nonanchor_hets = true; break;
            case kGcGafPad: opts.gaf_pad = parse_int_arg(optarg, "--gaf-pad"); break;
            case kGcSnarlAllelePhasing: opts.snarl_allele_phasing = true; break;
            case kGcSnarlKeepWhole: opts.snarl_keep_whole = true; opts.snarl_allele_phasing = true; break;
            case kGcSnarlTop2Frac: opts.snarl_top2_frac = parse_double_arg(optarg, "--snarl-top2-frac"); break;
            case 't': opts.threads = parse_int_arg(optarg, "--threads"); break;
            case 'q': opts.min_mapq = parse_int_arg(optarg, "--min-mapq"); break;
            case 'D': opts.min_depth = parse_int_arg(optarg, "--min-depth"); break;
            case kGcMinAltDepth:  opts.min_alt_depth = parse_int_arg(optarg, "--min-alt-depth"); break;
            case kGcMinAf:        opts.min_af = parse_double_arg(optarg, "--min-af"); break;
            case kGcMaxAf:        opts.max_af = parse_double_arg(optarg, "--max-af"); break;
            case kGcMinSvLen:     opts.min_sv_len = parse_int_arg(optarg, "--min-sv-len"); break;
            case kGcChunkSize:
                opts.chunk_size = parse_ll_arg(optarg, "--chunk-size");
                chunk_size_explicit = true;
                break;
            case 'r': opts.regions.push_back(optarg); break;
            case kGcRegionFile:   opts.region_file = optarg; break;
            case kGcAutosome:     opts.autosome = true; break;
            case kGcGafFile:      opts.gaf_file = optarg; break;
            case kGcGbzDb:        opts.gbz_db = optarg; break;
            case kGcGafDb:        opts.gaf_db = optarg; break;
            case kGcSample:       opts.graph_sample = optarg; break;
            case kGcHifi:         opts.read_technology = ReadTechnology::Hifi; break;
            case kGcOnt:          opts.read_technology = ReadTechnology::Ont; break;
            case kGcStrandBiasPval: opts.strand_bias_pval = parse_double_arg(optarg, "--strand-bias-pval"); break;
            case kGcRef:          opts.ref_fasta = optarg; break;
            case kGcSites:        opts.graph_sites_vcf = optarg; break;
            case kGcPgbamFile:      opts.pgbam_file = optarg; break;
            case kGcPgbamPrimaryMargin: opts.pgbam_primary_polarity_margin = parse_int_arg(optarg, "--pgbam-primary-margin"); break;
            case kGcPgbamPrimaryMinWinning: opts.pgbam_primary_min_winning_threads = parse_int_arg(optarg, "--pgbam-primary-min-winning"); break;
            case kGcNoPgbamCleanupPass: opts.pgbam_cleanup_pass = false; break;
            case kGcPgbamCleanupMargin: opts.pgbam_cleanup_polarity_margin = parse_int_arg(optarg, "--pgbam-cleanup-margin"); break;
            case kGcPgbamCleanupMinWinning: opts.pgbam_cleanup_min_winning_threads = parse_int_arg(optarg, "--pgbam-cleanup-min-winning"); break;
            case kGcNoPgbamRelaxedCleanupPass: opts.pgbam_relaxed_cleanup_pass = false; break;
            case kGcPgbamRelaxedCleanupMargin: opts.pgbam_relaxed_cleanup_polarity_margin = parse_int_arg(optarg, "--pgbam-relaxed-cleanup-margin"); break;
            case kGcPgbamRelaxedCleanupMinWinning: opts.pgbam_relaxed_cleanup_min_winning_threads = parse_int_arg(optarg, "--pgbam-relaxed-cleanup-min-winning"); break;
            case 'V': opts.verbose = parse_int_arg(optarg, "--verbose"); break;
            case 'h': print_graph_collect_help(); return 0;
            default:  print_graph_collect_help(); return 1;
        }
    }

    if (!chunk_size_explicit && !opts.primary_bam_file().empty())
        opts.chunk_size = kDefaultGraphRecoveryChunkSize;

    if (opts.ref_fasta.empty() || opts.graph_sites_vcf.empty()) {
        std::cerr << "Error: --ref and --sites are required\n";
        print_graph_collect_help();
        return 1;
    }

    require_graph_site_vcf_tabix_index(opts.graph_sites_vcf);

    if (opts.gaf_file.empty() && (opts.gbz_db.empty() || opts.gaf_db.empty())) {
        std::cerr << "Error: provide --gaf, or provide --gbz-db and --gaf-db\n";
        print_graph_collect_help();
        return 1;
    }

    if (opts.threads < 1 || opts.min_mapq < 0 || opts.min_depth < 0 ||
        opts.min_alt_depth < 0 || opts.min_af < 0.0 || opts.max_af < opts.min_af ||
        opts.min_sv_len < 0 || opts.chunk_size < 1 || opts.verbose < 0 ||
        opts.strand_bias_pval < 0.0 || opts.strand_bias_pval > 1.0) {
        std::cerr << "Error: invalid numeric threshold\n";
        return 1;
    }

    try {
        run_collect_graph_variation(opts);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
