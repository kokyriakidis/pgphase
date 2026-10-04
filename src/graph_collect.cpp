#include "graph_collect.hpp"
#include "collect_pipeline.hpp"

#include "arg_parse.hpp"
#include "collect_output.hpp"
#include "collect_phase.hpp"
#include "collect_phase_pgbam.hpp"
#include "collect_phase_noisy.hpp"
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
#include <cctype>
#include <atomic>
#include <cstdio>
#include <cstdint>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <memory>
#include <limits>
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
static constexpr double kDefaultSourceGraphMaxP = 0.01;
static bool source_graph_vote_supported(
        int same, int cross, int parity,
        double max_p = kDefaultSourceGraphMaxP) {
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
    return tail <= max_p;
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

// Two complete BAM source paths can meet inside a graph seam even when the
// outer graph-to-graph transaction abstains. Join their current blocks only
// through a verified deletion-to-SNP edge and a direct SNP-to-graph anchor
// link. The first SNP in the downstream source is chosen by coordinate before
// examining alleles; later noisy sites cannot be searched for a favorable vote.
static void stitch_complete_bam_sources_at_repeat_cut(GraphChunkBuildResult& gc) {
    PhasingChunk& chunk = gc.chunk;
    constexpr int kMinBridgeMapq = 30;
    constexpr double kMaxBridgeP = 0.01;
    for (const RecoverySeam& seam : gc.recovery_windows) {
        std::vector<const RecoverySourceSite*> sites;
        for (const RecoverySourceSite& site : gc.recovery_source_sites) {
            if (!site.can_adopt || site.candidate_index >= chunk.candidates.size())
                continue;
            const CandidateVariant& candidate =
                chunk.candidates[site.candidate_index];
            const hts_pos_t pos = candidate.key.sort_pos();
            if (pos < seam.beg || pos > seam.end ||
                candidate.phase_set <= 0 || !candidate.bam_injected)
                continue;
            const auto supported =
                gc.recovery_source_path_supported.find(site.phase_set);
            if (supported != gc.recovery_source_path_supported.end() &&
                supported->second)
                sites.push_back(&site);
        }
        std::sort(sites.begin(), sites.end(),
                  [](const RecoverySourceSite* a,
                     const RecoverySourceSite* b) {
                      return a->candidate_index < b->candidate_index;
                  });
        for (size_t i = 0; i + 1 < sites.size(); ++i) {
            const RecoverySourceSite& left_site = *sites[i];
            const RecoverySourceSite& next_site = *sites[i + 1];
            if (left_site.phase_set == next_site.phase_set) continue;
            const size_t left_i = left_site.candidate_index;
            const CandidateVariant& left = chunk.candidates[left_i];
            if (left.key.type != VariantType::Deletion ||
                !left.is_homopolymer_indel || !left.msa_verified ||
                !left.alignment_verified)
                continue;
            std::optional<size_t> right_i;
            for (size_t j = i + 1; j < sites.size() &&
                 sites[j]->phase_set == next_site.phase_set; ++j) {
                const size_t candidate_i = sites[j]->candidate_index;
                const CandidateVariant& candidate =
                    chunk.candidates[candidate_i];
                if (candidate.key.type == VariantType::Snp &&
                    candidate.msa_verified && candidate.alignment_verified) {
                    right_i = candidate_i;
                    break;
                }
            }
            if (!right_i) continue;
            const CandidateVariant& right = chunk.candidates[*right_i];
            const hts_pos_t left_ps = left.phase_set;
            const hts_pos_t right_ps = right.phase_set;
            if (left_ps == right_ps ||
                left.key.sort_pos() >= right.key.sort_pos())
                continue;
            std::optional<size_t> graph_anchor;
            for (size_t ci = *right_i + 1; ci < chunk.candidates.size(); ++ci) {
                const CandidateVariant& candidate = chunk.candidates[ci];
                if (candidate.phase_set == right_ps &&
                    !candidate.bam_injected &&
                    candidate.counts.category == VariantCategory::CleanHetSnp &&
                    candidate.hap_to_cons_alle[1] >= 0 &&
                    candidate.hap_to_cons_alle[2] >= 0 &&
                    candidate.hap_to_cons_alle[1] !=
                        candidate.hap_to_cons_alle[2]) {
                    graph_anchor = ci;
                    break;
                }
            }
            if (!graph_anchor) continue;
            const std::optional<bool> flip = local_run_boundary_flip(
                chunk, left_i, *right_i, kMinBridgeMapq, kMaxBridgeP);
            const std::optional<bool> anchored = local_run_boundary_flip(
                chunk, *right_i, *graph_anchor, kMinBridgeMapq, kMaxBridgeP);
            if (!flip || !anchored || *anchored) continue;
            if (merge_phase_sets_in_place(chunk, left_ps, right_ps, *flip))
                break;
        }
    }
}

// Graph catalog substitutions can span multiple bases. Their first changed
// base is the reference coordinate used by recovery and path ordering.
static std::optional<hts_pos_t> graph_substitution_start(
        const GraphSiteMeta& meta, const std::string& alt) {
    if (meta.ref.empty() || meta.ref.size() != alt.size())
        return std::nullopt;
    size_t first = 0;
    while (first < meta.ref.size() &&
           std::toupper(static_cast<unsigned char>(meta.ref[first])) ==
           std::toupper(static_cast<unsigned char>(alt[first])))
        ++first;
    if (first == meta.ref.size()) return std::nullopt;
    return meta.pos + static_cast<hts_pos_t>(first);
}

// Check the graph SNP path within each original block. The main BAM stitch
// may supply a certified edge between blocks when GAF has no reversal vote.
// A deletion bridge may also bypass one weak SNP using a significant direct
// edge between its neighbors; other internal edges need both haplotypes.
static bool graph_snp_path_supported(
        const GraphChunkBuildResult& gc, hts_pos_t phase_set,
        std::pair<hts_pos_t, hts_pos_t>* disconnected_pair = nullptr,
        const std::map<std::string, hts_pos_t>* original_graph_phase_sets = nullptr,
        const std::map<std::string, hts_pos_t>* stitched_graph_phase_sets = nullptr,
        bool allow_single_weak_site = false,
        bool allow_singleton = false,
        hts_pos_t min_pos = 0,
        bool report_one_hap_cut = false,
        bool allow_complete_bam_source = false,
        bool allow_source_first_edge = false) {
    if (disconnected_pair != nullptr) *disconnected_pair = {0, 0};
    constexpr int kMinPathMapq = 30;
    constexpr int kUnknownMapq = 255;
    constexpr int kMinSourceGraphAnchorVotes = 2;
    const PhasingChunk& chunk = gc.chunk;
    std::vector<std::pair<hts_pos_t, size_t>> sites;
    std::map<std::pair<hts_pos_t, std::string>, size_t> strongest_per_snarl;
    for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
        const CandidateVariant& candidate = chunk.candidates[ci];
        if (candidate.phase_set != phase_set || candidate.bam_injected ||
            candidate.counts.category != VariantCategory::CleanHetSnp ||
            candidate.hap_to_cons_alle[1] < 0 ||
            candidate.hap_to_cons_alle[2] < 0 ||
            candidate.hap_to_cons_alle[1] == candidate.hap_to_cons_alle[2] ||
            ci >= gc.site_meta.size())
            continue;
        const std::string* alt = selected_graph_candidate_alt(gc, ci);
        if (alt == nullptr) continue;
        const GraphSiteMeta& meta = gc.site_meta[ci];
        const VariantKey physical = vcf_to_variant_key(
            candidate.key.tid, meta.pos, meta.ref, *alt);
        if (physical.type != VariantType::Snp || physical.ref_len != 1 ||
            physical.alt.size() != 1 || physical.pos < min_pos ||
            ci >= gc.site_ids.size())
            continue;
        // Biallelic rows from one multiallelic snarl share a graph locus.
        // Requiring an independent two-haplotype edge between them invents a
        // boundary within that locus and can reject an otherwise sound block.
        std::string snarl_id = gc.site_ids[ci];
        const size_t colon = snarl_id.rfind(':');
        if (colon != std::string::npos && colon + 1 < snarl_id.size() &&
            snarl_id.find_first_not_of("0123456789", colon + 1) ==
                std::string::npos)
            snarl_id.resize(colon);
        const auto locus = std::make_pair(meta.pos, std::move(snarl_id));
        const auto [it, inserted] = strongest_per_snarl.emplace(locus, ci);
        if (!inserted && candidate.counts.total_cov >
                chunk.candidates[it->second].counts.total_cov)
            it->second = ci;
    }
    for (const auto& [locus, ci] : strongest_per_snarl) {
        (void)locus;
        const std::string* alt = selected_graph_candidate_alt(gc, ci);
        const GraphSiteMeta& meta = gc.site_meta[ci];
        const VariantKey physical = vcf_to_variant_key(
            chunk.candidates[ci].key.tid, meta.pos, meta.ref, *alt);
        sites.emplace_back(physical.pos, ci);
    }
    if (sites.size() < 2) {
        if (!allow_singleton || sites.size() != 1) return false;
        // A BAM source can extend one original graph SNP. There is no graph
        // edge within that singleton, but another graph row would need its
        // own path certificate before the whole phase set can be joined.
        return std::count_if(chunk.candidates.begin(), chunk.candidates.end(),
                             [phase_set](const CandidateVariant& candidate) {
                                 return candidate.phase_set == phase_set &&
                                     !candidate.bam_injected;
                             }) == 1;
    }
    std::sort(sites.begin(), sites.end());
    const auto pair_votes = [&](size_t left_i, size_t right_i) {
        std::array<int, 3> votes{};
        for (size_t ri = 0; ri < chunk.read_var_profile.size(); ++ri) {
            if (ri >= chunk.reads.size() || chunk.reads[ri].is_skipped ||
                chunk.reads[ri].mapq < kMinPathMapq ||
                chunk.reads[ri].mapq == kUnknownMapq)
                continue;
            const ReadVariantProfile& profile = chunk.read_var_profile[ri];
            if (profile.start_var_idx < 0 ||
                std::min(left_i, right_i) <
                    static_cast<size_t>(profile.start_var_idx) ||
                std::max(left_i, right_i) >
                    static_cast<size_t>(profile.end_var_idx))
                continue;
            const size_t left_offset =
                left_i - static_cast<size_t>(profile.start_var_idx);
            const size_t right_offset =
                right_i - static_cast<size_t>(profile.start_var_idx);
            if (std::max(left_offset, right_offset) >=
                profile.graph_alleles.size()) continue;
            const int left_allele = profile.graph_alleles[left_offset];
            const int right_allele = profile.graph_alleles[right_offset];
            if ((left_allele != 0 && left_allele != 1) ||
                (right_allele != 0 && right_allele != 1))
                continue;
            const int left_hap =
                left_allele == chunk.candidates[left_i].hap_to_cons_alle[1]
                    ? 1 : 2;
            const int right_hap =
                right_allele == chunk.candidates[right_i].hap_to_cons_alle[1]
                    ? 1 : 2;
            ++votes[left_hap == right_hap ? left_hap - 1 : 2];
        }
        return votes;
    };
    bool skipped_weak_site = false;
    bool supported_edge = false;
    for (size_t si = 1; si < sites.size(); ++si) {
        const auto [left_pos, left_i] = sites[si - 1];
        const auto [right_pos, right_i] = sites[si];
        if (left_pos >= right_pos) return false;
        const std::array<int, 3> votes = pair_votes(left_i, right_i);
        const int agree = votes[0] + votes[1];
        // At the start of a mixed block, an exact shared SNP pair can
        // anchor a complete BAM source even before the first GAF edge.
        // Require both rows to have the same source orientation and no cuts.
        if (allow_source_first_edge && agree == 0 && votes[2] == 0) {
            const RecoverySourceSite* left_source = nullptr;
            const RecoverySourceSite* right_source = nullptr;
            for (const RecoverySourceSite& site : gc.recovery_source_sites) {
                if (site.candidate_index == left_i) left_source = &site;
                if (site.candidate_index == right_i) right_source = &site;
            }
            if (left_source && right_source &&
                left_source->phase_set > 0 &&
                left_source->phase_set == right_source->phase_set &&
                left_source->clean_shared_snp &&
                right_source->clean_shared_snp &&
                left_source->can_adopt && right_source->can_adopt &&
                (left_source->hap1_allele ==
                    chunk.candidates[left_i].hap_to_cons_alle[1]) ==
                (right_source->hap1_allele ==
                    chunk.candidates[right_i].hap_to_cons_alle[1])) {
                const hts_pos_t source_ps = left_source->phase_set;
                const auto source_path =
                    gc.recovery_source_path_supported.find(source_ps);
                if (source_path != gc.recovery_source_path_supported.end() &&
                    source_path->second &&
                    (gc.recovery_source_weak_cuts.count(source_ps) == 0 ||
                     gc.recovery_source_weak_cuts.at(source_ps).empty()) &&
                    (gc.recovery_source_quality_cuts.count(source_ps) == 0 ||
                     gc.recovery_source_quality_cuts.at(source_ps).empty())) {
                    supported_edge = true;
                    continue;
                }
            }
        }
        // An independent complete BAM source may fill a missing graph edge
        // only after a two-haplotype graph edge anchors the source to this
        // block. A graph reversal or one-haplotype edge remains a veto.
        if (allow_complete_bam_source && supported_edge &&
            ((votes[0] >= kMinSourceGraphAnchorVotes &&
              votes[1] >= kMinSourceGraphAnchorVotes &&
              agree == votes[2]) ||
             (agree == 0 && votes[2] == 0)))
            continue;
        if (votes[0] >= kMinSourceGraphAnchorVotes &&
            votes[1] >= kMinSourceGraphAnchorVotes && agree > votes[2])
            supported_edge = true;
        if (votes[0] == 0 || votes[1] == 0 ||
            votes[0] + votes[1] <= votes[2]) {
            // One weak site can be bypassed only when its flanking graph
            // SNPs have a significant direct edge on both haplotypes. Never
            // bypass an edge whose reversal votes dominate agreement.
            if (allow_single_weak_site && !skipped_weak_site &&
                (votes[0] == 0 || votes[1] == 0) &&
                votes[0] + votes[1] > votes[2] &&
                si + 1 < sites.size() &&
                right_pos < sites[si + 1].first) {
                const std::array<int, 3> bypass =
                    pair_votes(left_i, sites[si + 1].second);
                if (bypass[0] >= 2 && bypass[1] >= 2 &&
                    source_graph_vote_supported(
                        bypass[0] + bypass[1], bypass[2], 0)) {
                    skipped_weak_site = true;
                    ++si;
                    continue;
                }
            }
            const bool no_graph_edge =
                votes[0] == 0 && votes[1] == 0 && votes[2] == 0;
            // The main BAM stitch may already have joined two original
            // graph blocks. Reuse that certified relation when GAF supplies
            // no contradictory vote, while checking the SNP path inside
            // each original block as usual.
            if (votes[2] == 0 && original_graph_phase_sets != nullptr &&
                stitched_graph_phase_sets != nullptr &&
                left_i < gc.site_ids.size() &&
                right_i < gc.site_ids.size()) {
                const std::string& left_id = gc.site_ids[left_i];
                const std::string& right_id = gc.site_ids[right_i];
                const auto left_original =
                    original_graph_phase_sets->find(left_id);
                const auto right_original =
                    original_graph_phase_sets->find(right_id);
                const auto left_stitched =
                    stitched_graph_phase_sets->find(left_id);
                const auto right_stitched =
                    stitched_graph_phase_sets->find(right_id);
                if (left_original != original_graph_phase_sets->end() &&
                    right_original != original_graph_phase_sets->end() &&
                    left_stitched != stitched_graph_phase_sets->end() &&
                    right_stitched != stitched_graph_phase_sets->end() &&
                    left_original->second > 0 && right_original->second > 0 &&
                    left_original->second != right_original->second &&
                    left_stitched->second > 0 &&
                    left_stitched->second == right_stitched->second)
                    continue;
            }
            if (disconnected_pair != nullptr &&
                (no_graph_edge || (report_one_hap_cut && agree > votes[2] &&
                                  (votes[0] == 0 || votes[1] == 0))))
                *disconnected_pair = {left_pos, right_pos};
            return false;
        }
    }
    return true;
}

// A weak BAM source may have an intact run after its last cut. Its clean SNP
// can bridge to the next graph block when both original graph blocks have
// continuous SNP paths and shared-source read gauges agree with that run.
static void stitch_bam_snp_runs_to_graph(GraphChunkBuildResult& gc) {
    PhasingChunk& chunk = gc.chunk;
    constexpr int kMinBridgeMapq = 30;
    constexpr double kMaxBridgeP = 0.01;
    constexpr size_t kMinSharedSitesPerFlank = 2;
    for (const RecoverySeam& seam : gc.recovery_windows) {
        std::optional<size_t> left_i;
        std::optional<size_t> right_i;
        for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
            const CandidateVariant& candidate = chunk.candidates[ci];
            const hts_pos_t pos = candidate.key.sort_pos();
            if (pos < seam.beg || pos > seam.end ||
                candidate.phase_set <= 0 ||
                candidate.hap_to_cons_alle[1] < 0 ||
                candidate.hap_to_cons_alle[2] < 0 ||
                candidate.hap_to_cons_alle[1] ==
                    candidate.hap_to_cons_alle[2])
                continue;
            if (candidate.phase_set == seam.left_phase_set)
                left_i = ci;
            else if (candidate.phase_set == seam.right_phase_set && !right_i)
                right_i = ci;
        }
        if (!left_i || !right_i || *left_i >= *right_i) continue;
        const CandidateVariant& left = chunk.candidates[*left_i];
        const CandidateVariant& right = chunk.candidates[*right_i];
        if (!left.bam_injected || right.bam_injected ||
            left.counts.category != VariantCategory::CleanHetSnp ||
            right.counts.category != VariantCategory::CleanHetSnp ||
            !graph_snp_path_supported(gc, seam.left_phase_set) ||
            !graph_snp_path_supported(gc, seam.right_phase_set))
            continue;
        const RecoverySourceSite* boundary_source = nullptr;
        for (const RecoverySourceSite& site : gc.recovery_source_sites)
            if (site.candidate_index == *left_i && site.can_adopt)
                boundary_source = &site;
        if (!boundary_source) continue;
        const hts_pos_t source_ps = boundary_source->phase_set;
        const auto cuts = gc.recovery_source_weak_cuts.find(source_ps);
        if (cuts == gc.recovery_source_weak_cuts.end() ||
            cuts->second.empty())
            continue;
        const hts_pos_t last_cut = cuts->second.back();
        const hts_pos_t left_pos = left.key.sort_pos();
        const hts_pos_t right_pos = right.key.sort_pos();
        if (last_cut >= left_pos || left_pos >= right_pos) continue;
        const int boundary_parity =
            left.hap_to_cons_alle[1] == boundary_source->hap1_allele
                ? 0 : left.hap_to_cons_alle[2] == boundary_source->hap1_allele
                    ? 1 : -1;
        if (boundary_parity < 0) continue;
        std::array<size_t, 2> shared_sites{};
        std::array<int, 2> flank_parity{{-1, -1}};
        bool conflict = false;
        for (const RecoverySourceSite& site : gc.recovery_source_sites) {
            if (site.phase_set != source_ps || !site.can_adopt ||
                site.candidate_index >= chunk.candidates.size() ||
                (site.graph_phase_set != seam.left_phase_set &&
                 site.graph_phase_set != seam.right_phase_set))
                continue;
            const CandidateVariant& candidate =
                chunk.candidates[site.candidate_index];
            const hts_pos_t pos = candidate.key.sort_pos();
            if (candidate.bam_injected || pos <= last_cut ||
                (site.graph_phase_set == seam.left_phase_set &&
                 pos > left_pos) ||
                candidate.counts.category != VariantCategory::CleanHetSnp ||
                candidate.phase_set != site.graph_phase_set)
                continue;
            const int parity =
                candidate.hap_to_cons_alle[1] == site.hap1_allele
                    ? 0 : candidate.hap_to_cons_alle[2] == site.hap1_allele
                        ? 1 : -1;
            const size_t flank =
                site.graph_phase_set == seam.left_phase_set ? 0 : 1;
            if (parity < 0 ||
                (flank_parity[flank] >= 0 && flank_parity[flank] != parity)) {
                conflict = true;
                break;
            }
            flank_parity[flank] = parity;
            ++shared_sites[flank];
        }
        if (conflict || shared_sites[0] < kMinSharedSitesPerFlank ||
            shared_sites[1] < kMinSharedSitesPerFlank ||
            flank_parity[0] != boundary_parity)
            continue;
        bool gauges_supported = false;
        for (const RecoveryPhaseGauge& gauge : gc.recovery_phase_gauges) {
            std::array<bool, 2> approved{};
            for (const RecoveryBlockGaugeVote& vote : gauge.block_votes) {
                if (vote.bam_phase_set != source_ps) continue;
                const size_t flank =
                    vote.graph_phase_set == seam.left_phase_set ? 0 :
                    vote.graph_phase_set == seam.right_phase_set ? 1 : 2;
                if (flank > 1) continue;
                const int same = vote.counts[0][0] + vote.counts[1][1];
                const int cross = vote.counts[0][1] + vote.counts[1][0];
                const bool both_haps =
                    vote.counts[0][0] + vote.counts[0][1] > 0 &&
                    vote.counts[1][0] + vote.counts[1][1] > 0;
                approved[flank] = both_haps &&
                    source_graph_vote_supported(
                        same, cross, flank_parity[flank]);
            }
            if (approved[0] && approved[1]) {
                gauges_supported = true;
                break;
            }
        }
        if (!gauges_supported) continue;
        const std::optional<bool> flip = local_run_boundary_flip(
            chunk, *left_i, *right_i, kMinBridgeMapq, kMaxBridgeP);
        if (!flip || *flip !=
                (flank_parity[0] != flank_parity[1]))
            continue;
        merge_phase_sets_in_place(
            chunk, seam.left_phase_set, seam.right_phase_set, *flip);
    }
}

// Join the boundary-side prefix only when the later graph path has an empty
// cut. An intermediate candidate or any crossing read would make that split
// discard an already established link, so both are hard vetoes.
static bool join_graph_prefix_at_disconnected_cut(
        GraphChunkBuildResult& gc, WorkerContext& context, int tid,
        hts_pos_t left_ps, hts_pos_t right_ps, bool flip,
        hts_pos_t boundary_pos,
        const std::pair<hts_pos_t, hts_pos_t>& cut) {
    if (cut.first < boundary_pos || cut.first >= cut.second) return false;
    std::unique_ptr<hts_itr_t, decltype(&hts_itr_destroy)> iterator(
        sam_itr_queryi(context.indexes.front().get(), tid,
                       cut.first - 1, cut.second), &hts_itr_destroy);
    if (!iterator) return false;
    std::unique_ptr<bam1_t, AlignmentDeleter> alignment(bam_init1());
    if (!alignment) return false;
    while (sam_itr_next(context.bams.front()->get(), iterator.get(),
                        alignment.get()) >= 0) {
        const bam1_t* read = alignment.get();
        if (read->core.flag & (BAM_FUNMAP | BAM_FSECONDARY |
                               BAM_FSUPPLEMENTARY))
            continue;
        if (read->core.pos <= cut.first - 1 &&
            bam_endpos(read) > cut.second - 1)
            return false;
    }

    PhasingChunk& chunk = gc.chunk;
    std::vector<bool> move(chunk.candidates.size(), false);
    bool has_prefix = false, has_tail = false;
    for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
        const CandidateVariant& candidate = chunk.candidates[ci];
        if (candidate.phase_set != right_ps) continue;
        hts_pos_t physical_pos = candidate.key.sort_pos();
        if (!candidate.bam_injected && ci < gc.site_meta.size()) {
            const std::string* alt = selected_graph_candidate_alt(gc, ci);
            if (alt == nullptr) return false;
            const GraphSiteMeta& meta = gc.site_meta[ci];
            physical_pos = vcf_to_variant_key(
                candidate.key.tid, meta.pos, meta.ref, *alt).pos;
        }
        // The cut must be empty of right-block candidates, including indels.
        if (physical_pos > cut.first && physical_pos < cut.second)
            return false;
        if (physical_pos <= cut.first) {
            move[ci] = true;
            has_prefix = true;
        } else {
            has_tail = true;
        }
    }
    if (!has_prefix || !has_tail) return false;
    std::vector<size_t> reads_to_move;
    for (size_t ri = 0; ri < chunk.reads.size(); ++ri) {
        if (ri >= chunk.read_var_profile.size()) continue;
        const ReadVariantProfile& profile = chunk.read_var_profile[ri];
        if (profile.start_var_idx < 0) continue;
        bool sees_prefix = false, sees_tail = false;
        for (size_t offset = 0; offset < profile.alleles.size(); ++offset) {
            const size_t ci =
                static_cast<size_t>(profile.start_var_idx) + offset;
            if (ci >= chunk.candidates.size()) break;
            if (chunk.candidates[ci].phase_set != right_ps) continue;
            const bool called = profile.alleles[offset] >= 0 ||
                (offset < profile.graph_alleles.size() &&
                 profile.graph_alleles[offset] >= 0) ||
                (offset < profile.bam_alleles.size() &&
                 profile.bam_alleles[offset] >= 0);
            if (!called) continue;
            if (move[ci]) sees_prefix = true;
            else sees_tail = true;
        }
        if (sees_prefix && sees_tail) return false;
        if (sees_prefix && ri < chunk.phase_sets.size() &&
            chunk.phase_sets[ri] == right_ps)
            reads_to_move.push_back(ri);
    }
    const auto flip_hap = [](int hap) {
        return hap == 1 ? 2 : hap == 2 ? 1 : hap;
    };
    for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
        if (!move[ci]) continue;
        CandidateVariant& candidate = chunk.candidates[ci];
        if (flip) {
            std::swap(candidate.hap_to_cons_alle[1],
                      candidate.hap_to_cons_alle[2]);
            std::swap(candidate.hap_to_alle_profile[1],
                      candidate.hap_to_alle_profile[2]);
            candidate.hap_alt = flip_hap(candidate.hap_alt);
            candidate.hap_ref = flip_hap(candidate.hap_ref);
        }
        candidate.phase_set = left_ps;
    }
    for (const size_t ri : reads_to_move) {
        chunk.phase_sets[ri] = left_ps;
        if (flip && ri < chunk.haps.size())
            chunk.haps[ri] = flip_hap(chunk.haps[ri]);
    }
    return true;
}

// A physical allele bridge may certify a suffix even when an earlier
// graph or BAM source edge is unsupported. Split at the least traversed
// candidate boundary inside that edge, and leave cross-boundary reads unphased:
// one read cannot carry two independent phase-set tags.
static bool join_graph_suffix_after_weak_edge(
        GraphChunkBuildResult& gc, WorkerContext& context, int tid,
        hts_pos_t left_ps, hts_pos_t right_ps, bool flip,
        size_t source_i, const std::pair<hts_pos_t, hts_pos_t>& cut,
        const std::map<std::string, hts_pos_t>* original_graph_phase_sets,
        const std::map<std::string, hts_pos_t>* stitched_graph_phase_sets,
        bool certified_suffix = false) {
    if (cut.first <= 0 || cut.first >= cut.second) return false;
    std::unique_ptr<hts_itr_t, decltype(&hts_itr_destroy)> iterator(
        sam_itr_queryi(context.indexes.front().get(), tid,
                       cut.first - 1, cut.second), &hts_itr_destroy);
    if (!iterator) return false;
    std::unique_ptr<bam1_t, AlignmentDeleter> alignment(bam_init1());
    if (!alignment) return false;
    // A molecule spanning both graph SNPs would provide a direct edge; in
    // that case this is not an evidence-free boundary to split.
    while (sam_itr_next(context.bams.front()->get(), iterator.get(),
                         alignment.get()) >= 0) {
        const bam1_t* read = alignment.get();
        if (read->core.flag & (BAM_FUNMAP | BAM_FSECONDARY |
                               BAM_FSUPPLEMENTARY))
            continue;
        if (read->core.pos <= cut.first - 1 &&
            bam_endpos(read) > cut.second - 1 && !certified_suffix)
            return false;
    }
    if (!graph_snp_path_supported(
            gc, left_ps, nullptr, original_graph_phase_sets,
            stitched_graph_phase_sets, false, false, cut.second) &&
        !certified_suffix)
        return false;

    PhasingChunk& chunk = gc.chunk;
    std::vector<hts_pos_t> positions(chunk.candidates.size(), 0);
    std::vector<hts_pos_t> ends(chunk.candidates.size(), 0);
    std::vector<hts_pos_t> boundaries{cut.second};
    for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
        const CandidateVariant& candidate = chunk.candidates[ci];
        if (candidate.phase_set != left_ps) continue;
        VariantKey physical = candidate.key;
        if (!candidate.bam_injected && ci < gc.site_meta.size()) {
            const std::string* alt = selected_graph_candidate_alt(gc, ci);
            if (alt == nullptr) return false;
            const GraphSiteMeta& meta = gc.site_meta[ci];
            physical = vcf_to_variant_key(
                candidate.key.tid, meta.pos, meta.ref, *alt);
        }
        positions[ci] = physical.sort_pos();
        ends[ci] = physical.pos + std::max(1, physical.ref_len) - 1;
        if (positions[ci] > cut.first && positions[ci] < cut.second)
            boundaries.push_back(positions[ci]);
    }
    if (source_i >= positions.size() || positions[source_i] < cut.second)
        return false;

    // Each tagged read is represented by its first and last called sites.
    // Counting ranges crossing each candidate boundary finds the split that
    // discards the fewest established read assignments.
    std::vector<std::pair<hts_pos_t, hts_pos_t>> read_ranges(
        chunk.reads.size(), {0, 0});
    for (size_t ri = 0; ri < chunk.reads.size(); ++ri) {
        if (ri >= chunk.read_var_profile.size() ||
            ri >= chunk.phase_sets.size() || chunk.phase_sets[ri] != left_ps)
            continue;
        const ReadVariantProfile& profile = chunk.read_var_profile[ri];
        if (profile.start_var_idx < 0) continue;
        hts_pos_t first = std::numeric_limits<hts_pos_t>::max();
        hts_pos_t last = 0;
        for (size_t offset = 0; offset < profile.alleles.size(); ++offset) {
            const size_t ci = static_cast<size_t>(profile.start_var_idx) + offset;
            if (ci >= positions.size()) break;
            if (chunk.candidates[ci].phase_set != left_ps) continue;
            const bool called = profile.alleles[offset] >= 0 ||
                (offset < profile.graph_alleles.size() &&
                 profile.graph_alleles[offset] >= 0) ||
                (offset < profile.bam_alleles.size() &&
                 profile.bam_alleles[offset] >= 0);
            if (!called) continue;
            first = std::min(first, positions[ci]);
            last = std::max(last, positions[ci]);
        }
        if (last > 0) read_ranges[ri] = {first, last};
    }
    size_t fewest_crossing = std::numeric_limits<size_t>::max();
    hts_pos_t boundary = 0;
    for (const hts_pos_t candidate_boundary : boundaries) {
        bool splits_candidate = false;
        bool has_prefix = false, has_suffix = false;
        for (size_t ci = 0; ci < positions.size(); ++ci) {
            if (chunk.candidates[ci].phase_set != left_ps) continue;
            if (positions[ci] < candidate_boundary) {
                has_prefix = true;
                if (ends[ci] >= candidate_boundary) splits_candidate = true;
            } else {
                has_suffix = true;
            }
        }
        if (splits_candidate || !has_prefix || !has_suffix) continue;
        size_t crossing = 0;
        for (const auto& [first, last] : read_ranges)
            if (first > 0 && first < candidate_boundary &&
                last >= candidate_boundary)
                ++crossing;
        if (crossing < fewest_crossing) {
            fewest_crossing = crossing;
            boundary = candidate_boundary;
        }
    }
    if (boundary == 0) return false;
    const auto flip_hap = [](int hap) {
        return hap == 1 ? 2 : hap == 2 ? 1 : hap;
    };
    for (size_t ci = 0; ci < positions.size(); ++ci) {
        if (chunk.candidates[ci].phase_set != left_ps ||
            positions[ci] < boundary)
            continue;
        CandidateVariant& candidate = chunk.candidates[ci];
        if (flip) {
            std::swap(candidate.hap_to_cons_alle[1],
                      candidate.hap_to_cons_alle[2]);
            std::swap(candidate.hap_to_alle_profile[1],
                      candidate.hap_to_alle_profile[2]);
            candidate.hap_alt = flip_hap(candidate.hap_alt);
            candidate.hap_ref = flip_hap(candidate.hap_ref);
        }
        candidate.phase_set = right_ps;
    }
    for (size_t ri = 0; ri < read_ranges.size(); ++ri) {
        const auto [first, last] = read_ranges[ri];
        if (first == 0 || last < boundary) continue;
        if (first < boundary) {
            chunk.phase_sets[ri] = 0;
            if (ri < chunk.haps.size()) chunk.haps[ri] = 0;
        } else {
            chunk.phase_sets[ri] = right_ps;
            if (flip && ri < chunk.haps.size())
                chunk.haps[ri] = flip_hap(chunk.haps[ri]);
        }
    }
    return true;
}

// Call an entire biallelic substitution from aligned bases. A multi-base
// graph allele is one site: a read carrying a mixed REF/ALT sequence cannot
// vote for either haplotype. The minimum quality covers every changed base.
static int physical_substitution_call(
        const bam1_t* read, hts_pos_t pos, const std::string& ref,
        const std::string& alt, int* min_quality) {
    if (ref.empty() || ref.size() != alt.size()) return -1;
    int allele = -1;
    int quality = std::numeric_limits<int>::max();
    for (size_t offset = 0; offset < ref.size(); ++offset) {
        int base_quality = 0;
        const int call = physical_snp_call(
            read, pos + static_cast<hts_pos_t>(offset), ref[offset],
            alt[offset], &base_quality);
        if (call != 0 && call != 2) return -1;
        quality = std::min(quality, base_quality);
        if (std::toupper(static_cast<unsigned char>(ref[offset])) ==
            std::toupper(static_cast<unsigned char>(alt[offset])))
            continue;
        if (allele >= 0 && allele != call) return -1;
        allele = call;
    }
    if (allele < 0) return -1;
    if (min_quality != nullptr) *min_quality = quality;
    return allele;
}

// Two separate BAM deletion rows can describe the two alleles at an
// overlapping locus. Keep both rows and use their exact CIGAR ALT calls to
// orient the next graph SNP; a REF call at either row is ambiguous when the
// other deletion removes its reference span.
static int physical_equivalent_deletion_call(
        const bam1_t* read, const CandidateVariant& deletion,
        WorkerContext& context, int tid, int min_baseq);

static bool stitch_complementary_deletions_to_snp(
        GraphChunkBuildResult& gc, WorkerContext& context, int tid,
        const RecoverySeam& seam, hts_pos_t left_ps, hts_pos_t right_ps,
        hts_pos_t right_pos, const std::string& right_ref,
        const std::string& right_alt,
        size_t right_i,
        const std::map<std::string, hts_pos_t>* original_graph_phase_sets,
        const std::map<std::string, hts_pos_t>* stitched_graph_phase_sets) {
    constexpr int kMinMapq = 30;
    constexpr int kMinBaseq = 30;
    constexpr int kUnknownQuality = 255;
    constexpr double kMaxWrongParity = 0.001;
    PhasingChunk& chunk = gc.chunk;
    std::optional<std::pair<size_t, size_t>> pair;
    hts_pos_t pair_end = 0;
    for (size_t i = 0; i < gc.recovery_source_sites.size(); ++i) {
        const RecoverySourceSite& first = gc.recovery_source_sites[i];
        if (first.candidate_index >= chunk.candidates.size()) continue;
        const CandidateVariant& a = chunk.candidates[first.candidate_index];
        if (a.phase_set != left_ps || a.key.type != VariantType::Deletion ||
            (a.counts.category != VariantCategory::CleanHetIndel &&
             a.counts.category != VariantCategory::NoisyCandHet) ||
            !a.msa_verified || a.hap_to_cons_alle[1] < 0 ||
            a.hap_to_cons_alle[1] > 1 ||
            a.hap_to_cons_alle[2] != 1 - a.hap_to_cons_alle[1] ||
            a.key.sort_pos() < seam.beg || a.key.pos >= right_pos)
            continue;
        for (size_t j = i + 1; j < gc.recovery_source_sites.size(); ++j) {
            const RecoverySourceSite& second = gc.recovery_source_sites[j];
            if (second.phase_set != first.phase_set ||
                second.candidate_index >= chunk.candidates.size())
                continue;
            const CandidateVariant& b =
                chunk.candidates[second.candidate_index];
            if (b.phase_set != left_ps ||
                b.key.type != VariantType::Deletion ||
                (b.counts.category != VariantCategory::CleanHetIndel &&
                 b.counts.category != VariantCategory::NoisyCandHet) ||
                !b.msa_verified || b.hap_to_cons_alle[1] < 0 ||
                b.hap_to_cons_alle[1] > 1 ||
                b.hap_to_cons_alle[2] != 1 - b.hap_to_cons_alle[1] ||
                a.hap_to_cons_alle[1] == b.hap_to_cons_alle[1] ||
                std::max(a.key.pos, b.key.pos) >=
                    std::min(a.key.pos + a.key.ref_len,
                             b.key.pos + b.key.ref_len))
                continue;
            const hts_pos_t end = std::max(a.key.pos + a.key.ref_len,
                                           b.key.pos + b.key.ref_len);
            if (end >= right_pos || end <= pair_end) continue;
            pair = {first.candidate_index, second.candidate_index};
            pair_end = end;
        }
    }
    if (!pair ||
        !graph_snp_path_supported(gc, left_ps, nullptr,
                                  original_graph_phase_sets,
                                  stitched_graph_phase_sets) ||
        !graph_snp_path_supported(gc, right_ps, nullptr,
                                  original_graph_phase_sets,
                                  stitched_graph_phase_sets))
        return false;
    const CandidateVariant& a = chunk.candidates[pair->first];
    const CandidateVariant& b = chunk.candidates[pair->second];
    std::unique_ptr<hts_itr_t, decltype(&hts_itr_destroy)> iterator(
        sam_itr_queryi(context.indexes.front().get(), tid,
                       std::min(a.key.pos, b.key.pos) - 1, right_pos),
        &hts_itr_destroy);
    if (!iterator) return false;
    std::unique_ptr<bam1_t, AlignmentDeleter> alignment(bam_init1());
    if (!alignment) return false;
    std::unordered_set<std::string> seen;
    std::array<int, 2> allele_support{};
    std::array<int, 2> relation_support{};
    double log_odds = 0.0;
    while (sam_itr_next(context.bams.front()->get(), iterator.get(),
                        alignment.get()) >= 0) {
        const bam1_t* read = alignment.get();
        if ((read->core.flag & (BAM_FUNMAP | BAM_FSECONDARY |
                                BAM_FSUPPLEMENTARY | BAM_FDUP | BAM_FQCFAIL)) ||
            read->core.qual < kMinMapq ||
            read->core.qual == kUnknownQuality ||
            seen.count(bam_get_qname(read)) != 0)
            continue;
        int right_quality = 0;
        const int right_call = physical_substitution_call(
            read, right_pos, right_ref, right_alt, &right_quality);
        if ((right_call != 0 && right_call != 2) ||
            right_quality < kMinBaseq || right_quality == kUnknownQuality)
            continue;
        // Score only a verified ALT of exactly one row. Shift equivalence
        // recovers the same deletion; another length never supplies its REF.
        const int a_call = physical_equivalent_deletion_call(
            read, a, context, tid, kMinBaseq);
        const int b_call = physical_equivalent_deletion_call(
            read, b, context, tid, kMinBaseq);
        if ((a_call == 1) == (b_call == 1)) continue;
        const CandidateVariant& deletion = a_call == 1 ? a : b;
        // Equivalence already verifies the actual placement's Q30 flanks.
        // A flank of the original key may lie inside the shifted deletion.
        seen.insert(bam_get_qname(read));
        ++allele_support[a_call == 1 ? 0 : 1];
        const bool left_hap1 = deletion.hap_to_cons_alle[1] == 1;
        const bool right_hap1 = (right_call == 2) ==
            (chunk.candidates[right_i].hap_to_cons_alle[1] == 1);
        const bool flip = left_hap1 != right_hap1;
        ++relation_support[flip ? 1 : 0];
        const double p = std::pow(10.0, -kMinBaseq / 10.0) +
            std::pow(10.0, -right_quality / 10.0) +
            2.0 * std::pow(10.0, -read->core.qual / 10.0);
        if (p > 0.0 && p < 0.5)
            log_odds += (flip ? 1.0 : -1.0) *
                std::log((1.0 - p) / p);
    }
    const double threshold = std::log((1.0 - kMaxWrongParity) /
                                      kMaxWrongParity);
    if (allele_support[0] == 0 || allele_support[1] == 0 ||
        (relation_support[0] > 0 && relation_support[1] > 0) ||
        std::abs(log_odds) < threshold)
        return false;
    return merge_phase_sets_in_place(chunk, left_ps, right_ps,
                                     log_odds > 0.0);
}

// CIGAR can place the same deletion at different bases in a short repeat.
// Compare the resulting reference strings; retain the candidate's original
// coordinates and never turn a nearby, different indel into a REF call.
static int physical_equivalent_deletion_call(
        const bam1_t* read, const CandidateVariant& deletion,
        WorkerContext& context, int tid, int min_baseq) {
    return bam_equivalent_deletion_allele(
        read, deletion, context.ref, tid, context.primary_header(), min_baseq);
}

// Verify that a boundary deletion retains the orientation of its nearest
// clean graph SNP. A weak, distant graph SNP edge cannot certify this local
// link, but independent Q30 primary calls can.
static std::optional<hts_pos_t> graph_snp_deletion_link_supported(
        const GraphChunkBuildResult& gc, WorkerContext& context, int tid,
        hts_pos_t phase_set, const CandidateVariant& deletion,
        bool snp_on_left) {
    constexpr int kMinMapq = 30;
    constexpr int kMinBaseq = 30;
    constexpr int kUnknownQuality = 255;
    constexpr int kMinPairedReads = 2;
    constexpr double kMaxWrongParity = 0.001;
    const PhasingChunk& chunk = gc.chunk;
    std::optional<size_t> nearest_i;
    VariantKey nearest_key;
    for (size_t ci = 0; ci < chunk.candidates.size() &&
                        ci < gc.site_meta.size(); ++ci) {
        const CandidateVariant& candidate = chunk.candidates[ci];
        if (candidate.phase_set != phase_set || candidate.bam_injected ||
            candidate.counts.category != VariantCategory::CleanHetSnp ||
            candidate.hap_to_cons_alle[1] < 0 ||
            candidate.hap_to_cons_alle[1] > 1 ||
            candidate.hap_to_cons_alle[2] !=
                1 - candidate.hap_to_cons_alle[1])
            continue;
        const std::string* alt = selected_graph_candidate_alt(gc, ci);
        if (alt == nullptr) continue;
        const GraphSiteMeta& meta = gc.site_meta[ci];
        const VariantKey key = vcf_to_variant_key(
            candidate.key.tid, meta.pos, meta.ref, *alt);
        if (key.type != VariantType::Snp || key.alt.size() != 1 ||
            (snp_on_left ? key.pos >= deletion.key.pos :
                           key.pos <= deletion.key.pos))
            continue;
        if (!nearest_i || (snp_on_left ?
                key.pos > nearest_key.pos : key.pos < nearest_key.pos)) {
            nearest_i = ci;
            nearest_key = key;
        }
    }
    if (!nearest_i) return std::nullopt;
    const char ref = context.ref.base(
        tid, nearest_key.pos, context.primary_header());
    if (ref == 'N') return std::nullopt;
    std::unique_ptr<hts_itr_t, decltype(&hts_itr_destroy)> iterator(
        sam_itr_queryi(context.indexes.front().get(), tid,
                       std::min(nearest_key.pos, deletion.key.pos) - 1,
                       std::max(nearest_key.pos, deletion.key.pos +
                                              deletion.key.ref_len)),
        &hts_itr_destroy);
    if (!iterator) return std::nullopt;
    std::unique_ptr<bam1_t, AlignmentDeleter> alignment(bam_init1());
    if (!alignment) return std::nullopt;
    std::unordered_set<std::string> seen;
    int paired = 0;
    double log_odds = 0.0;
    while (sam_itr_next(context.bams.front()->get(), iterator.get(),
                        alignment.get()) >= 0) {
        const bam1_t* read = alignment.get();
        if ((read->core.flag & (BAM_FUNMAP | BAM_FSECONDARY |
                                BAM_FSUPPLEMENTARY | BAM_FDUP | BAM_FQCFAIL)) ||
            read->core.qual < kMinMapq ||
            read->core.qual == kUnknownQuality ||
            seen.count(bam_get_qname(read)) != 0)
            continue;
        int snp_quality = 0;
        const int snp_call = physical_snp_call(
            read, nearest_key.pos, ref, nearest_key.alt[0], &snp_quality);
        const int deletion_call = physical_equivalent_deletion_call(
            read, deletion, context, tid, kMinBaseq);
        if ((snp_call != 0 && snp_call != 2) ||
            deletion_call < 0 || snp_quality < kMinBaseq ||
            snp_quality == kUnknownQuality)
            continue;
        seen.insert(bam_get_qname(read));
        const bool snp_hap1 = (snp_call == 2) ==
            (chunk.candidates[*nearest_i].hap_to_cons_alle[1] == 1);
        const bool deletion_hap1 = (deletion_call == 1) ==
            (deletion.hap_to_cons_alle[1] == 1);
        if (snp_hap1 != deletion_hap1) return std::nullopt;
        const double p = std::pow(10.0, -snp_quality / 10.0) +
            std::pow(10.0, -kMinBaseq / 10.0) +
            2.0 * std::pow(10.0, -read->core.qual / 10.0);
        if (p <= 0.0 || p >= 0.5) continue;
        log_odds += std::log((1.0 - p) / p);
        ++paired;
    }
    const double threshold = std::log((1.0 - kMaxWrongParity) /
                                      kMaxWrongParity);
    return paired >= kMinPairedReads && log_odds >= threshold ?
        std::optional<hts_pos_t>(nearest_key.pos) : std::nullopt;
}

// Exact clean graph deletion boundaries can be linked by the same primary
// molecules even when neither side has a callable SNP at the seam. A weak
// upstream graph edge leaves the left prefix separate from the joined suffix.
static bool stitch_graph_deletion_pair(
        GraphChunkBuildResult& gc, WorkerContext& context, int tid,
        const RecoverySeam& seam,
        const std::map<std::string, hts_pos_t>* original_graph_phase_sets,
        const std::map<std::string, hts_pos_t>* stitched_graph_phase_sets) {
    constexpr int kMinMapq = 30;
    constexpr int kMinBaseq = 30;
    constexpr int kUnknownQuality = 255;
    constexpr int kMinPairedReads = 3;
    PhasingChunk& chunk = gc.chunk;
    std::array<std::optional<size_t>, 2> boundary;
    std::array<VariantKey, 2> physical_keys;
    for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
        const CandidateVariant& candidate = chunk.candidates[ci];
        const int hap1 = candidate.hap_to_cons_alle[1];
        if (candidate.phase_set <= 0 || hap1 < 0 || hap1 > 1 ||
            candidate.hap_to_cons_alle[2] != 1 - hap1)
            continue;
        VariantKey physical = candidate.key;
        if (ci < gc.site_meta.size() && !candidate.bam_injected) {
            const std::string* alt = selected_graph_candidate_alt(gc, ci);
            if (alt == nullptr) continue;
            const GraphSiteMeta& meta = gc.site_meta[ci];
            physical = vcf_to_variant_key(candidate.key.tid,
                                          meta.pos, meta.ref, *alt);
        }
        const int side = candidate.phase_set == seam.left_phase_set &&
                physical.sort_pos() == seam.beg ? 0 :
            candidate.phase_set == seam.right_phase_set &&
                physical.sort_pos() == seam.end ? 1 : -1;
        if (side < 0) continue;
        // Another phased row at the same locus would make the REF call
        // ambiguous, so neither allele can certify a whole-block join.
        if (boundary[side]) return false;
        if (ci >= gc.site_meta.size() || candidate.bam_injected ||
            physical.type != VariantType::Deletion ||
            candidate.counts.category != VariantCategory::CleanHetIndel)
            return false;
        boundary[side] = ci;
        physical_keys[side] = std::move(physical);
    }
    if (!boundary[0] || !boundary[1]) return false;
    std::pair<hts_pos_t, hts_pos_t> left_cut, right_cut;
    const bool left_path = graph_snp_path_supported(
        gc, seam.left_phase_set, &left_cut, original_graph_phase_sets,
        stitched_graph_phase_sets, false, false, 0, true);
    const bool right_path = graph_snp_path_supported(
        gc, seam.right_phase_set, &right_cut, original_graph_phase_sets,
        stitched_graph_phase_sets, false, false, 0, true);
    CandidateVariant left = chunk.candidates[*boundary[0]];
    CandidateVariant right = chunk.candidates[*boundary[1]];
    left.key = std::move(physical_keys[0]);
    right.key = std::move(physical_keys[1]);
    std::unique_ptr<hts_itr_t, decltype(&hts_itr_destroy)> iterator(
        sam_itr_queryi(context.indexes.front().get(), tid,
                       left.key.pos - 1, right.key.pos), &hts_itr_destroy);
    if (!iterator) return false;
    std::unique_ptr<bam1_t, AlignmentDeleter> alignment(bam_init1());
    if (!alignment) return false;
    std::unordered_set<std::string> seen;
    std::array<int, 2> parity_votes{};
    std::array<std::array<int, 2>, 2> allele_support{};
    while (sam_itr_next(context.bams.front()->get(), iterator.get(),
                        alignment.get()) >= 0) {
        const bam1_t* read = alignment.get();
        if ((read->core.flag & (BAM_FUNMAP | BAM_FSECONDARY |
                                BAM_FSUPPLEMENTARY | BAM_FDUP | BAM_FQCFAIL)) ||
            read->core.qual < kMinMapq ||
            read->core.qual == kUnknownQuality ||
            seen.count(bam_get_qname(read)) != 0)
            continue;
        const int left_call = physical_equivalent_deletion_call(
            read, left, context, tid, kMinBaseq);
        const int right_call = physical_equivalent_deletion_call(
            read, right, context, tid, kMinBaseq);
        if (left_call < 0 || right_call < 0) continue;
        seen.insert(bam_get_qname(read));
        ++allele_support[0][static_cast<size_t>(left_call)];
        ++allele_support[1][static_cast<size_t>(right_call)];
        const bool left_hap1 =
            (left_call == 1) == (left.hap_to_cons_alle[1] == 1);
        const bool right_hap1 =
            (right_call == 1) == (right.hap_to_cons_alle[1] == 1);
        ++parity_votes[left_hap1 != right_hap1 ? 1 : 0];
    }
    if (parity_votes[0] + parity_votes[1] < kMinPairedReads ||
        (parity_votes[0] > 0 && parity_votes[1] > 0) ||
        allele_support[0][0] == 0 || allele_support[0][1] == 0 ||
        allele_support[1][0] == 0 || allele_support[1][1] == 0)
        return false;
    const bool flip = parity_votes[1] > 0;
    if (left_path && right_path)
        return merge_phase_sets_in_place(chunk, seam.left_phase_set,
                                         seam.right_phase_set, flip);
    // A missing or one-haplotype graph SNP edge can precede the boundary
    // deletion. Join only its directly certified suffix. The right graph
    // block can have a later weak edge, provided its local prefix is intact
    // and the first SNP is physically linked to the deletion.
    if (left_path || left_cut.first <= 0 ||
        left_cut.second >= left.key.sort_pos() ||
        (!right_path &&
         (right_cut.first <= right.key.sort_pos() ||
          right_cut.first >= right_cut.second)))
        return false;
    const auto left_snp = graph_snp_deletion_link_supported(
        gc, context, tid, seam.left_phase_set, left, true);
    const auto right_snp = graph_snp_deletion_link_supported(
        gc, context, tid, seam.right_phase_set, right, false);
    if (!left_snp || *left_snp != left_cut.second ||
        !right_snp || (!right_path && *right_snp > right_cut.first))
        return false;
    return join_graph_suffix_after_weak_edge(
        gc, context, tid, seam.left_phase_set, seam.right_phase_set,
        flip, *boundary[0], left_cut, original_graph_phase_sets,
        stitched_graph_phase_sets, true);
}

// Use an established graph read haplotype as the left observation when its
// BAM SNP base is low quality. On the right, only an MSA-verified deletion ALT
// is informative: overlapping deletion rows make their REF calls ambiguous.
static bool stitch_graph_snp_to_complementary_deletions(
        GraphChunkBuildResult& gc, WorkerContext& context, int tid,
        const RecoverySeam& seam, hts_pos_t left_ps, hts_pos_t right_ps,
        hts_pos_t left_pos, size_t left_i,
        const std::map<std::string, hts_pos_t>* original_graph_phase_sets,
        const std::map<std::string, hts_pos_t>* stitched_graph_phase_sets,
        std::optional<bool> required_flip = std::nullopt) {
    constexpr int kMinMapq = 30;
    // Exact allele sequences can retain moderate-quality repeat bases. The
    // parity likelihood below charges every molecule the conservative Q20
    // error bound and still requires two independent unanimous witnesses.
    constexpr int kMinDeletionBaseq = 20;
    constexpr int kUnknownQuality = 255;
    constexpr int kMinPairedReads = 2;
    constexpr double kMaxWrongParity = 0.001;
    PhasingChunk& chunk = gc.chunk;
    if (chunk.candidates[left_i].bam_injected ||
        !graph_snp_path_supported(gc, left_ps, nullptr,
                                  original_graph_phase_sets,
                                  stitched_graph_phase_sets) ||
        // A conditional snarl SNP may lack one allele's observations. Only
        // a decisive direct edge between both flanking SNPs can bypass it.
        !graph_snp_path_supported(gc, right_ps, nullptr,
                                  original_graph_phase_sets,
                                  stitched_graph_phase_sets, true))
        return false;
    std::optional<std::pair<size_t, size_t>> pair;
    hts_pos_t nearest = std::numeric_limits<hts_pos_t>::max();
    for (size_t i = 0; i < gc.recovery_source_sites.size(); ++i) {
        const RecoverySourceSite& first = gc.recovery_source_sites[i];
        if (first.candidate_index >= chunk.candidates.size()) continue;
        const CandidateVariant& a = chunk.candidates[first.candidate_index];
        if (a.phase_set != right_ps || a.key.type != VariantType::Deletion ||
            (a.counts.category != VariantCategory::CleanHetIndel &&
             a.counts.category != VariantCategory::NoisyCandHet) ||
            !a.msa_verified || a.hap_to_cons_alle[1] < 0 ||
            a.hap_to_cons_alle[1] > 1 ||
            a.hap_to_cons_alle[2] != 1 - a.hap_to_cons_alle[1] ||
            a.key.sort_pos() <= left_pos || a.key.sort_pos() > seam.end)
            continue;
        for (size_t j = i + 1; j < gc.recovery_source_sites.size(); ++j) {
            const RecoverySourceSite& second = gc.recovery_source_sites[j];
            if (second.phase_set != first.phase_set ||
                second.candidate_index >= chunk.candidates.size())
                continue;
            const CandidateVariant& b =
                chunk.candidates[second.candidate_index];
            if (b.phase_set != right_ps || b.key.type != VariantType::Deletion ||
                (b.counts.category != VariantCategory::CleanHetIndel &&
                 b.counts.category != VariantCategory::NoisyCandHet) ||
                !b.msa_verified || b.hap_to_cons_alle[1] < 0 ||
                b.hap_to_cons_alle[1] > 1 ||
                b.hap_to_cons_alle[2] != 1 - b.hap_to_cons_alle[1] ||
                a.hap_to_cons_alle[1] == b.hap_to_cons_alle[1] ||
                std::max(a.key.pos, b.key.pos) >=
                    std::min(a.key.pos + a.key.ref_len,
                             b.key.pos + b.key.ref_len))
                continue;
            const hts_pos_t pos = std::min(a.key.sort_pos(),
                                           b.key.sort_pos());
            if (pos >= nearest) continue;
            nearest = pos;
            pair = {first.candidate_index, second.candidate_index};
        }
    }
    if (!pair) return false;
    const CandidateVariant& a = chunk.candidates[pair->first];
    const CandidateVariant& b = chunk.candidates[pair->second];
    std::unordered_map<std::string, std::pair<int, int>> graph_haps;
    for (size_t ri = 0; ri < chunk.read_var_profile.size() &&
                        ri < chunk.reads.size() &&
                        ri < chunk.phase_sets.size() &&
                        ri < chunk.haps.size(); ++ri) {
        const ReadRecord& read = chunk.reads[ri];
        if (read.is_skipped || read.mapq < kMinMapq ||
            read.mapq == kUnknownQuality ||
            chunk.phase_sets[ri] != left_ps ||
            (chunk.haps[ri] != 1 && chunk.haps[ri] != 2))
            continue;
        const ReadVariantProfile& profile = chunk.read_var_profile[ri];
        if (profile.start_var_idx < 0 ||
            left_i < static_cast<size_t>(profile.start_var_idx) ||
            left_i > static_cast<size_t>(profile.end_var_idx))
            continue;
        const size_t offset = left_i -
            static_cast<size_t>(profile.start_var_idx);
        if (offset >= profile.alleles.size() ||
            profile.alleles[offset] !=
                chunk.candidates[left_i].hap_to_cons_alle[chunk.haps[ri]])
            continue;
        auto [it, inserted] = graph_haps.emplace(
            read.qname, std::make_pair(chunk.haps[ri], read.mapq));
        if (!inserted && it->second.first != chunk.haps[ri])
            it->second.first = 0;
    }
    if (graph_haps.empty()) return false;
    std::unique_ptr<hts_itr_t, decltype(&hts_itr_destroy)> iterator(
        sam_itr_queryi(context.indexes.front().get(), tid, left_pos - 1,
                       std::max(a.key.pos + a.key.ref_len,
                                b.key.pos + b.key.ref_len)),
        &hts_itr_destroy);
    if (!iterator) return false;
    std::unique_ptr<bam1_t, AlignmentDeleter> alignment(bam_init1());
    if (!alignment) return false;
    std::unordered_set<std::string> seen;
    std::array<int, 2> relation_support{};
    double log_odds = 0.0;
    while (sam_itr_next(context.bams.front()->get(), iterator.get(),
                        alignment.get()) >= 0) {
        const bam1_t* read = alignment.get();
        if ((read->core.flag & (BAM_FUNMAP | BAM_FSECONDARY |
                                BAM_FSUPPLEMENTARY | BAM_FDUP | BAM_FQCFAIL)) ||
            read->core.qual < kMinMapq ||
            read->core.qual == kUnknownQuality ||
            read->core.pos > left_pos - 1 ||
            bam_endpos(read) <= nearest ||
            seen.count(bam_get_qname(read)) != 0)
            continue;
        const auto graph = graph_haps.find(bam_get_qname(read));
        if (graph == graph_haps.end() || graph->second.first == 0) continue;
        const bool a_matches = bam_matches_deletion_sequence(
            read, a, context.ref, tid, context.primary_header(), kMinDeletionBaseq);
        const bool b_matches = bam_matches_deletion_sequence(
            read, b, context.ref, tid, context.primary_header(), kMinDeletionBaseq);
        // Failure to certify one ALT is not a REF observation of that row.
        if (a_matches == b_matches) continue;
        const CandidateVariant& deletion = a_matches ? a : b;
        const int right_hap = deletion.hap_to_cons_alle[1] == 1 ? 1 : 2;
        const bool flip = graph->second.first != right_hap;
        seen.insert(bam_get_qname(read));
        ++relation_support[flip ? 1 : 0];
        const double p = std::pow(10.0, -kMinDeletionBaseq / 10.0) +
            std::pow(10.0, -read->core.qual / 10.0) +
            std::pow(10.0, -graph->second.second / 10.0);
        if (p > 0.0 && p < 0.5)
            log_odds += (flip ? 1.0 : -1.0) *
                std::log((1.0 - p) / p);
    }
    const double threshold = std::log((1.0 - kMaxWrongParity) /
                                      kMaxWrongParity);
    if (relation_support[0] + relation_support[1] < kMinPairedReads ||
        (relation_support[0] > 0 && relation_support[1] > 0) ||
        std::abs(log_odds) < threshold ||
        (required_flip && *required_flip != (log_odds > 0.0)))
        return false;
    return merge_phase_sets_in_place(chunk, left_ps, right_ps,
                                     log_odds > 0.0);
}

// CIGAR can shift an insertion within a repeat. Compare the edited
// reference strings and require clean bases across both placements.
static constexpr size_t kLongPhysicalInsertionThreshold = 64;
static int physical_equivalent_insertion_call(
        const bam1_t* read, const CandidateVariant& insertion,
        WorkerContext& context, int tid, int min_baseq,
        int* allele_quality = nullptr, int min_insertion_baseq = -1) {
    constexpr size_t kMaxInsertionLength = 128;
    constexpr int kUnknownQuality = 255;
    const hts_pos_t target_pos = insertion.key.pos;
    const std::string& target_bases = insertion.key.alt;
    constexpr size_t kShortInsertionMaxLength = 2;
    constexpr hts_pos_t kMaxShortEquivalentShift = 64;
    constexpr hts_pos_t kMaxLongEquivalentShift = 16;
    const hts_pos_t kMaxEquivalentShift =
        target_bases.size() <= kShortInsertionMaxLength ?
            kMaxShortEquivalentShift : kMaxLongEquivalentShift;
    if (target_pos <= kMaxEquivalentShift || target_bases.empty() ||
        target_bases.size() > kMaxInsertionLength) return -1;
    const hts_pos_t window_beg = target_pos - kMaxEquivalentShift;
    const hts_pos_t window_end = target_pos + kMaxEquivalentShift;
    hts_pos_t ref_pos = read->core.pos + 1;
    int query_pos = 0;
    hts_pos_t observed_pos = -1;
    std::string observed_bases;
    int observed_quality = kUnknownQuality;
    const int insertion_baseq = min_insertion_baseq < 0 ?
        min_baseq : min_insertion_baseq;
    int nearby_indels = 0;
    const hts_pos_t interference_flank =
        target_bases.size() > kLongPhysicalInsertionThreshold ?
            static_cast<hts_pos_t>(target_bases.size()) :
            kMaxEquivalentShift;
    const uint32_t* cigar = bam_get_cigar(read);
    const uint8_t* sequence = bam_get_seq(read);
    const uint8_t* qualities = bam_get_qual(read);
    for (uint32_t ci = 0; ci < read->core.n_cigar; ++ci) {
        const int op = bam_cigar_op(cigar[ci]);
        const int length = bam_cigar_oplen(cigar[ci]);
        // A nearby different indel can make a long insertion's apparent REF
        // call an alternate placement of the same complex repeat allele.
        if (target_bases.size() > kLongPhysicalInsertionThreshold &&
            ((op == BAM_CINS &&
              ref_pos >= target_pos - interference_flank &&
              ref_pos <= target_pos + interference_flank &&
              (ref_pos < window_beg || ref_pos > window_end)) ||
             (op == BAM_CDEL &&
              ref_pos < target_pos + interference_flank &&
              ref_pos + length > target_pos - interference_flank &&
              (ref_pos + length <= window_beg || ref_pos >= window_end))))
            return -1;
        if (op == BAM_CINS && ref_pos >= window_beg &&
            ref_pos <= window_end) {
            ++nearby_indels;
            observed_pos = ref_pos;
            for (int qi = query_pos; qi < query_pos + length; ++qi) {
                if (qualities[qi] < insertion_baseq ||
                    qualities[qi] == kUnknownQuality) return -1;
                observed_quality = std::min(observed_quality,
                                            static_cast<int>(qualities[qi]));
                observed_bases.push_back(static_cast<char>(std::toupper(
                    seq_nt16_str[bam_seqi(sequence, qi)])));
            }
        } else if (op == BAM_CDEL && ref_pos < window_end &&
                   ref_pos + length > window_beg) {
            ++nearby_indels;
        }
        if (bam_cigar_type(op) & 1) query_pos += length;
        if (bam_cigar_type(op) & 2) ref_pos += length;
    }
    if (nearby_indels > 1 || (nearby_indels == 1 && observed_pos < 0))
        return -1;
    const hts_pos_t beg = nearby_indels == 1 ?
        std::min(target_pos, observed_pos) : target_pos;
    const hts_pos_t end = nearby_indels == 1 ?
        std::max(target_pos, observed_pos) : target_pos;
    if (nearby_indels == 1) {
        if (observed_bases.size() != target_bases.size() ||
            observed_pos < window_beg || observed_pos > window_end)
            return -1;
        std::string ref;
        for (hts_pos_t pos = beg; pos < end; ++pos) {
            const char base = context.ref.base(
                tid, pos, context.primary_header());
            if (base == 'N') return -1;
            ref.push_back(static_cast<char>(std::toupper(base)));
        }
        std::string expected = ref;
        expected.insert(static_cast<size_t>(target_pos - beg), target_bases);
        ref.insert(static_cast<size_t>(observed_pos - beg), observed_bases);
        if (ref != expected) return -1;
    }
    for (hts_pos_t pos = beg - 1; pos <= end; ++pos) {
        const char ref = context.ref.base(
            tid, pos, context.primary_header());
        int quality = 0;
        if (ref == 'N' ||
            physical_snp_call(read, pos, ref, 'N', &quality) != 0 ||
            quality < min_baseq || quality == kUnknownQuality)
            return -1;
    }
    // Preserve the existing nearby-indel interference check. Only a verified
    // shifted ALT can replace an apparent REF outside that bounded search;
    // unrelated edits elsewhere in a long repeat do not weaken a callable ALT.
    if (nearby_indels == 0 && target_bases.size() > kShortInsertionMaxLength) {
        const int query_index = bam_shifted_repeat_insertion_query_index(
            read, insertion, context.ref, tid, context.primary_header(), min_baseq);
        if (query_index >= 0) {
            if (allele_quality != nullptr) {
                int quality = kUnknownQuality;
                for (size_t offset = 0; offset < target_bases.size(); ++offset)
                    quality = std::min(quality, static_cast<int>(
                        qualities[query_index + offset]));
                *allele_quality = quality;
            }
            return 1;
        }
    }
    if (allele_quality != nullptr)
        *allele_quality = nearby_indels == 1 ? observed_quality : min_baseq;
    return nearby_indels == 1 ? 1 : 0;
}

// A recovered deletion at the end of one block can be oriented by reads
// already assigned to the next block. At a complementary deletion locus,
// calling either row REF is ambiguous, so use the right block's observed HP
// gauge and move only the left boundary row across the weak source edge.
static bool stitch_deletion_to_phased_right_reads(
        GraphChunkBuildResult& gc, WorkerContext& context, int tid,
        const RecoverySeam& seam) {
    constexpr int kMinMapq = 30;
    constexpr int kMinBaseq = 30;
    constexpr int kUnknownQuality = 255;
    constexpr int kMinSupportPerAllele = 2;
    constexpr hts_pos_t kMaxGaugeSiteDistance = 5000;
    PhasingChunk& chunk = gc.chunk;
    std::optional<size_t> left_i;
    std::vector<size_t> right_rows;
    hts_pos_t right_source = 0;
    for (const RecoverySourceSite& source : gc.recovery_source_sites) {
        const size_t ci = source.candidate_index;
        if (ci >= chunk.candidates.size()) continue;
        const CandidateVariant& candidate = chunk.candidates[ci];
        if (!candidate.bam_injected || !candidate.msa_verified ||
            candidate.key.type != VariantType::Deletion ||
            (candidate.counts.category != VariantCategory::CleanHetIndel &&
             candidate.counts.category != VariantCategory::NoisyCandHet) ||
            candidate.hap_to_cons_alle[1] < 0 ||
            candidate.hap_to_cons_alle[1] > 1 ||
            candidate.hap_to_cons_alle[2] !=
                1 - candidate.hap_to_cons_alle[1])
            continue;
        if (candidate.phase_set == seam.left_phase_set &&
            candidate.key.sort_pos() == seam.beg) {
            if (left_i) return false;
            left_i = ci;
        }
        if (candidate.phase_set == seam.right_phase_set &&
            candidate.key.sort_pos() == seam.end) {
            if (right_source != 0 && right_source != source.phase_set)
                return false;
            right_source = source.phase_set;
            right_rows.push_back(ci);
        }
    }
    if (!left_i || right_rows.size() != 2 ||
        chunk.candidates[right_rows[0]].key.ref_len ==
            chunk.candidates[right_rows[1]].key.ref_len ||
        chunk.candidates[right_rows[0]].hap_to_cons_alle[1] ==
            chunk.candidates[right_rows[1]].hap_to_cons_alle[1])
        return false;
    const CandidateVariant& left = chunk.candidates[*left_i];
    // A boundary row may leave its old block only when no other oriented row
    // follows it. This preserves the earlier block and all its read labels.
    for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
        if (ci == *left_i) continue;
        const CandidateVariant& candidate = chunk.candidates[ci];
        if (candidate.phase_set == seam.left_phase_set &&
            candidate.key.sort_pos() >= left.key.sort_pos() &&
            candidate.hap_to_cons_alle[1] >= 0 &&
            candidate.hap_to_cons_alle[1] !=
                candidate.hap_to_cons_alle[2])
            return false;
    }
    hts_pos_t right_end = seam.end;
    for (const size_t ci : right_rows)
        right_end = std::max(right_end,
            chunk.candidates[ci].key.pos +
            chunk.candidates[ci].key.ref_len);
    std::unordered_map<std::string, size_t> right_reads;
    for (size_t ri = 0; ri < chunk.reads.size(); ++ri)
        if (ri < chunk.phase_sets.size() && ri < chunk.haps.size() &&
            chunk.phase_sets[ri] == seam.right_phase_set &&
            (chunk.haps[ri] == 1 || chunk.haps[ri] == 2))
            right_reads.emplace(chunk.reads[ri].qname, ri);
    if (right_reads.empty()) return false;
    std::unique_ptr<hts_itr_t, decltype(&hts_itr_destroy)> iterator(
        sam_itr_queryi(context.indexes.front().get(), tid,
                       left.key.pos - 1, right_end), &hts_itr_destroy);
    if (!iterator) return false;
    std::unique_ptr<bam1_t, AlignmentDeleter> alignment(bam_init1());
    if (!alignment) return false;
    std::unordered_set<std::string> seen;
    std::array<int, 2> parity_votes{};
    std::array<int, 2> allele_support{};
    while (sam_itr_next(context.bams.front()->get(), iterator.get(),
                        alignment.get()) >= 0) {
        const bam1_t* read = alignment.get();
        if ((read->core.flag & (BAM_FUNMAP | BAM_FSECONDARY |
                                BAM_FSUPPLEMENTARY | BAM_FDUP | BAM_FQCFAIL)) ||
            read->core.qual < kMinMapq ||
            read->core.qual == kUnknownQuality ||
            bam_endpos(read) < right_end ||
            seen.count(bam_get_qname(read)) != 0)
            continue;
        const auto assigned = right_reads.find(bam_get_qname(read));
        if (assigned == right_reads.end()) continue;
        const size_t ri = assigned->second;
        if (ri >= chunk.read_var_profile.size()) continue;
        const ReadVariantProfile& profile = chunk.read_var_profile[ri];
        if (profile.start_var_idx < 0) continue;
        bool observes_right = false;
        for (size_t offset = 0; offset < profile.alleles.size(); ++offset) {
            const size_t ci =
                static_cast<size_t>(profile.start_var_idx) + offset;
            if (ci >= chunk.candidates.size()) break;
            const CandidateVariant& site = chunk.candidates[ci];
            const hts_pos_t pos = site.key.sort_pos();
            if (site.phase_set != seam.right_phase_set ||
                pos < seam.end ||
                pos > seam.end + kMaxGaugeSiteDistance ||
                site.hap_to_cons_alle[1] < 0 ||
                site.hap_to_cons_alle[1] ==
                    site.hap_to_cons_alle[2])
                continue;
            if (profile.alleles[offset] >= 0 ||
                (offset < profile.bam_alleles.size() &&
                 profile.bam_alleles[offset] >= 0) ||
                (offset < profile.graph_alleles.size() &&
                 profile.graph_alleles[offset] >= 0)) {
                observes_right = true;
                break;
            }
        }
        if (!observes_right) continue;
        const int left_call = physical_equivalent_deletion_call(
            read, left, context, tid, kMinBaseq);
        if (left_call < 0) continue;
        seen.insert(bam_get_qname(read));
        ++allele_support[static_cast<size_t>(left_call)];
        const bool left_hap1 = (left_call == 1) ==
            (left.hap_to_cons_alle[1] == 1);
        const bool right_hap1 = chunk.haps[ri] == 1;
        ++parity_votes[left_hap1 != right_hap1 ? 1 : 0];
    }
    if (allele_support[0] < kMinSupportPerAllele ||
        allele_support[1] < kMinSupportPerAllele ||
        (parity_votes[0] > 0 && parity_votes[1] > 0) ||
        !source_graph_vote_supported(
            parity_votes[0], parity_votes[1],
            parity_votes[1] > parity_votes[0] ? 1 : 0))
        return false;
    CandidateVariant& moved = chunk.candidates[*left_i];
    if (parity_votes[1] > 0) {
        std::swap(moved.hap_to_cons_alle[1],
                  moved.hap_to_cons_alle[2]);
        std::swap(moved.hap_to_alle_profile[1],
                  moved.hap_to_alle_profile[2]);
        const auto flip_hap = [](int hap) {
            return hap == 1 ? 2 : hap == 2 ? 1 : hap;
        };
        moved.hap_alt = flip_hap(moved.hap_alt);
        moved.hap_ref = flip_hap(moved.hap_ref);
    }
    moved.phase_set = seam.right_phase_set;
    return true;
}

// Two overlapping BAM deletion rows can encode opposite haplotypes. Their
// ALT calls certify a local source run even when the preceding graph SNP has
// no read edge into it; ambiguous REF calls at that locus never cast a vote.
static bool stitch_complementary_left_deletions(
        GraphChunkBuildResult& gc, WorkerContext& context, int tid,
        const RecoverySeam& seam,
        const std::map<std::string, hts_pos_t>* original_graph_phase_sets,
        const std::map<std::string, hts_pos_t>* stitched_graph_phase_sets) {
    constexpr int kMinAlleleClassSupport = 2;
    constexpr double kMaxWrongParity = 0.001;
    constexpr int kMinMapq = 30;
    constexpr int kUnknownMapq = 255;
    PhasingChunk& chunk = gc.chunk;
    const RecoverySourceSite* boundary = nullptr;
    const RecoverySourceSite* right = nullptr;
    for (const RecoverySourceSite& source : gc.recovery_source_sites) {
        if (source.candidate_index >= chunk.candidates.size()) continue;
        const CandidateVariant& site = chunk.candidates[source.candidate_index];
        if (!site.bam_injected || !site.msa_verified ||
            site.key.type != VariantType::Deletion ||
            site.hap_to_cons_alle[1] < 0 ||
            site.hap_to_cons_alle[1] > 1 ||
            site.hap_to_cons_alle[2] != 1 - site.hap_to_cons_alle[1])
            continue;
        if (site.phase_set == seam.left_phase_set &&
            site.key.sort_pos() == seam.beg) {
            if (boundary) return false;
            boundary = &source;
        }
        if (site.phase_set == seam.right_phase_set &&
            site.key.sort_pos() == seam.end) {
            if (right) return false;
            right = &source;
        }
    }
    if (!boundary || !right) return false;
    const RecoverySourceSite* first = nullptr;
    const RecoverySourceSite* second = nullptr;
    hts_pos_t pair_pos = 0;
    for (size_t i = 0; i < gc.recovery_source_sites.size(); ++i) {
        const RecoverySourceSite& a = gc.recovery_source_sites[i];
        if (a.phase_set != boundary->phase_set ||
            a.candidate_index >= chunk.candidates.size()) continue;
        const CandidateVariant& left_a = chunk.candidates[a.candidate_index];
        if (left_a.phase_set != seam.left_phase_set ||
            !left_a.bam_injected || !left_a.msa_verified ||
            left_a.key.type != VariantType::Deletion ||
            left_a.key.sort_pos() >= seam.beg ||
            left_a.key.sort_pos() <= pair_pos)
            continue;
        for (size_t j = i + 1; j < gc.recovery_source_sites.size(); ++j) {
            const RecoverySourceSite& b = gc.recovery_source_sites[j];
            if (b.phase_set != a.phase_set ||
                b.candidate_index >= chunk.candidates.size()) continue;
            const CandidateVariant& left_b = chunk.candidates[b.candidate_index];
            if (left_b.phase_set != seam.left_phase_set ||
                !left_b.bam_injected || !left_b.msa_verified ||
                left_b.key.type != VariantType::Deletion ||
                left_a.key.sort_pos() != left_b.key.sort_pos() ||
                left_a.key.ref_len == left_b.key.ref_len ||
                left_a.hap_to_cons_alle[1] == left_b.hap_to_cons_alle[1] ||
                std::max(left_a.key.pos, left_b.key.pos) >=
                    std::min(left_a.key.pos + left_a.key.ref_len,
                             left_b.key.pos + left_b.key.ref_len))
                continue;
            first = &a;
            second = &b;
            pair_pos = left_a.key.sort_pos();
        }
    }
    if (!first || !second) return false;
    size_t rows_at_pair = 0;
    for (const RecoverySourceSite& source : gc.recovery_source_sites) {
        if (source.phase_set != boundary->phase_set ||
            source.candidate_index >= chunk.candidates.size()) continue;
        const CandidateVariant& site = chunk.candidates[source.candidate_index];
        if (site.phase_set == seam.left_phase_set &&
            site.bam_injected && site.msa_verified &&
            site.key.type == VariantType::Deletion &&
            site.key.sort_pos() == pair_pos)
            ++rows_at_pair;
    }
    if (rows_at_pair != 2) return false;
    const auto has_internal_cut = [&](const auto& cuts) {
        const auto found = cuts.find(boundary->phase_set);
        return found != cuts.end() &&
            std::any_of(found->second.begin(), found->second.end(),
                        [pair_pos, &seam](hts_pos_t pos) {
                            return pos > pair_pos && pos <= seam.beg;
                        });
    };
    if (has_internal_cut(gc.recovery_source_weak_cuts) ||
        has_internal_cut(gc.recovery_source_quality_cuts))
        return false;
    hts_pos_t preceding_pos = 0;
    bool has_graph_prefix = false;
    for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
        const CandidateVariant& site = chunk.candidates[ci];
        if (site.phase_set != seam.left_phase_set ||
            site.hap_to_cons_alle[1] < 0 ||
            site.hap_to_cons_alle[1] == site.hap_to_cons_alle[2])
            continue;
        VariantKey physical = site.key;
        if (!site.bam_injected) {
            if (ci >= gc.site_meta.size()) return false;
            const std::string* alt = selected_graph_candidate_alt(gc, ci);
            if (alt == nullptr) return false;
            const GraphSiteMeta& meta = gc.site_meta[ci];
            physical = vcf_to_variant_key(
                site.key.tid, meta.pos, meta.ref, *alt);
            if (physical.sort_pos() >= pair_pos) return false;
            has_graph_prefix = true;
        }
        if (physical.sort_pos() < pair_pos)
            preceding_pos = std::max(preceding_pos, physical.sort_pos());
    }
    if (!has_graph_prefix || preceding_pos == 0) return false;
    const CandidateVariant& a = chunk.candidates[first->candidate_index];
    const CandidateVariant& b = chunk.candidates[second->candidate_index];
    const CandidateVariant& r = chunk.candidates[right->candidate_index];
    std::unordered_set<std::string> seen;
    std::unordered_map<std::string, int> bridge_haps;
    std::array<int, 2> allele_support{};
    std::array<int, 2> parity_votes{};
    for (size_t ri = 0; ri < chunk.read_var_profile.size() &&
                        ri < chunk.reads.size(); ++ri) {
        const ReadRecord& read = chunk.reads[ri];
        const ReadVariantProfile& profile = chunk.read_var_profile[ri];
        // These are BAM-channel allele calls. A low GAF MAPQ must not
        // discard an independently well-mapped BAM alignment.
        if (read.is_skipped || profile.bam_mapq < kMinMapq ||
            profile.bam_mapq == kUnknownMapq ||
            seen.count(read.qname) != 0)
            continue;
        if (profile.start_var_idx < 0) continue;
        const auto call = [&](size_t ci) {
            if (ci < static_cast<size_t>(profile.start_var_idx)) return -1;
            const size_t offset = ci - static_cast<size_t>(profile.start_var_idx);
            return offset < profile.bam_alleles.size() ?
                static_cast<int>(profile.bam_alleles[offset]) : -1;
        };
        const int a_call = call(first->candidate_index);
        const int b_call = call(second->candidate_index);
        const int r_call = call(right->candidate_index);
        if (r_call < 0 || r_call > 1 ||
            !((a_call == 1 && b_call == 0) ||
              (a_call == 0 && b_call == 1)))
            continue;
        seen.insert(read.qname);
        const CandidateVariant& left = a_call == 1 ? a : b;
        const bool left_hap1 = left.hap_to_cons_alle[1] == 1;
        const bool right_hap1 = r_call == r.hap_to_cons_alle[1];
        bridge_haps.emplace(read.qname, right_hap1 ? 1 : 2);
        ++allele_support[static_cast<size_t>(left_hap1)];
        ++parity_votes[left_hap1 != right_hap1 ? 1 : 0];
    }
    const int paired = parity_votes[0] + parity_votes[1];
    // With unanimous calls the one-sided fair-parity tail is exactly 2^-n.
    if (allele_support[0] < kMinAlleleClassSupport ||
        allele_support[1] < kMinAlleleClassSupport ||
        (parity_votes[0] > 0 && parity_votes[1] > 0) ||
        std::ldexp(1.0, -paired) > kMaxWrongParity)
        return false;
    // The paired deletions may sit inside a longer supported BAM run. Its
    // nearest clean SNP can certify the earlier part without absorbing the
    // preceding catalog block. An ALT-absence call alone is ambiguous at a
    // two-deletion locus; only mutually exclusive ALT observations can extend it.
    hts_pos_t transfer_beg = pair_pos;
    hts_pos_t transfer_prefix = preceding_pos;
    hts_pos_t graph_end = 0;
    hts_pos_t graph_pos = 0;
    for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
        const CandidateVariant& candidate = chunk.candidates[ci];
        if (candidate.phase_set != seam.left_phase_set ||
            !is_phase_set_anchor(candidate)) continue;
        if (!candidate.bam_injected) {
            if (ci >= gc.site_meta.size()) continue;
            const std::string* alt = selected_graph_candidate_alt(gc, ci);
            if (alt == nullptr) continue;
            const GraphSiteMeta& meta = gc.site_meta[ci];
            const VariantKey physical = vcf_to_variant_key(
                candidate.key.tid, meta.pos, meta.ref, *alt);
            if (physical.sort_pos() >= pair_pos) continue;
            graph_pos = std::max(graph_pos, physical.sort_pos());
            graph_end = std::max(graph_end, physical.pos +
                std::max(1, physical.ref_len) - 1);
        }
    }
    const auto prefix = bam_prefix_before_deletion_pair(
        gc, first->candidate_index, second->candidate_index, graph_end);
    if (prefix && graph_pos > 0 && graph_pos < *prefix) {
        transfer_beg = *prefix;
        transfer_prefix = graph_pos;
    }
    if (!join_graph_suffix_after_weak_edge(
            gc, context, tid, seam.left_phase_set, seam.right_phase_set,
            parity_votes[1] > 0, boundary->candidate_index,
            {transfer_prefix, transfer_beg}, original_graph_phase_sets,
            stitched_graph_phase_sets, true))
        return false;
    // The generic split clears reads calling both sides of its cut. A read
    // that also calls the certified deletion pair and right deletion already
    // has an independent right-block allele gauge, so keep that assignment.
    for (size_t ri = 0; ri < chunk.reads.size() &&
                        ri < chunk.phase_sets.size() &&
                        ri < chunk.haps.size(); ++ri) {
        if (chunk.phase_sets[ri] != 0 || chunk.haps[ri] != 0) continue;
        const auto found = bridge_haps.find(chunk.reads[ri].qname);
        if (found == bridge_haps.end()) continue;
        chunk.phase_sets[ri] = seam.right_phase_set;
        chunk.haps[ri] = found->second;
    }
    return true;
}

// A BAM source can contain a well-linked insertion run after a weak cut while
// the neighboring graph block extends past a second unsupported edge. Orient
// just that run from reads already phased on the left; never absorb the graph
// block or let the inherited BAM PS stand in for an observed allele edge.
static bool stitch_bam_insertion_run_to_left_reads(
        GraphChunkBuildResult& gc, const RecoverySeam& seam,
        const std::map<std::string, hts_pos_t>* original_graph_phase_sets,
        const std::map<std::string, hts_pos_t>* stitched_graph_phase_sets) {
    constexpr int kMinMapq = 30;
    constexpr int kUnknownMapq = 255;
    constexpr int kMinAlleleClassSupport = 2;
    constexpr double kMaxInternalEdgeP = 0.01;
    constexpr double kMaxReadGaugeP = 0.05;
    PhasingChunk& chunk = gc.chunk;
    std::optional<size_t> insertion_i;
    hts_pos_t source_ps = 0;
    for (const RecoverySourceSite& source : gc.recovery_source_sites) {
        if (source.candidate_index >= chunk.candidates.size()) continue;
        const CandidateVariant& candidate =
            chunk.candidates[source.candidate_index];
        if (candidate.phase_set != seam.right_phase_set ||
            candidate.key.sort_pos() != seam.end ||
            candidate.key.type != VariantType::Insertion ||
            !candidate.bam_injected || !candidate.msa_verified ||
            !candidate.alignment_verified ||
            candidate.hap_to_cons_alle[1] < 0 ||
            candidate.hap_to_cons_alle[1] > 1 ||
            candidate.hap_to_cons_alle[2] !=
                1 - candidate.hap_to_cons_alle[1])
            continue;
        if (insertion_i) return false;
        insertion_i = source.candidate_index;
        source_ps = source.phase_set;
    }
    if (!insertion_i || source_ps <= 0) return false;
    const auto cuts = gc.recovery_source_weak_cuts.find(source_ps);
    if (cuts == gc.recovery_source_weak_cuts.end() ||
        !std::binary_search(cuts->second.begin(), cuts->second.end(), seam.beg))
        return false;
    hts_pos_t first_graph_pos = std::numeric_limits<hts_pos_t>::max();
    for (const CandidateVariant& candidate : chunk.candidates)
        if (candidate.phase_set == seam.right_phase_set &&
            !candidate.bam_injected)
            first_graph_pos = std::min(first_graph_pos,
                                       candidate.key.sort_pos());
    if (first_graph_pos <= seam.end ||
        first_graph_pos == std::numeric_limits<hts_pos_t>::max())
        return false;
    bool left_source_anchor = false;
    std::vector<size_t> following;
    for (const RecoverySourceSite& source : gc.recovery_source_sites) {
        if (source.phase_set != source_ps ||
            source.candidate_index >= chunk.candidates.size()) continue;
        const CandidateVariant& candidate =
            chunk.candidates[source.candidate_index];
        const hts_pos_t pos = candidate.key.sort_pos();
        left_source_anchor |= pos == seam.beg &&
            candidate.phase_set == seam.left_phase_set &&
            candidate.hap_to_cons_alle[1] >= 0 &&
            candidate.hap_to_cons_alle[1] <= 1 &&
            candidate.hap_to_cons_alle[2] ==
                1 - candidate.hap_to_cons_alle[1];
        if (pos <= seam.end || pos >= first_graph_pos) continue;
        if (candidate.phase_set != seam.right_phase_set ||
            !candidate.bam_injected ||
            candidate.hap_to_cons_alle[1] < 0 ||
            candidate.hap_to_cons_alle[2] !=
                1 - candidate.hap_to_cons_alle[1])
            return false;
        following.push_back(source.candidate_index);
    }
    if (!left_source_anchor || following.empty() ||
        !graph_snp_path_supported(gc, seam.left_phase_set, nullptr,
                                  original_graph_phase_sets,
                                  stitched_graph_phase_sets, false, true))
        return false;
    std::sort(following.begin(), following.end(),
              [&chunk](size_t a, size_t b) {
                  const hts_pos_t a_pos = chunk.candidates[a].key.sort_pos();
                  const hts_pos_t b_pos = chunk.candidates[b].key.sort_pos();
                  return a_pos == b_pos ? a < b : a_pos < b_pos;
              });
    std::vector<size_t> run{*insertion_i};
    for (const size_t ci : following) {
        // This transfers observed BAM rows, without calling an allele from
        // its length or absorbing the next graph site. A long insertion needs
        // the same paired-allele certificate as every other edge in the run.
        const std::optional<bool> internal_flip = local_run_boundary_flip(
            chunk, run.back(), ci, kMinMapq, kMaxInternalEdgeP);
        if (!internal_flip || *internal_flip) break;
        run.push_back(ci);
    }
    if (run.size() < 2) return false;
    const CandidateVariant& insertion = chunk.candidates[*insertion_i];
    std::array<int, 2> allele_support{};
    std::array<int, 2> parity_votes{};
    std::unordered_set<std::string> seen;
    for (size_t ri = 0; ri < chunk.reads.size(); ++ri) {
        if (ri >= chunk.phase_sets.size() || ri >= chunk.haps.size() ||
            ri >= chunk.read_var_profile.size() ||
            chunk.phase_sets[ri] != seam.left_phase_set ||
            (chunk.haps[ri] != 1 && chunk.haps[ri] != 2) ||
            chunk.reads[ri].is_skipped ||
            chunk.reads[ri].mapq < kMinMapq ||
            chunk.reads[ri].mapq == kUnknownMapq ||
            seen.count(chunk.reads[ri].qname) != 0)
            continue;
        const ReadVariantProfile& profile = chunk.read_var_profile[ri];
        if (profile.start_var_idx < 0 ||
            *insertion_i < static_cast<size_t>(profile.start_var_idx))
            continue;
        const size_t offset = *insertion_i -
            static_cast<size_t>(profile.start_var_idx);
        if (offset >= profile.alleles.size()) continue;
        const int allele = profile.alleles[offset];
        if (allele != 0 && allele != 1) continue;
        seen.insert(chunk.reads[ri].qname);
        ++allele_support[allele];
        const bool source_hap1 = allele == insertion.hap_to_cons_alle[1];
        const bool left_hap1 = chunk.haps[ri] == 1;
        ++parity_votes[source_hap1 != left_hap1 ? 1 : 0];
    }
    const bool flip = parity_votes[1] > parity_votes[0];
    if (allele_support[0] < kMinAlleleClassSupport ||
        allele_support[1] < kMinAlleleClassSupport ||
        !source_graph_vote_supported(parity_votes[0], parity_votes[1],
                                     flip ? 1 : 0, kMaxReadGaugeP))
        return false;
    const auto flip_hap = [](int hap) {
        return hap == 1 ? 2 : hap == 2 ? 1 : hap;
    };
    for (const size_t ci : run) {
        CandidateVariant& candidate = chunk.candidates[ci];
        if (flip) {
            std::swap(candidate.hap_to_cons_alle[1],
                      candidate.hap_to_cons_alle[2]);
            std::swap(candidate.hap_to_alle_profile[1],
                      candidate.hap_to_alle_profile[2]);
            candidate.hap_alt = flip_hap(candidate.hap_alt);
            candidate.hap_ref = flip_hap(candidate.hap_ref);
        }
        candidate.phase_set = seam.left_phase_set;
    }
    return true;
}

// Compare both insertion edits against the BAM source's common background.
static bool insertion_edits_match(const VariantKey& left, const VariantKey& right,
                                  const PhasingChunk& chunk,
                                  WorkerContext& context, int tid) {
    if (left.type != VariantType::Insertion || right.type != VariantType::Insertion ||
        left.tid != tid || right.tid != tid || left.ref_len != 0 || right.ref_len != 0 ||
        left.pos > right.pos || left.alt.empty() || left.alt.size() != right.alt.size() ||
        right.pos - left.pos > std::numeric_limits<int>::max())
        return false;
    const size_t shift = static_cast<size_t>(right.pos - left.pos);
    std::string background = shift == 0 ? std::string() :
        context.ref.subseq(tid, left.pos, static_cast<int>(shift), context.primary_header());
    if (background.size() != shift) return false;
    // A homozygous source SNP belongs to both haplotypes. Omitting it can
    // falsely distinguish a repeat-shifted edit from its graph representation.
    // Conflicting common alleles cannot certify equivalence.
    std::map<hts_pos_t, char> common;
    for (const CandidateVariant& candidate : chunk.candidates) {
        if (!candidate.bam_injected || candidate.key.type != VariantType::Snp ||
            candidate.key.ref_len != 1 || candidate.key.alt.size() != 1 ||
            (candidate.counts.category != VariantCategory::CleanHom &&
             !(candidate.counts.category == VariantCategory::NoisyCandHom && candidate.msa_verified)) ||
            candidate.hap_to_cons_alle[1] != 1 || candidate.hap_to_cons_alle[2] != 1 ||
            candidate.key.pos < left.pos || candidate.key.pos >= right.pos) continue;
        const char alt = candidate.key.alt[0];
        const auto [it, inserted] = common.emplace(candidate.key.pos, alt);
        if (!inserted && base_to_nt4(it->second) != base_to_nt4(alt)) return false;
        background[static_cast<size_t>(candidate.key.pos - left.pos)] = alt;
    }
    return insertion_edits_are_equivalent(left.pos, left.alt, right.pos,
                                          right.alt, background);
}

// Equivalent edits certify allele identity; independent paired calls certify
// the two phase gauges. Keep both candidate rows and all original observations.
static std::optional<std::pair<size_t, size_t>> find_equivalent_source_insertion_join(
        const GraphChunkBuildResult& gc, WorkerContext& context, int tid,
        const RecoverySeam& seam,
        const std::map<std::string, hts_pos_t>* original_graph_phase_sets,
        const std::map<std::string, hts_pos_t>* stitched_graph_phase_sets) {
    constexpr int kMinMapq = 30;
    constexpr int kUnknownMapq = 255;
    constexpr int kMinPerAllele = 2;
    constexpr double kMaxParityP = 0.001;
    const PhasingChunk& chunk = gc.chunk;
    std::optional<size_t> left_i, right_i;
    VariantKey right_key;
    for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
        const CandidateVariant& candidate = chunk.candidates[ci];
        if (!is_phase_set_anchor(candidate) || candidate.key.type != VariantType::Insertion)
            continue;
        if (candidate.phase_set == seam.left_phase_set &&
            candidate.key.sort_pos() == seam.beg && candidate.bam_injected &&
            candidate.msa_verified && candidate.alignment_verified &&
            candidate.msa_insertion_alts.empty() &&
            candidate.hap_to_cons_alle[1] >= 0 && candidate.hap_to_cons_alle[1] <= 1 &&
            candidate.hap_to_cons_alle[2] == 1 - candidate.hap_to_cons_alle[1]) {
            if (left_i) return std::nullopt;
            left_i = ci;
        }
        if (candidate.phase_set != seam.right_phase_set || candidate.bam_injected ||
            ci >= gc.site_meta.size() || candidate.hap_to_cons_alle[1] < 0 ||
            candidate.hap_to_cons_alle[1] > 1 ||
            candidate.hap_to_cons_alle[2] != 1 - candidate.hap_to_cons_alle[1]) continue;
        const std::string* alt = selected_graph_candidate_alt(gc, ci);
        if (alt == nullptr) continue;
        const GraphSiteMeta& meta = gc.site_meta[ci];
        if (ci >= gc.site_allele_orig_idx.size() || gc.site_allele_orig_idx[ci].size() != 2 ||
            gc.site_allele_orig_idx[ci][0] != 0) continue;
        const VariantKey physical = vcf_to_variant_key(tid, meta.pos, meta.ref, *alt);
        if (physical.type != VariantType::Insertion || physical.sort_pos() != seam.end) continue;
        if (right_i) return std::nullopt;
        right_i = ci;
        right_key = physical;
    }
    if (!left_i || !right_i) return std::nullopt;
    if (!insertion_edits_match(chunk.candidates[*left_i].key, right_key, chunk, context, tid)) return std::nullopt;
    std::array<int, 2> counts{};
    std::unordered_set<std::string> seen;
    for (size_t ri = 0; ri < chunk.read_var_profile.size() && ri < chunk.reads.size(); ++ri) {
        const ReadRecord& read = chunk.reads[ri];
        if (read.is_skipped || read.mapq < kMinMapq || read.mapq == kUnknownMapq ||
            seen.count(read.qname) != 0) continue;
        const ReadVariantProfile& profile = chunk.read_var_profile[ri];
        if (profile.bam_mapq < kMinMapq || profile.bam_mapq == kUnknownMapq) continue;
        if (profile.start_var_idx < 0 || *left_i < static_cast<size_t>(profile.start_var_idx) ||
            *right_i < static_cast<size_t>(profile.start_var_idx)) continue;
        const size_t lo = *left_i - static_cast<size_t>(profile.start_var_idx);
        const size_t ro = *right_i - static_cast<size_t>(profile.start_var_idx);
        if (lo >= profile.bam_alleles.size() || ro >= profile.graph_alleles.size()) continue;
        const int a = profile.bam_alleles[lo], b = profile.graph_alleles[ro];
        if ((a != 0 && a != 1) || (b != 0 && b != 1)) continue;
        if (a != b) return std::nullopt;
        seen.insert(read.qname);
        ++counts[a];
    }
    if (counts[0] < kMinPerAllele || counts[1] < kMinPerAllele ||
        !source_graph_vote_supported(counts[0] + counts[1], 0, 0, kMaxParityP) ||
        !graph_snp_path_supported(gc, seam.left_phase_set, nullptr,
                                  original_graph_phase_sets, stitched_graph_phase_sets, false, true) ||
        !graph_snp_path_supported(gc, seam.right_phase_set, nullptr,
                                  original_graph_phase_sets, stitched_graph_phase_sets, true))
        return std::nullopt;
    std::optional<size_t> first_snp;
    hts_pos_t first_pos = std::numeric_limits<hts_pos_t>::max();
    for (size_t ci = 0; ci < gc.site_meta.size() && ci < chunk.candidates.size(); ++ci) {
        const CandidateVariant& candidate = chunk.candidates[ci];
        if (candidate.phase_set != seam.right_phase_set || candidate.bam_injected ||
            candidate.counts.category != VariantCategory::CleanHetSnp) continue;
        const std::string* alt = selected_graph_candidate_alt(gc, ci);
        if (alt == nullptr) continue;
        const GraphSiteMeta& meta = gc.site_meta[ci];
        const VariantKey key = vcf_to_variant_key(candidate.key.tid, meta.pos, meta.ref, *alt);
        if (key.type != VariantType::Snp || key.ref_len != 1 || key.alt.size() != 1 ||
            key.pos <= right_key.pos || key.pos >= first_pos) continue;
        first_snp = ci;
        first_pos = key.pos;
    }
    if (!first_snp) return std::nullopt;
    const auto edge = local_run_boundary_flip(chunk, *right_i, *first_snp, kMinMapq, kMaxParityP);
    if (!edge || *edge) return std::nullopt;
    return std::make_pair(*left_i, *right_i);
}

// A repeat's net CIGAR length can separate insertion haplotypes when its
// individual bases have low quality. This never changes the candidate allele.
static int physical_repeat_insertion_call(
        const bam1_t* read, const CandidateVariant& insertion,
        WorkerContext& context, int tid,
        int* observed_length = nullptr) {
    constexpr hts_pos_t kRepeatFlank = 16;
    const std::string& alt = insertion.key.alt;
    if (alt.empty() || alt.size() > kLongPhysicalInsertionThreshold)
        return -1;
    size_t period = 1;
    for (; period <= 2; ++period)
        if (alt.size() % period == 0 &&
            std::equal(alt.begin() + period, alt.end(), alt.begin()))
            break;
    if (period > 2) return -1;
    // The reference must contain at least two copies of this motif; length
    // alone is not an allele observation outside a tandem repeat.
    for (size_t offset = 0; offset < 2 * period; ++offset) {
        const char base = context.ref.base(
            tid, insertion.key.pos + static_cast<hts_pos_t>(offset),
            context.primary_header());
        if (std::toupper(static_cast<unsigned char>(base)) !=
            std::toupper(static_cast<unsigned char>(alt[offset % period])))
            return -1;
    }
    const hts_pos_t beg = insertion.key.pos - kRepeatFlank;
    const hts_pos_t end = insertion.key.pos + kRepeatFlank;
    if (read->core.pos + 1 > beg || bam_endpos(read) < end)
        return -1;
    hts_pos_t ref_pos = read->core.pos + 1;
    int net_length = 0;
    int indel_events = 0;
    const uint32_t* cigar = bam_get_cigar(read);
    for (uint32_t ci = 0; ci < read->core.n_cigar; ++ci) {
        const int op = bam_cigar_op(cigar[ci]);
        const int length = bam_cigar_oplen(cigar[ci]);
        if (op == BAM_CINS && ref_pos >= beg && ref_pos <= end) {
            net_length += length;
            ++indel_events;
        }
        if (op == BAM_CDEL && ref_pos < end &&
            ref_pos + length > beg) {
            net_length -= static_cast<int>(
                std::min(ref_pos + length, end) -
                std::max(ref_pos, beg));
            ++indel_events;
        }
        if (bam_cigar_type(op) & 2) ref_pos += length;
    }
    if (indel_events > 1) return -1;
    if (observed_length != nullptr) *observed_length = net_length;
    const int alt_length = static_cast<int>(alt.size());
    if (2 * net_length == alt_length) return -1;
    return 2 * net_length < alt_length ? 0 : 1;
}

static bool source_snp_insertion_link_supported(
        const GraphChunkBuildResult& gc, size_t snp_i, size_t insertion_i,
        WorkerContext& context, int tid, bool allow_one_hap = false);

// Two separated BAM insertion blocks may have no SNP-spanning molecule.
// Use the repeat-length class of each boundary on the same primary reads.
// Require an exact zero-length observation at each boundary as well as two
// agreeing molecules; shorter slippage events never become REF variant rows.
static bool stitch_repeat_insertion_pair(
        GraphChunkBuildResult& gc, WorkerContext& context, int tid,
        const RecoverySeam& seam, hts_pos_t left_ps, hts_pos_t right_ps) {
    constexpr int kMinMapq = 30;
    constexpr int kUnknownQuality = 255;
    constexpr int kMinPairedReads = 2;
    const PhasingChunk& chunk = gc.chunk;
    const auto complete_source = [&gc](hts_pos_t ps) {
        const auto path = gc.recovery_source_path_supported.find(ps);
        return path != gc.recovery_source_path_supported.end() &&
            path->second &&
            (gc.recovery_source_weak_cuts.count(ps) == 0 ||
             gc.recovery_source_weak_cuts.at(ps).empty()) &&
            (gc.recovery_source_quality_cuts.count(ps) == 0 ||
             gc.recovery_source_quality_cuts.at(ps).empty());
    };
    std::optional<size_t> left_i, right_i;
    for (const RecoverySourceSite& source : gc.recovery_source_sites) {
        const size_t ci = source.candidate_index;
        if (ci >= chunk.candidates.size()) continue;
        const CandidateVariant& candidate = chunk.candidates[ci];
        const hts_pos_t pos = candidate.key.sort_pos();
        if (candidate.key.type != VariantType::Insertion ||
            !candidate.msa_verified ||
            candidate.hap_to_cons_alle[1] < 0 ||
            candidate.hap_to_cons_alle[1] > 1 ||
            candidate.hap_to_cons_alle[2] !=
                1 - candidate.hap_to_cons_alle[1] ||
            pos < seam.beg || pos > seam.end)
            continue;
        if (candidate.phase_set == left_ps &&
            (!left_i || pos > chunk.candidates[*left_i].key.sort_pos()))
            left_i = ci;
        if (candidate.phase_set == right_ps &&
            (!right_i || pos < chunk.candidates[*right_i].key.sort_pos()))
            right_i = ci;
    }
    if (!left_i || !right_i) return false;
    const CandidateVariant& left = chunk.candidates[*left_i];
    const CandidateVariant& right = chunk.candidates[*right_i];
    if (left.key.sort_pos() != seam.beg ||
        right.key.sort_pos() != seam.end ||
        left.key.sort_pos() >= right.key.sort_pos())
        return false;
    // A second phased allele at either repeat locus makes a non-target
    // length class ambiguous. In particular, an overlapping deletion cannot
    // be treated as the REF side of an insertion's diploid phase vote.
    for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
        if (ci == *left_i || ci == *right_i) continue;
        const CandidateVariant& other = chunk.candidates[ci];
        const int hap1 = other.hap_to_cons_alle[1];
        if ((other.phase_set != left_ps &&
             other.phase_set != right_ps) ||
            hap1 < 0 || hap1 > 1 ||
            other.hap_to_cons_alle[2] != 1 - hap1)
            continue;
        const CandidateVariant& boundary =
            other.phase_set == left_ps ? left : right;
        const hts_pos_t anchor = boundary.key.sort_pos();
        const hts_pos_t other_end = other.key.pos +
            std::max(1, other.key.ref_len) - 1;
        if (other.key.sort_pos() == anchor ||
            (other.key.type == VariantType::Deletion &&
             other.key.pos <= boundary.key.pos &&
             other_end >= boundary.key.pos))
            return false;
    }
    std::unique_ptr<hts_itr_t, decltype(&hts_itr_destroy)> iterator(
        sam_itr_queryi(context.indexes.front().get(), tid,
                       left.key.pos - 1, right.key.pos + 1),
        &hts_itr_destroy);
    if (!iterator) return false;
    std::unique_ptr<bam1_t, AlignmentDeleter> alignment(bam_init1());
    if (!alignment) return false;
    std::unordered_set<std::string> seen;
    std::unordered_set<std::string> right_callable;
    std::array<int, 2> votes{};
    std::array<bool, 2> exact_ref{};
    while (sam_itr_next(context.bams.front()->get(), iterator.get(),
                        alignment.get()) >= 0) {
        const bam1_t* read = alignment.get();
        if ((read->core.flag & (BAM_FUNMAP | BAM_FSECONDARY |
                                BAM_FSUPPLEMENTARY | BAM_FDUP | BAM_FQCFAIL)) ||
            read->core.qual < kMinMapq ||
            read->core.qual == kUnknownQuality ||
            seen.count(bam_get_qname(read)) != 0)
            continue;
        int left_length = 0, right_length = 0;
        const int right_call = physical_repeat_insertion_call(
            read, right, context, tid, &right_length);
        if (right_call < 0) continue;
        right_callable.insert(bam_get_qname(read));
        const int left_call = physical_repeat_insertion_call(
            read, left, context, tid, &left_length);
        if (left_call < 0) continue;
        seen.insert(bam_get_qname(read));
        exact_ref[0] = exact_ref[0] || left_length == 0;
        exact_ref[1] = exact_ref[1] || right_length == 0;
        const bool left_hap1 =
            (left_call == 1) == (left.hap_to_cons_alle[1] == 1);
        const bool right_hap1 =
            (right_call == 1) == (right.hap_to_cons_alle[1] == 1);
        ++votes[left_hap1 != right_hap1 ? 1 : 0];
    }
    if (votes[0] + votes[1] < kMinPairedReads ||
        (votes[0] > 0 && votes[1] > 0) ||
        !exact_ref[0] || !exact_ref[1])
        return false;
    if (complete_source(left_ps) && complete_source(right_ps))
        return merge_phase_sets_in_place(gc.chunk, left_ps, right_ps,
                                         votes[1] > 0);

    // Only the boundary-side source run is justified by these two reads.
    // A weak edge earlier in the left source must retain its old phase set.
    const auto left_cuts = gc.recovery_source_weak_cuts.find(left_ps);
    const auto right_cuts = gc.recovery_source_weak_cuts.find(right_ps);
    if (left_cuts == gc.recovery_source_weak_cuts.end() ||
        left_cuts->second.empty() ||
        right_cuts == gc.recovery_source_weak_cuts.end() ||
        right_cuts->second.size() != 1 ||
        right_cuts->second.front() != right.key.sort_pos() ||
        (gc.recovery_source_quality_cuts.count(left_ps) != 0 &&
         !gc.recovery_source_quality_cuts.at(left_ps).empty()) ||
        (gc.recovery_source_quality_cuts.count(right_ps) != 0 &&
         !gc.recovery_source_quality_cuts.at(right_ps).empty()))
        return false;
    const hts_pos_t cut = *std::max_element(
        left_cuts->second.begin(), left_cuts->second.end());
    if (cut >= left.key.sort_pos()) return false;
    size_t right_sites = 0;
    hts_pos_t suffix_beg = left.key.sort_pos();
    for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
        const CandidateVariant& candidate = chunk.candidates[ci];
        const int hap1 = candidate.hap_to_cons_alle[1];
        const int hap2 = candidate.hap_to_cons_alle[2];
        const bool oriented = hap1 >= 0 && hap1 <= 1 &&
            hap2 == 1 - hap1;
        if (candidate.phase_set == right_ps && oriented)
            ++right_sites;
        if (candidate.phase_set != left_ps ||
            candidate.key.sort_pos() <= cut || !oriented)
            continue;
        const bool belongs_to_source = std::any_of(
            gc.recovery_source_sites.begin(), gc.recovery_source_sites.end(),
            [ci, left_ps](const RecoverySourceSite& source) {
                return source.candidate_index == ci &&
                       source.phase_set == left_ps;
            });
        if (!belongs_to_source) return false;
        suffix_beg = std::min(suffix_beg, candidate.key.sort_pos());
    }
    if (right_sites != 1 || suffix_beg <= cut) return false;
    // Graph-projected read starts need not match the primary BAM alignment.
    // Keep only reads with an actual callable boundary allele in the joined
    // gauge. Others retain their HP assignment in an independent read-only
    // phase set so the weak right edge cannot become a whole-block switch.
    std::vector<size_t> right_tail_reads;
    for (size_t ri = 0; ri < chunk.reads.size(); ++ri)
        if (ri < chunk.phase_sets.size() &&
            chunk.phase_sets[ri] == right_ps &&
            right_callable.count(chunk.reads[ri].qname) == 0)
            right_tail_reads.push_back(ri);
    hts_pos_t detached_ps = right.key.pos + 1;
    const auto used_ps = [&chunk](hts_pos_t ps) {
        return std::find(chunk.phase_sets.begin(), chunk.phase_sets.end(),
                         ps) != chunk.phase_sets.end() ||
            std::find(chunk.gap_phase_sets.begin(),
                      chunk.gap_phase_sets.end(), ps) !=
                chunk.gap_phase_sets.end() ||
            std::any_of(chunk.candidates.begin(), chunk.candidates.end(),
                        [ps](const CandidateVariant& candidate) {
                            return candidate.phase_set == ps;
                        });
    };
    while (used_ps(detached_ps)) ++detached_ps;
    // The source may mark a one-haplotype edge weak even when two clean,
    // independent molecules directly confirm its inherited orientation.
    // Certifying that edge keeps the established left read block intact.
    bool whole_left_certified = false;
    if (left_cuts->second.size() == 1) {
        std::optional<size_t> cut_snp, suffix_insertion;
        for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
            const CandidateVariant& candidate = chunk.candidates[ci];
            if (candidate.phase_set != left_ps) continue;
            if (candidate.key.type == VariantType::Snp &&
                candidate.key.sort_pos() == cut &&
                candidate.counts.category == VariantCategory::CleanHetSnp)
                cut_snp = ci;
            if (candidate.key.type == VariantType::Insertion &&
                candidate.key.sort_pos() == suffix_beg)
                suffix_insertion = ci;
        }
        whole_left_certified = cut_snp && suffix_insertion &&
            source_snp_insertion_link_supported(
                gc, *cut_snp, *suffix_insertion, context, tid, true);
    }
    const bool joined = whole_left_certified ?
        merge_phase_sets_in_place(gc.chunk, left_ps, right_ps,
                                  votes[1] > 0) :
        join_graph_suffix_after_weak_edge(
            gc, context, tid, left_ps, right_ps, votes[1] > 0,
            *left_i, {cut, suffix_beg}, nullptr, nullptr, true);
    if (!joined) return false;
    for (const size_t ri : right_tail_reads)
        gc.chunk.phase_sets[ri] = detached_ps;

    return true;
}

// A clean left SNP can orient a right BAM insertion when the seam contains
// no clean right SNP. Long insertions may have only ALT-spanning molecules;
// require two such reads, accurate inserted bases, and a complete BAM source.
static bool stitch_snp_to_msa_insertion(
        GraphChunkBuildResult& gc, WorkerContext& context, int tid,
        const RecoverySeam& seam, hts_pos_t left_ps, hts_pos_t right_ps,
        hts_pos_t left_pos, char left_ref, char left_alt, size_t left_i,
        const std::map<std::string, hts_pos_t>* original_graph_phase_sets,
        const std::map<std::string, hts_pos_t>* stitched_graph_phase_sets) {
    constexpr int kMinMapq = 30;
    constexpr int kLongInsertionMinMapq = 5;
    constexpr int kMinBaseq = 30;
    constexpr int kLongInsertionMinBaseq = 20;
    constexpr int kUnknownQuality = 255;
    constexpr int kMinSupportPerAllele = 2;
    constexpr double kMaxWrongParity = 0.001;
    constexpr double kLongInsertionMaxWrongParity = 0.01;
    PhasingChunk& chunk = gc.chunk;
    std::optional<size_t> insertion_i;
    hts_pos_t insertion_pos = std::numeric_limits<hts_pos_t>::max();
    for (const RecoverySourceSite& source : gc.recovery_source_sites) {
        if (source.candidate_index >= chunk.candidates.size()) continue;
        const CandidateVariant& candidate =
            chunk.candidates[source.candidate_index];
        if (candidate.phase_set != right_ps ||
            candidate.key.type != VariantType::Insertion ||
            (candidate.counts.category != VariantCategory::CleanHetIndel &&
             candidate.counts.category != VariantCategory::NoisyCandHet) ||
            !candidate.msa_verified || candidate.hap_to_cons_alle[1] < 0 ||
            candidate.hap_to_cons_alle[1] > 1 ||
            candidate.hap_to_cons_alle[2] !=
                1 - candidate.hap_to_cons_alle[1] ||
            candidate.key.pos <= left_pos ||
            candidate.key.sort_pos() > seam.end ||
            candidate.key.pos >= insertion_pos)
            continue;
        insertion_i = source.candidate_index;
        insertion_pos = candidate.key.pos;
    }
    if (!insertion_i) return false;
    const CandidateVariant& insertion = chunk.candidates[*insertion_i];
    const bool long_insertion =
        insertion.key.alt.size() > kLongPhysicalInsertionThreshold;
    const int min_mapq = long_insertion ? kLongInsertionMinMapq : kMinMapq;
    std::unique_ptr<hts_itr_t, decltype(&hts_itr_destroy)> iterator(
        sam_itr_queryi(context.indexes.front().get(), tid,
                       left_pos - 1, insertion_pos),
        &hts_itr_destroy);
    if (!iterator) return false;
    std::unique_ptr<bam1_t, AlignmentDeleter> alignment(bam_init1());
    if (!alignment) return false;
    std::unordered_set<std::string> seen;
    std::array<int, 2> allele_support{};
    double log_odds = 0.0;
    while (sam_itr_next(context.bams.front()->get(), iterator.get(),
                        alignment.get()) >= 0) {
        const bam1_t* read = alignment.get();
        if ((read->core.flag & (BAM_FUNMAP | BAM_FSECONDARY |
                                BAM_FSUPPLEMENTARY | BAM_FDUP | BAM_FQCFAIL)) ||
            read->core.qual < min_mapq ||
            read->core.qual == kUnknownQuality ||
            seen.count(bam_get_qname(read)) != 0)
            continue;
        int left_quality = 0;
        const int left_call = physical_snp_call(
            read, left_pos, left_ref, left_alt, &left_quality);
        if ((left_call != 0 && left_call != 2) ||
            left_quality < kMinBaseq || left_quality == kUnknownQuality)
            continue;
        int insertion_quality = kMinBaseq;
        const int insertion_call = physical_equivalent_insertion_call(
            read, insertion, context, tid, kMinBaseq,
            &insertion_quality,
            long_insertion ? kLongInsertionMinBaseq : kMinBaseq);
        if (insertion_call != 0 && insertion_call != 1) continue;
        seen.insert(bam_get_qname(read));
        ++allele_support[insertion_call];
        const bool left_hap1 = (left_call == 2) ==
            (chunk.candidates[left_i].hap_to_cons_alle[1] == 1);
        const bool right_hap1 = (insertion_call == 1) ==
            (insertion.hap_to_cons_alle[1] == 1);
        const bool flip = left_hap1 != right_hap1;
        const double p = std::pow(10.0, -left_quality / 10.0) +
            std::pow(10.0, -insertion_quality / 10.0) +
            2.0 * std::pow(10.0, -read->core.qual / 10.0);
        if (p > 0.0 && p < 0.5)
            log_odds += (flip ? 1.0 : -1.0) *
                std::log((1.0 - p) / p);
    }
    const double max_wrong_parity = long_insertion ?
        kLongInsertionMaxWrongParity : kMaxWrongParity;
    const double threshold = std::log((1.0 - max_wrong_parity) /
                                      max_wrong_parity);
    // Earlier transfers can relabel the insertion. Only its original BAM
    // source and a consistent current allele gauge certify that source path.
    const bool complete_right_source = long_insertion &&
        bam_source_site_path_supported(gc, *insertion_i);
    if ((!long_insertion && allele_support[0] < kMinSupportPerAllele) ||
        allele_support[1] < kMinSupportPerAllele ||
        std::abs(log_odds) < threshold ||
        (!complete_right_source &&
         !graph_snp_path_supported(gc, right_ps, nullptr,
                                   original_graph_phase_sets,
                                   stitched_graph_phase_sets)))
        return false;
    // A repeat insertion can carry a misleading physical length class.
    // Prefer a decisive clean-SNP pair when the BAM source has a supported
    // path from this insertion to that SNP. A weak cut or disagreeing SNPs
    // prevent the substitution of one relation for the other.
    constexpr hts_pos_t kMaxConfirmingSnpDistance = 10000;
    hts_pos_t insertion_source_ps = 0;
    for (const RecoverySourceSite& site : gc.recovery_source_sites)
        if (site.candidate_index == *insertion_i) {
            insertion_source_ps = site.phase_set;
            break;
        }
    std::optional<bool> clean_snp_flip;
    for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
        const CandidateVariant& snp = chunk.candidates[ci];
        if (snp.phase_set != right_ps || !snp.bam_injected ||
            !snp.alignment_verified ||
            snp.counts.category != VariantCategory::CleanHetSnp ||
            snp.key.type != VariantType::Snp ||
            snp.key.pos <= insertion_pos ||
            snp.key.pos - insertion_pos > kMaxConfirmingSnpDistance)
            continue;
        const auto source = std::find_if(
            gc.recovery_source_sites.begin(),
            gc.recovery_source_sites.end(),
            [ci](const RecoverySourceSite& site) {
                return site.candidate_index == ci;
            });
        if (insertion_source_ps <= 0 ||
            source == gc.recovery_source_sites.end() ||
            source->phase_set != insertion_source_ps)
            continue;
        const auto cuts = gc.recovery_source_weak_cuts.find(
            insertion_source_ps);
        if (cuts == gc.recovery_source_weak_cuts.end() ||
            std::any_of(cuts->second.begin(), cuts->second.end(),
                        [&insertion, &snp](hts_pos_t cut) {
                            return insertion.key.sort_pos() <= cut &&
                                cut < snp.key.sort_pos();
                        }))
            continue;
        const std::optional<bool> flip = local_run_boundary_flip(
            chunk, left_i, ci, kMinMapq, kMaxWrongParity);
        if (!flip) continue;
        if (clean_snp_flip && *clean_snp_flip != *flip) return false;
        clean_snp_flip = *flip;
    }
    const bool selected_flip = clean_snp_flip.value_or(log_odds > 0.0);
    if (graph_snp_path_supported(gc, left_ps, nullptr,
                                 original_graph_phase_sets,
                                 stitched_graph_phase_sets, false, true))
        return merge_phase_sets_in_place(chunk, left_ps, right_ps,
                                         selected_flip);
    if (!long_insertion) return false;
    std::pair<hts_pos_t, hts_pos_t> cut;
    graph_snp_path_supported(gc, left_ps, &cut,
                             original_graph_phase_sets,
                             stitched_graph_phase_sets, false, true, 0, true);
    // The left phase set was already joined before this seam. Verify its
    // boundary-side SNP suffix, then orient only the new right source; splitting
    // the established left block would undo its previously supported join.
    // The suffix check does not change the selected clean-SNP parity.
    return cut.first > 0 && cut.first < cut.second &&
        graph_snp_path_supported(
            gc, left_ps, nullptr, original_graph_phase_sets,
            stitched_graph_phase_sets, false, false, cut.second) &&
        merge_phase_sets_in_place(chunk, left_ps, right_ps,
                                  selected_flip);
}

// A recovered deletion can be the right boundary of a graph seam even when
// no right SNP is available. Score the physical left SNP and the exact BAM
// deletion allele together; a single read is usable only when its combined
// allele and mapping error passes the wrong-parity bound.
static bool stitch_snp_to_msa_deletion(
        GraphChunkBuildResult& gc, WorkerContext& context, int tid,
        const RecoverySeam& seam, hts_pos_t left_ps, hts_pos_t right_ps,
        hts_pos_t left_pos, char left_ref, char left_alt, size_t left_i,
        const std::map<std::string, hts_pos_t>* original_graph_phase_sets,
        const std::map<std::string, hts_pos_t>* stitched_graph_phase_sets) {
    constexpr int kMinMapq = 30;
    constexpr int kMinSnpBaseq = 10;
    constexpr int kMinDeletionBaseq = 30;
    constexpr int kUnknownQuality = 255;
    constexpr double kMaxWrongParity = 0.05;
    PhasingChunk& chunk = gc.chunk;
    std::optional<size_t> deletion_i;
    for (const RecoverySourceSite& source : gc.recovery_source_sites) {
        if (source.candidate_index >= chunk.candidates.size()) continue;
        const CandidateVariant& candidate =
            chunk.candidates[source.candidate_index];
        if (candidate.phase_set != right_ps ||
            candidate.key.type != VariantType::Deletion ||
            candidate.key.sort_pos() <= left_pos ||
            candidate.key.sort_pos() > seam.end ||
            (candidate.counts.category != VariantCategory::CleanHetIndel &&
             candidate.counts.category != VariantCategory::NoisyCandHet) ||
            !candidate.msa_verified || !candidate.alignment_verified ||
            candidate.hap_to_cons_alle[1] < 0 ||
            candidate.hap_to_cons_alle[1] > 1 ||
            candidate.hap_to_cons_alle[2] !=
                1 - candidate.hap_to_cons_alle[1])
            continue;
        if (deletion_i) return false;
        deletion_i = source.candidate_index;
    }
    if (!deletion_i) return false;
    const CandidateVariant& deletion = chunk.candidates[*deletion_i];
    // Transfer can move this deletion into a different BAM block. Its new
    // PS does not identify its original path or certify its allele gauge.
    if (!bam_source_site_path_supported(gc, *deletion_i) ||
        !graph_snp_path_supported(gc, left_ps, nullptr,
                                  original_graph_phase_sets,
                                  stitched_graph_phase_sets, false, true) ||
        !graph_snp_path_supported(gc, right_ps, nullptr,
                                  original_graph_phase_sets,
                                  stitched_graph_phase_sets, false, true))
        return false;
    for (const RecoverySourceSite& source : gc.recovery_source_sites) {
        if (source.candidate_index == *deletion_i ||
            source.candidate_index >= chunk.candidates.size()) continue;
        const CandidateVariant& other =
            chunk.candidates[source.candidate_index];
        if (other.phase_set == right_ps &&
            other.key.type == VariantType::Deletion &&
            other.msa_verified &&
            other.hap_to_cons_alle[1] ==
                1 - deletion.hap_to_cons_alle[1] &&
            std::max(other.key.pos, deletion.key.pos) <
                std::min(other.key.pos + other.key.ref_len,
                         deletion.key.pos + deletion.key.ref_len))
            return false;
    }
    std::unique_ptr<hts_itr_t, decltype(&hts_itr_destroy)> iterator(
        sam_itr_queryi(context.indexes.front().get(), tid,
                       left_pos - 1,
                       deletion.key.pos + deletion.key.ref_len),
        &hts_itr_destroy);
    if (!iterator) return false;
    std::unique_ptr<bam1_t, AlignmentDeleter> alignment(bam_init1());
    if (!alignment) return false;
    std::unordered_set<std::string> seen;
    double log_odds = 0.0;
    while (sam_itr_next(context.bams.front()->get(), iterator.get(),
                        alignment.get()) >= 0) {
        const bam1_t* read = alignment.get();
        if ((read->core.flag & (BAM_FUNMAP | BAM_FSECONDARY |
                                BAM_FSUPPLEMENTARY | BAM_FDUP | BAM_FQCFAIL)) ||
            read->core.qual < kMinMapq ||
            read->core.qual == kUnknownQuality ||
            seen.count(bam_get_qname(read)) != 0)
            continue;
        int left_quality = 0;
        const int left_call = physical_snp_call(
            read, left_pos, left_ref, left_alt, &left_quality);
        if ((left_call != 0 && left_call != 2) ||
            left_quality < kMinSnpBaseq ||
            left_quality == kUnknownQuality)
            continue;
        const int deletion_call = physical_equivalent_deletion_call(
            read, deletion, context, tid, kMinDeletionBaseq);
        if (deletion_call != 0 && deletion_call != 1) continue;
        seen.insert(bam_get_qname(read));
        const bool left_hap1 = (left_call == 2) ==
            (chunk.candidates[left_i].hap_to_cons_alle[1] == 1);
        const bool right_hap1 = (deletion_call == 1) ==
            (deletion.hap_to_cons_alle[1] == 1);
        const double p = std::pow(10.0, -left_quality / 10.0) +
            std::pow(10.0, -kMinDeletionBaseq / 10.0) +
            2.0 * std::pow(10.0, -read->core.qual / 10.0);
        if (p > 0.0 && p < 0.5)
            log_odds += (left_hap1 != right_hap1 ? 1.0 : -1.0) *
                std::log((1.0 - p) / p);
    }
    const double threshold = std::log((1.0 - kMaxWrongParity) /
                                      kMaxWrongParity);
    return !seen.empty() && std::abs(log_odds) >= threshold &&
        merge_phase_sets_in_place(chunk, left_ps, right_ps,
                                  log_odds > 0.0);
}

// A direct molecule link can certify a BAM source's last weak cut against
// a clean SNP in the retained suffix. Ordinarily require both insertion
// alleles; an explicit one-haplotype certificate needs two independent
// agreeing molecules and the same stringent parity bound.
static bool source_snp_insertion_link_supported(
        const GraphChunkBuildResult& gc, size_t snp_i, size_t insertion_i,
        WorkerContext& context, int tid, bool allow_one_hap) {
    constexpr int kMinMapq = 30;
    constexpr int kMinBaseq = 30;
    constexpr int kUnknownQuality = 255;
    constexpr int kMinSupportPerAllele = 2;
    constexpr double kMaxWrongParity = 0.001;
    const PhasingChunk& chunk = gc.chunk;
    const CandidateVariant& snp = chunk.candidates[snp_i];
    const CandidateVariant& insertion = chunk.candidates[insertion_i];
    VariantKey physical_snp = snp.key;
    if (!snp.bam_injected && snp_i < gc.site_meta.size()) {
        const std::string* alt = selected_graph_candidate_alt(gc, snp_i);
        if (alt == nullptr) return false;
        const GraphSiteMeta& meta = gc.site_meta[snp_i];
        physical_snp = vcf_to_variant_key(snp.key.tid,
                                           meta.pos, meta.ref, *alt);
    }
    if (physical_snp.type != VariantType::Snp ||
        physical_snp.alt.size() != 1 ||
        snp.counts.category != VariantCategory::CleanHetSnp ||
        snp.hap_to_cons_alle[1] < 0 || snp.hap_to_cons_alle[1] > 1 ||
        snp.hap_to_cons_alle[2] != 1 - snp.hap_to_cons_alle[1] ||
        physical_snp.pos >= insertion.key.pos)
        return false;
    const char ref = context.ref.base(
        tid, physical_snp.pos, context.primary_header());
    if (ref == 'N') return false;
    std::unique_ptr<hts_itr_t, decltype(&hts_itr_destroy)> iterator(
        sam_itr_queryi(context.indexes.front().get(), tid,
                       physical_snp.pos - 1, insertion.key.pos), &hts_itr_destroy);
    if (!iterator) return false;
    std::unique_ptr<bam1_t, AlignmentDeleter> alignment(bam_init1());
    if (!alignment) return false;
    std::unordered_set<std::string> seen;
    std::array<int, 2> insertion_support{};
    double log_odds = 0.0;
    while (sam_itr_next(context.bams.front()->get(), iterator.get(),
                        alignment.get()) >= 0) {
        const bam1_t* read = alignment.get();
        if ((read->core.flag & (BAM_FUNMAP | BAM_FSECONDARY |
                                BAM_FSUPPLEMENTARY | BAM_FDUP | BAM_FQCFAIL)) ||
            read->core.qual < kMinMapq || read->core.qual == kUnknownQuality ||
            seen.count(bam_get_qname(read)) != 0)
            continue;
        int snp_quality = 0;
        const int snp_call = physical_snp_call(
            read, physical_snp.pos, ref, physical_snp.alt[0],
            &snp_quality);
        if ((snp_call != 0 && snp_call != 2) ||
            snp_quality < kMinBaseq || snp_quality == kUnknownQuality)
            continue;
        const int insertion_call = physical_equivalent_insertion_call(
            read, insertion, context, tid, kMinBaseq);
        if (insertion_call != 0 && insertion_call != 1) continue;
        seen.insert(bam_get_qname(read));
        ++insertion_support[static_cast<size_t>(insertion_call)];
        const bool snp_hap1 = (snp_call == 2) ==
            (snp.hap_to_cons_alle[1] == 1);
        const bool insertion_hap1 = (insertion_call == 1) ==
            (insertion.hap_to_cons_alle[1] == 1);
        const double p = std::pow(10.0, -snp_quality / 10.0) +
            std::pow(10.0, -kMinBaseq / 10.0) +
            2.0 * std::pow(10.0, -read->core.qual / 10.0);
        if (p > 0.0 && p < 0.5)
            log_odds += (snp_hap1 != insertion_hap1 ? 1.0 : -1.0) *
                std::log((1.0 - p) / p);
    }
    const double threshold = std::log((1.0 - kMaxWrongParity) /
                                      kMaxWrongParity);
    if (allow_one_hap)
        return ((insertion_support[0] >= kMinSupportPerAllele &&
                 insertion_support[1] == 0) ||
                (insertion_support[1] >= kMinSupportPerAllele &&
                 insertion_support[0] == 0)) &&
            log_odds <= -threshold;
    return insertion_support[0] >= kMinSupportPerAllele &&
           insertion_support[1] >= kMinSupportPerAllele &&
           log_odds <= -threshold;
}

// Certify one weak graph SNP edge from the original BAM bases. A one-haplotype
// GAF edge can still be reliable when two independent molecules call both
// physical SNPs cleanly and agree with the block's existing orientation.
static std::optional<size_t> physical_graph_snp_edge_supported(
        const GraphChunkBuildResult& gc, WorkerContext& context, int tid,
        hts_pos_t phase_set, const std::pair<hts_pos_t, hts_pos_t>& cut,
        const std::map<std::string, hts_pos_t>* original_graph_phase_sets,
        const std::map<std::string, hts_pos_t>* stitched_graph_phase_sets,
        bool require_significant_votes = false) {
    constexpr int kMinMapq = 30;
    constexpr int kMinBaseq = 30;
    constexpr int kMinPairedReads = 2;
    constexpr int kUnknownQuality = 255;
    constexpr double kMaxWrongParity = 0.001;
    const PhasingChunk& chunk = gc.chunk;
    std::optional<size_t> left_i, right_i;
    VariantKey left_key, right_key;
    int suffix_snps = 0;
    for (size_t ci = 0; ci < chunk.candidates.size() &&
                        ci < gc.site_meta.size(); ++ci) {
        const CandidateVariant& candidate = chunk.candidates[ci];
        if (candidate.phase_set != phase_set || candidate.bam_injected ||
            candidate.counts.category != VariantCategory::CleanHetSnp ||
            candidate.hap_to_cons_alle[1] < 0 ||
            candidate.hap_to_cons_alle[1] > 1 ||
            candidate.hap_to_cons_alle[2] !=
                1 - candidate.hap_to_cons_alle[1])
            continue;
        const std::string* alt = selected_graph_candidate_alt(gc, ci);
        if (alt == nullptr) continue;
        const GraphSiteMeta& meta = gc.site_meta[ci];
        const VariantKey key = vcf_to_variant_key(
            candidate.key.tid, meta.pos, meta.ref, *alt);
        if (key.type != VariantType::Snp || key.ref_len != 1 ||
            key.alt.size() != 1) continue;
        if (key.pos >= cut.second) ++suffix_snps;
        if (key.pos == cut.first) {
            if (left_i) return std::nullopt;
            left_i = ci;
            left_key = key;
        }
        if (key.pos == cut.second) {
            if (right_i) return std::nullopt;
            right_i = ci;
            right_key = key;
        }
    }
    if (!left_i || !right_i || suffix_snps == 0 ||
        (suffix_snps > 1 && !graph_snp_path_supported(
            gc, phase_set, nullptr, original_graph_phase_sets,
            stitched_graph_phase_sets, false, false, cut.second)))
        return std::nullopt;
    std::unique_ptr<hts_itr_t, decltype(&hts_itr_destroy)> iterator(
        sam_itr_queryi(context.indexes.front().get(), tid,
                       cut.first - 1, cut.second), &hts_itr_destroy);
    if (!iterator) return std::nullopt;
    std::unique_ptr<bam1_t, AlignmentDeleter> alignment(bam_init1());
    if (!alignment) return std::nullopt;
    std::unordered_set<std::string> seen;
    int paired = 0;
    int same = 0, cross = 0;
    std::array<int, 2> supporting_haps{};
    double log_odds = 0.0;
    while (sam_itr_next(context.bams.front()->get(), iterator.get(),
                        alignment.get()) >= 0) {
        const bam1_t* read = alignment.get();
        if ((read->core.flag & (BAM_FUNMAP | BAM_FSECONDARY |
                                BAM_FSUPPLEMENTARY | BAM_FDUP | BAM_FQCFAIL)) ||
            read->core.qual < kMinMapq ||
            read->core.qual == kUnknownQuality ||
            seen.count(bam_get_qname(read)) != 0)
            continue;
        int left_quality = 0, right_quality = 0;
        const int left_call = physical_snp_call(
            read, cut.first, context.ref.base(
                tid, cut.first, context.primary_header()),
            left_key.alt[0], &left_quality);
        const int right_call = physical_snp_call(
            read, cut.second, context.ref.base(
                tid, cut.second, context.primary_header()),
            right_key.alt[0], &right_quality);
        if ((left_call != 0 && left_call != 2) ||
            (right_call != 0 && right_call != 2) ||
            left_quality < kMinBaseq || right_quality < kMinBaseq ||
            left_quality == kUnknownQuality ||
            right_quality == kUnknownQuality)
            continue;
        seen.insert(bam_get_qname(read));
        const bool left_hap1 = (left_call == 2) ==
            (chunk.candidates[*left_i].hap_to_cons_alle[1] == 1);
        const bool right_hap1 = (right_call == 2) ==
            (chunk.candidates[*right_i].hap_to_cons_alle[1] == 1);
        if (left_hap1 != right_hap1) {
            ++cross;
        } else {
            ++same;
            ++supporting_haps[left_hap1 ? 0 : 1];
        }
        const double p = std::pow(10.0, -left_quality / 10.0) +
            std::pow(10.0, -right_quality / 10.0) +
            2.0 * std::pow(10.0, -read->core.qual / 10.0);
        if (p <= 0.0 || p >= 0.5) continue;
        log_odds += (left_hap1 == right_hap1 ? 1.0 : -1.0) *
            std::log((1.0 - p) / p);
        ++paired;
    }
    const double threshold = std::log((1.0 - kMaxWrongParity) /
                                      kMaxWrongParity);
    // A lone high-quality cohort is insufficient when this check supplies
    // a new clean-SNP path. Require both haplotypes and a decisive count tail;
    // signed quality odds keep opposing molecules in the error bound.
    return paired >= kMinPairedReads && log_odds >= threshold &&
        ((!require_significant_votes && cross == 0) ||
         (supporting_haps[0] > 0 && supporting_haps[1] > 0 &&
          source_graph_vote_supported(same, cross, 0, kMaxWrongParity))) ?
        right_i : std::nullopt;
}

// An MSA-verified BAM indel can supply physical boundary evidence when the
// selected clean-SNP pair has no callable reads. Keep its original site row
// while orienting it to the next graph SNP from the BAM alignment.
static bool stitch_msa_indel_to_snp(
        GraphChunkBuildResult& gc, WorkerContext& context, int tid,
        const RecoverySeam& seam, hts_pos_t left_ps, hts_pos_t right_ps,
        hts_pos_t right_pos, const std::string& right_ref,
        const std::string& right_alt,
        size_t right_i,
        const std::map<std::string, hts_pos_t>* original_graph_phase_sets,
        const std::map<std::string, hts_pos_t>* stitched_graph_phase_sets) {
    constexpr int kMinMapq = 30;
    constexpr int kMinBaseq = 30;
    constexpr int kMinInsertionBaseq = 10;
    constexpr int kUnknownQuality = 255;
    constexpr double kMaxWrongParity = 0.001;
    constexpr int kMinSupportPerAllele = 2;
    PhasingChunk& chunk = gc.chunk;
    std::optional<size_t> indel_i;
    hts_pos_t indel_pos = 0;
    for (const RecoverySourceSite& source : gc.recovery_source_sites) {
        if (source.candidate_index >= chunk.candidates.size()) continue;
        const CandidateVariant& candidate =
            chunk.candidates[source.candidate_index];
        if (candidate.phase_set != left_ps ||
            (candidate.key.type != VariantType::Insertion &&
             candidate.key.type != VariantType::Deletion) ||
            (candidate.counts.category != VariantCategory::CleanHetIndel &&
             candidate.counts.category != VariantCategory::NoisyCandHet) ||
            !candidate.msa_verified || candidate.hap_to_cons_alle[1] < 0 ||
            candidate.hap_to_cons_alle[1] > 1 ||
            candidate.hap_to_cons_alle[2] !=
                1 - candidate.hap_to_cons_alle[1] ||
            candidate.key.sort_pos() < seam.beg ||
            candidate.key.pos >= right_pos ||
            candidate.key.pos <= indel_pos)
            continue;
        indel_i = source.candidate_index;
        indel_pos = candidate.key.pos;
    }
    if (!indel_i) return false;
    const CandidateVariant& indel = chunk.candidates[*indel_i];
    // A graph chunk can end directly after the first SNP of the right block.
    // There is then no internal edge to validate, but the sole row can still
    // be oriented by molecules crossing this boundary.
    const bool right_singleton =
        !chunk.candidates[right_i].bam_injected &&
        std::count_if(chunk.candidates.begin(), chunk.candidates.end(),
                      [right_ps](const CandidateVariant& candidate) {
                          return candidate.phase_set == right_ps;
                      }) == 1;
    if (indel.key.type == VariantType::Deletion) {
        // Opposite-haplotype overlapping deletion rows are separate alleles;
        // their REF calls are ambiguous and use the pair bridge above.
        for (const RecoverySourceSite& source : gc.recovery_source_sites) {
            if (source.candidate_index == *indel_i ||
                source.candidate_index >= chunk.candidates.size())
                continue;
            const CandidateVariant& other =
                chunk.candidates[source.candidate_index];
            if (other.phase_set == left_ps &&
                other.key.type == VariantType::Deletion &&
                other.msa_verified &&
                other.hap_to_cons_alle[1] ==
                    1 - indel.hap_to_cons_alle[1] &&
                std::max(other.key.pos, indel.key.pos) <
                    std::min(other.key.pos + other.key.ref_len,
                             indel.key.pos + indel.key.ref_len))
                return false;
        }
    }
    std::unique_ptr<hts_itr_t, decltype(&hts_itr_destroy)> iterator(
        sam_itr_queryi(context.indexes.front().get(), tid,
                       indel.key.pos - 1, right_pos),
        &hts_itr_destroy);
    if (!iterator) return false;
    std::unique_ptr<bam1_t, AlignmentDeleter> alignment(bam_init1());
    if (!alignment) return false;
    std::unordered_set<std::string> seen;
    std::array<int, 2> allele_support{};
    double log_odds = 0.0;
    while (sam_itr_next(context.bams.front()->get(), iterator.get(),
                        alignment.get()) >= 0) {
        const bam1_t* read = alignment.get();
        if ((read->core.flag & (BAM_FUNMAP | BAM_FSECONDARY |
                                BAM_FSUPPLEMENTARY | BAM_FDUP | BAM_FQCFAIL)) ||
            read->core.qual < kMinMapq ||
            read->core.qual == kUnknownQuality ||
            seen.count(bam_get_qname(read)) != 0)
            continue;
        int right_quality = 0;
        const int right_call = physical_substitution_call(
            read, right_pos, right_ref, right_alt, &right_quality);
        if ((right_call != 0 && right_call != 2) ||
            right_quality < kMinBaseq ||
            right_quality == kUnknownQuality)
            continue;
        int indel_call = -1;
        int indel_quality = kMinBaseq;
        if (indel.key.type == VariantType::Deletion) {
            indel_call = physical_equivalent_deletion_call(
                read, indel, context, tid, kMinBaseq);
        } else {
            indel_call = physical_equivalent_insertion_call(
                read, indel, context, tid, kMinBaseq,
                &indel_quality, kMinInsertionBaseq);
        }
        if (indel_call != 0 && indel_call != 1) continue;
        seen.insert(bam_get_qname(read));
        ++allele_support[indel_call];
        const bool left_hap1 = (indel_call == 1) ==
            (indel.hap_to_cons_alle[1] == 1);
        const bool right_hap1 = (right_call == 2) ==
            (chunk.candidates[right_i].hap_to_cons_alle[1] == 1);
        const bool flip = left_hap1 != right_hap1;
        const double p = std::pow(10.0, -indel_quality / 10.0) +
            std::pow(10.0, -right_quality / 10.0) +
            2.0 * std::pow(10.0, -read->core.qual / 10.0);
        if (p > 0.0 && p < 0.5)
            log_odds += (flip ? 1.0 : -1.0) *
                std::log((1.0 - p) / p);
    }
    const double threshold = std::log((1.0 - kMaxWrongParity) /
                                      kMaxWrongParity);
    if ((allele_support[0] == 0 || allele_support[1] == 0) ||
        std::abs(log_odds) < threshold)
        return false;
    std::pair<hts_pos_t, hts_pos_t> right_cut;
    if (!right_singleton && !graph_snp_path_supported(
            gc, right_ps, &right_cut, original_graph_phase_sets,
            stitched_graph_phase_sets, false, false, 0, true)) {
        // The boundary bridge does not certify an internal graph edge. Apply
        // the same independent physical SNP check used on the left; a
        // dominant GAF reversal never reports an eligible cut.
        if (right_cut.first < right_pos ||
            !physical_graph_snp_edge_supported(
                gc, context, tid, right_ps, right_cut,
                original_graph_phase_sets, stitched_graph_phase_sets))
            return false;
    }
    const bool full_allele_support =
        allele_support[0] >= kMinSupportPerAllele &&
        allele_support[1] >= kMinSupportPerAllele;
    std::pair<hts_pos_t, hts_pos_t> weak_cut;
    // Attachment can relabel a cut-free BAM component with a catalog PS.
    // Shared catalog rows retain exact BAM provenance; a later cut elsewhere
    // in their original source must not invalidate this component's extent.
    if (graph_snp_path_supported(
            gc, left_ps, &weak_cut, original_graph_phase_sets,
            stitched_graph_phase_sets,
            indel.key.type == VariantType::Deletion) ||
        (indel.key.type == VariantType::Insertion &&
         bam_source_run_supported(gc, left_ps, true)))
        return full_allele_support &&
            merge_phase_sets_in_place(chunk, left_ps, right_ps,
                                      log_odds > 0.0);
    // Certify the missing graph edge separately from the SNP-to-indel link.
    // Both insertion and deletion boundaries need this complete left path;
    // their decisive cross-seam vote alone cannot certify the earlier edge.
    if (full_allele_support) {
        std::pair<hts_pos_t, hts_pos_t> graph_cut;
        graph_snp_path_supported(gc, left_ps, &graph_cut,
                                 original_graph_phase_sets,
                                 stitched_graph_phase_sets, false, false,
                                 0, true);
        if (graph_cut.first > 0) {
            const std::optional<size_t> suffix_snp =
                physical_graph_snp_edge_supported(
                    gc, context, tid, left_ps, graph_cut,
                    original_graph_phase_sets, stitched_graph_phase_sets);
            if (suffix_snp) {
                const auto deletion_snp = indel.key.type == VariantType::Deletion
                    ? graph_snp_deletion_link_supported(
                        gc, context, tid, left_ps, indel, true)
                    : std::nullopt;
                const bool boundary_link = indel.key.type == VariantType::Insertion
                    ? source_snp_insertion_link_supported(
                        gc, *suffix_snp, *indel_i, context, tid)
                    : deletion_snp && *deletion_snp >= graph_cut.second;
                if (boundary_link)
                    return merge_phase_sets_in_place(
                        chunk, left_ps, right_ps, log_odds > 0.0);
            }
        }
    }
    // A BAM-only block can have a weak edge before its final insertion.
    // Transfer its suffix only when a clean SNP inside that suffix directly
    // confirms the insertion's inherited orientation across the final cut.
    if (indel.key.type == VariantType::Insertion) {
        const auto cuts_it = gc.recovery_source_weak_cuts.find(left_ps);
        if (cuts_it == gc.recovery_source_weak_cuts.end() ||
            cuts_it->second.empty())
            return false;
        const std::vector<hts_pos_t>& cuts = cuts_it->second;
        const hts_pos_t last_cut = *std::max_element(cuts.begin(), cuts.end());
        const hts_pos_t insertion_pos = indel.key.sort_pos();
        if (last_cut >= insertion_pos) return false;
        std::optional<size_t> snp_i;
        hts_pos_t snp_pos = 0;
        for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
            const CandidateVariant& candidate = chunk.candidates[ci];
            const hts_pos_t pos = candidate.key.sort_pos();
            if (candidate.phase_set == left_ps && candidate.bam_injected &&
                candidate.key.type == VariantType::Snp &&
                candidate.counts.category == VariantCategory::CleanHetSnp &&
                pos < last_cut && pos > snp_pos) {
                snp_i = ci;
                snp_pos = pos;
            }
        }
        if (!snp_i || !source_snp_insertion_link_supported(
                gc, *snp_i, *indel_i, context, tid))
            return false;
        hts_pos_t preceding_cut = 0;
        for (const hts_pos_t cut : cuts)
            if (cut < snp_pos) preceding_cut = std::max(preceding_cut, cut);
        if (preceding_cut <= 0) return false;
        return join_graph_suffix_after_weak_edge(
            gc, context, tid, left_ps, right_ps, log_odds > 0.0,
            *indel_i, {preceding_cut, snp_pos},
            original_graph_phase_sets, stitched_graph_phase_sets, true);
    }
    if (indel.key.type != VariantType::Deletion) return false;
    hts_pos_t suffix_start = 0;
    std::pair<hts_pos_t, hts_pos_t> boundary_cut;
    for (size_t attempt = 0; attempt < chunk.candidates.size(); ++attempt) {
        std::pair<hts_pos_t, hts_pos_t> next_cut;
        if (graph_snp_path_supported(
                gc, left_ps, &next_cut, original_graph_phase_sets,
                stitched_graph_phase_sets, false, false,
                suffix_start, true))
            break;
        if (next_cut.first <= 0 || next_cut.second <= suffix_start)
            return false;
        boundary_cut = next_cut;
        suffix_start = next_cut.second;
    }
    if (boundary_cut.first <= 0 || !full_allele_support)
        return false;
    return join_graph_suffix_after_weak_edge(
        gc, context, tid, left_ps, right_ps, log_odds > 0.0,
        *indel_i, boundary_cut, original_graph_phase_sets,
        stitched_graph_phase_sets);
}

// A source block can retain phased variant rows but no read tags after its
// reads attach to a neighboring graph block. Use direct physical allele calls
// to give those rows the same gauge only when the source path and the graph
// block's one weak SNP edge are independently certified.
static void attach_readless_insertion_source_blocks(
        GraphChunkBuildResult& gc, WorkerContext& context, int tid,
        const std::map<std::string, hts_pos_t>* original_graph_phase_sets,
        const std::map<std::string, hts_pos_t>* stitched_graph_phase_sets) {
    constexpr hts_pos_t kMaxSnpContext = 2000;
    PhasingChunk& chunk = gc.chunk;
    std::unordered_set<hts_pos_t> read_phase_sets(
        chunk.phase_sets.begin(), chunk.phase_sets.end());
    for (const RecoverySourceSite& source : gc.recovery_source_sites) {
        if (source.candidate_index >= chunk.candidates.size()) continue;
        const CandidateVariant& insertion =
            chunk.candidates[source.candidate_index];
        const hts_pos_t source_ps = insertion.phase_set;
        const auto complete =
            gc.recovery_source_path_supported.find(source.phase_set);
        if (source_ps <= 0 || read_phase_sets.count(source_ps) != 0 ||
            insertion.key.type != VariantType::Insertion ||
            !insertion.msa_verified ||
            (insertion.counts.category != VariantCategory::CleanHetIndel &&
             insertion.counts.category != VariantCategory::NoisyCandHet) ||
            insertion.hap_to_cons_alle[1] < 0 ||
            insertion.hap_to_cons_alle[1] > 1 ||
            insertion.hap_to_cons_alle[2] !=
                1 - insertion.hap_to_cons_alle[1] ||
            complete == gc.recovery_source_path_supported.end() ||
            !complete->second ||
            (gc.recovery_source_weak_cuts.count(source.phase_set) != 0 &&
             !gc.recovery_source_weak_cuts.at(source.phase_set).empty()) ||
            (gc.recovery_source_quality_cuts.count(source.phase_set) != 0 &&
             !gc.recovery_source_quality_cuts.at(source.phase_set).empty()))
            continue;
        for (size_t ci = 0; ci < chunk.candidates.size() &&
                            ci < gc.site_meta.size(); ++ci) {
            const CandidateVariant& snp = chunk.candidates[ci];
            if (snp.phase_set <= 0 || snp.phase_set == source_ps ||
                snp.bam_injected ||
                snp.counts.category != VariantCategory::CleanHetSnp ||
                snp.key.pos >= insertion.key.pos ||
                insertion.key.pos - snp.key.pos > kMaxSnpContext ||
                !source_snp_insertion_link_supported(
                    gc, ci, source.candidate_index, context, tid))
                continue;
            std::pair<hts_pos_t, hts_pos_t> graph_cut;
            if (graph_snp_path_supported(
                    gc, snp.phase_set, &graph_cut,
                    original_graph_phase_sets,
                    stitched_graph_phase_sets, false, false, 0, true) ||
                graph_cut.first <= 0 ||
                !physical_graph_snp_edge_supported(
                    gc, context, tid, snp.phase_set, graph_cut,
                    original_graph_phase_sets, stitched_graph_phase_sets))
                continue;
            merge_phase_sets_in_place(chunk, snp.phase_set, source_ps, false);
            break;
        }
    }
}

// Direct BAM bases restore the clean-SNP bridge evidence that a targeted
// sub-solve can lose when one boundary is a graph-only or demoted BAM site.
// Whole-block joins require a continuous SNP path on both sides. A detached
// BAM-only run qualifies only within its original source's weak-cut bounds.
static void stitch_physical_allele_seams(
        GraphChunkBuildResult& gc, WorkerContext& context,
        const char* contig,
        const std::map<std::string, hts_pos_t>* original_graph_phase_sets = nullptr,
        std::map<std::string, hts_pos_t>* stitched_graph_phase_sets = nullptr) {
    PhasingChunk& chunk = gc.chunk;
    constexpr int kMinMapq = 30;
    constexpr int kMinBaseq = 30;
    constexpr int kUnknownQuality = 255;
    constexpr double kMaxWrongParity = 0.001;
    struct PhysicalSubstitution {
        hts_pos_t pos = 0;
        std::string ref, alt;
    };
    const auto physical_substitution = [&](size_t ci, bool allow_noisy = false)
            -> std::optional<PhysicalSubstitution> {
        if (ci >= chunk.candidates.size()) return std::nullopt;
        const CandidateVariant& candidate = chunk.candidates[ci];
        if (candidate.phase_set <= 0 ||
            candidate.key.type != VariantType::Snp ||
            (candidate.counts.category != VariantCategory::CleanHetSnp &&
             !(allow_noisy && candidate.counts.category ==
                 VariantCategory::NoisyCandHet)) ||
            candidate.hap_to_cons_alle[1] < 0 ||
            candidate.hap_to_cons_alle[1] > 1 ||
            candidate.hap_to_cons_alle[2] !=
                1 - candidate.hap_to_cons_alle[1])
            return std::nullopt;
        if (ci < gc.site_meta.size() && !gc.site_meta[ci].ref.empty()) {
            const std::string* alt = selected_graph_candidate_alt(gc, ci);
            if (alt == nullptr) return std::nullopt;
            const GraphSiteMeta& meta = gc.site_meta[ci];
            if (meta.ref.size() != alt->size()) return std::nullopt;
            size_t prefix = 0, suffix = 0;
            while (prefix < meta.ref.size() &&
                   std::toupper(static_cast<unsigned char>(meta.ref[prefix])) ==
                   std::toupper(static_cast<unsigned char>((*alt)[prefix])))
                ++prefix;
            if (prefix == meta.ref.size()) return std::nullopt;
            while (suffix + prefix < meta.ref.size() &&
                   std::toupper(static_cast<unsigned char>(
                       meta.ref[meta.ref.size() - suffix - 1])) ==
                   std::toupper(static_cast<unsigned char>(
                       (*alt)[alt->size() - suffix - 1])))
                ++suffix;
            const size_t length = meta.ref.size() - prefix - suffix;
            return PhysicalSubstitution{
                meta.pos + static_cast<hts_pos_t>(prefix),
                meta.ref.substr(prefix, length), alt->substr(prefix, length)};
        }
        const hts_pos_t offset = candidate.key.pos - chunk.ref_beg;
        if (candidate.key.alt.size() != 1 || offset < 0 ||
            static_cast<size_t>(offset) >= chunk.ref_seq.size())
            return std::nullopt;
        return PhysicalSubstitution{
            candidate.key.pos,
            std::string(1, chunk.ref_seq[static_cast<size_t>(offset)]),
            candidate.key.alt};
    };
    const int tid = contig == nullptr ? -1 :
        sam_hdr_name2tid(context.primary_header(), contig);
    if (tid < 0 || context.bams.empty() || context.indexes.empty()) return;
    // A graph block with no read assignments can still have phased SNP rows.
    // Use independent molecules calling multiple clean SNPs on both flanks to
    // orient it, and require a physical path through all its clean SNPs.
    const auto stitch_orphan_snp_block = [&](const RecoverySeam& seam) {
        constexpr int kMinIndependentMolecules = 3;
        constexpr int kMinSnpCallsPerSide = 2;
        constexpr double kMaxWrongParityProbability = 0.001;
        const hts_pos_t left_ps = seam.left_phase_set;
        const hts_pos_t right_ps = seam.right_phase_set;
        if (left_ps <= 0 || right_ps <= 0 || left_ps == right_ps ||
            std::find(chunk.phase_sets.begin(), chunk.phase_sets.end(),
                      left_ps) == chunk.phase_sets.end() ||
            std::find(chunk.phase_sets.begin(), chunk.phase_sets.end(),
                      right_ps) != chunk.phase_sets.end())
            return false;
        struct BridgeSnp {
            hts_pos_t pos;
            char ref;
            char alt;
            int hap1_allele;
        };
        std::array<std::vector<BridgeSnp>, 2> sites;
        for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
            const CandidateVariant& candidate = chunk.candidates[ci];
            const int side = candidate.phase_set == left_ps ? 0 :
                candidate.phase_set == right_ps ? 1 : -1;
            if (side < 0) continue;
            const auto physical = physical_substitution(ci);
            if (!physical || physical->ref.size() != 1 ||
                physical->alt.size() != 1 ||
                (side == 0 && physical->pos > seam.beg) ||
                (side == 1 && physical->pos < seam.end))
                continue;
            sites[side].push_back(BridgeSnp{
                physical->pos, physical->ref[0], physical->alt[0],
                candidate.hap_to_cons_alle[1]});
        }
        for (auto& side_sites : sites) {
            std::sort(side_sites.begin(), side_sites.end(),
                      [](const BridgeSnp& a, const BridgeSnp& b) {
                          return a.pos < b.pos;
                      });
            std::vector<BridgeSnp> unique;
            for (size_t i = 0; i < side_sites.size(); ++i) {
                if ((i > 0 && side_sites[i - 1].pos == side_sites[i].pos) ||
                    (i + 1 < side_sites.size() &&
                     side_sites[i + 1].pos == side_sites[i].pos))
                    continue;
                unique.push_back(side_sites[i]);
            }
            side_sites = std::move(unique);
        }
        if (sites[0].size() < kMinSnpCallsPerSide ||
            sites[1].size() < kMinSnpCallsPerSide ||
            !std::any_of(gc.recovery_windows.begin(),
                         gc.recovery_windows.end(),
                         [&](const RecoverySeam& target) {
                             return target.beg <= seam.beg &&
                                 sites[1].back().pos <= target.end;
                         }))
            return false;
        std::unique_ptr<hts_itr_t, decltype(&hts_itr_destroy)> iterator(
            sam_itr_queryi(context.indexes.front().get(), tid,
                           seam.beg - 1, sites[1].back().pos),
            &hts_itr_destroy);
        if (!iterator) return false;
        std::unique_ptr<bam1_t, AlignmentDeleter> alignment(bam_init1());
        if (!alignment) return false;
        std::unordered_set<std::string> seen;
        std::vector<std::array<int, 2>> path_delta(sites[1].size() + 1);
        std::array<int, 2> parity{};
        double wrong_parity_bound = 1.0;
        while (sam_itr_next(context.bams.front()->get(), iterator.get(),
                            alignment.get()) >= 0) {
            const bam1_t* read = alignment.get();
            if ((read->core.flag & (BAM_FUNMAP | BAM_FSECONDARY |
                                    BAM_FSUPPLEMENTARY | BAM_FDUP |
                                    BAM_FQCFAIL)) ||
                read->core.qual < kMinMapq ||
                read->core.qual == kUnknownQuality ||
                !seen.insert(bam_get_qname(read)).second)
                continue;
            std::array<int, 3> left_votes{}, right_votes{};
            std::vector<std::pair<size_t, int>> right_calls;
            double read_error_bound =
                std::pow(10.0, -static_cast<double>(read->core.qual) / 10.0);
            for (size_t side = 0; side < sites.size(); ++side) {
                const auto& side_sites = sites[side];
                for (size_t si = 0; si < side_sites.size(); ++si) {
                    const BridgeSnp& site = side_sites[si];
                    if (site.pos <= read->core.pos ||
                        site.pos > bam_endpos(read)) continue;
                    int quality = 0;
                    const int call = physical_snp_call(
                        read, site.pos, site.ref, site.alt, &quality);
                    if ((call != 0 && call != 2) || quality < kMinBaseq ||
                        quality == kUnknownQuality) continue;
                    const int hap = (call == 2) ==
                        (site.hap1_allele == 1) ? 1 : 2;
                    (side == 0 ? left_votes : right_votes)[hap]++;
                    if (side == 1) right_calls.emplace_back(si, hap);
                    read_error_bound +=
                        std::pow(10.0, -static_cast<double>(quality) / 10.0);
                }
            }
            for (size_t i = 1; i < right_calls.size(); ++i) {
                const size_t previous = right_calls[i - 1].first;
                const size_t current = right_calls[i].first;
                const size_t vote = right_calls[i - 1].second ==
                    right_calls[i].second ? 0 : 1;
                ++path_delta[previous][vote];
                --path_delta[current][vote];
            }
            if (read->core.pos + 1 > seam.beg ||
                bam_endpos(read) < seam.end ||
                left_votes[1] + left_votes[2] < kMinSnpCallsPerSide ||
                right_votes[1] + right_votes[2] < kMinSnpCallsPerSide ||
                (left_votes[1] != 0 && left_votes[2] != 0) ||
                (right_votes[1] != 0 && right_votes[2] != 0))
                continue;
            const bool flip = (left_votes[1] != 0) !=
                (right_votes[1] != 0);
            ++parity[flip ? 1 : 0];
            wrong_parity_bound *= read_error_bound;
        }
        int consistent = 0, conflicting = 0;
        for (size_t i = 0; i + 1 < sites[1].size(); ++i) {
            consistent += path_delta[i][0];
            conflicting += path_delta[i][1];
            if (consistent == 0 || conflicting != 0) return false;
        }
        if (parity[0] + parity[1] < kMinIndependentMolecules ||
            (parity[0] != 0 && parity[1] != 0) ||
            wrong_parity_bound > kMaxWrongParityProbability)
            return false;
        return merge_phase_sets_in_place(chunk, left_ps, right_ps,
                                         parity[1] != 0);
    };
    std::vector<RecoverySeam> stitch_seams = gc.recovery_windows;
    std::set<std::pair<hts_pos_t, hts_pos_t>> seen_pairs;
    std::set<std::pair<hts_pos_t, hts_pos_t>> insertion_seams;
    std::set<std::pair<hts_pos_t, hts_pos_t>> deletion_snp_seams;
    std::set<std::pair<hts_pos_t, hts_pos_t>> deferred_clean_snp_pairs;
    for (const RecoverySeam& seam : stitch_seams)
        seen_pairs.emplace(seam.left_phase_set, seam.right_phase_set);
    for (const RecoverySeam& seam : collect_phase_set_seams(gc)) {
        const bool inside_target = std::any_of(
            gc.recovery_windows.begin(), gc.recovery_windows.end(),
            [&seam](const RecoverySeam& original) {
                return original.beg <= seam.beg && seam.end <= original.end;
            });
        if (!inside_target ||
            seen_pairs.count({seam.left_phase_set, seam.right_phase_set}) != 0)
            continue;
        bool right_mnp = false;
        bool right_insertion = false;
        bool right_deletion = false;
        bool left_deletion = false;
        for (const RecoverySourceSite& source : gc.recovery_source_sites) {
            if (source.candidate_index >= chunk.candidates.size()) continue;
            const CandidateVariant& candidate =
                chunk.candidates[source.candidate_index];
            if (!candidate.msa_verified) continue;
            if (candidate.phase_set == seam.right_phase_set &&
                candidate.key.type == VariantType::Insertion &&
                candidate.key.alt.size() >
                    kLongPhysicalInsertionThreshold &&
                candidate.key.sort_pos() == seam.end)
                right_insertion = true;
            if (candidate.phase_set == seam.right_phase_set &&
                candidate.key.type == VariantType::Deletion &&
                candidate.key.sort_pos() == seam.end)
                right_deletion = true;
            if (candidate.phase_set == seam.left_phase_set &&
                candidate.key.type == VariantType::Deletion &&
                candidate.key.sort_pos() == seam.beg)
                left_deletion = true;
            if (right_insertion && left_deletion) break;
        }
        for (size_t ci = 0; ci < gc.site_meta.size() &&
                            ci < chunk.candidates.size(); ++ci) {
            if (chunk.candidates[ci].phase_set != seam.right_phase_set ||
                chunk.candidates[ci].bam_injected ||
                gc.site_meta[ci].ref.size() <= 1)
                continue;
            const std::string* alt = selected_graph_candidate_alt(gc, ci);
            const auto pos = alt == nullptr ? std::nullopt :
                graph_substitution_start(gc.site_meta[ci], *alt);
            if (pos && *pos == seam.end) {
                right_mnp = true;
                break;
            }
        }
        // Transfer may expose a right deletion opposite a mixed graph/BAM
        // block. Retry that pair using nearby clean SNPs outside the narrow
        // seam, with the source path and graph edges checked below.
        // A recovered deletion can become the new left boundary after the
        // original graph seam was captured. Revisit its physical allele link
        // to the right SNP under the existing read and block-path checks.
        if ((right_mnp || right_insertion || right_deletion || left_deletion) &&
            seen_pairs.emplace(seam.left_phase_set,
                               seam.right_phase_set).second) {
            stitch_seams.push_back(seam);
            if (right_insertion)
                insertion_seams.emplace(seam.left_phase_set,
                                        seam.right_phase_set);
            if (right_deletion && !right_mnp && !right_insertion &&
                !left_deletion)
                deletion_snp_seams.emplace(seam.left_phase_set,
                                           seam.right_phase_set);
        }
    }
    // A physical bridge can join two original graph blocks after the BAM
    // stitch map was captured. Carry that one proven boundary edge forward so
    // the next seam's SNP-path check does not mistake it for a missing link.
    struct PendingPhysicalEdge {
        size_t left_i;
        size_t right_i;
        hts_pos_t left_ps;
        hts_pos_t right_ps;
    };
    const auto record_physical_edge = [&](const std::optional<PendingPhysicalEdge>& edge) {
        if (!edge || stitched_graph_phase_sets == nullptr) return;
        const CandidateVariant& left = chunk.candidates[edge->left_i];
        const CandidateVariant& right = chunk.candidates[edge->right_i];
        if (edge->left_ps == edge->right_ps || left.phase_set <= 0 ||
            left.phase_set != right.phase_set) return;
        const hts_pos_t joined = left.phase_set;
        (*stitched_graph_phase_sets)[gc.site_ids[edge->left_i]] = joined;
        (*stitched_graph_phase_sets)[gc.site_ids[edge->right_i]] = joined;
    };
    std::optional<PendingPhysicalEdge> pending_edge;
    const auto complete_mixed_left_path = [&](hts_pos_t phase_set,
                                              size_t bridge_i) {
        const auto source = gc.recovery_source_path_supported.find(phase_set);
        if (source == gc.recovery_source_path_supported.end() ||
            !source->second ||
            (gc.recovery_source_weak_cuts.count(phase_set) != 0 &&
             !gc.recovery_source_weak_cuts.at(phase_set).empty()) ||
            (gc.recovery_source_quality_cuts.count(phase_set) != 0 &&
             !gc.recovery_source_quality_cuts.at(phase_set).empty()))
            return false;
        bool shared_bridge = std::any_of(
            gc.recovery_source_sites.begin(), gc.recovery_source_sites.end(),
            [&](const RecoverySourceSite& site) {
                return site.candidate_index == bridge_i &&
                    site.phase_set == phase_set && site.clean_shared_snp &&
                    site.can_adopt &&
                    site.hap1_allele ==
                        chunk.candidates[bridge_i].hap_to_cons_alle[1];
            });
        // A graph SNP can be unphased before transfer yet become the clean
        // boundary of a complete BAM source. Require its exact source allele
        // and a direct, non-reversing edge to the preceding graph SNP.
        if (!shared_bridge && original_graph_phase_sets != nullptr &&
            chunk.candidates[bridge_i].counts.category ==
                VariantCategory::CleanHetSnp &&
            chunk.candidates[bridge_i].alignment_verified) {
            const CandidateVariant& bridge = chunk.candidates[bridge_i];
            const bool exact_source = std::any_of(
                gc.recovery_source_sites.begin(), gc.recovery_source_sites.end(),
                [&](const RecoverySourceSite& site) {
                    return site.candidate_index == bridge_i &&
                        site.phase_set == phase_set && site.can_adopt &&
                        site.graph_phase_set == 0 &&
                        site.hap1_allele == bridge.hap_to_cons_alle[1] &&
                        site.hap2_allele == bridge.hap_to_cons_alle[2];
                });
            std::optional<size_t> previous;
            for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
                const CandidateVariant& candidate = chunk.candidates[ci];
                if (candidate.phase_set != phase_set || candidate.bam_injected ||
                    candidate.counts.category != VariantCategory::CleanHetSnp ||
                    candidate.key.pos >= bridge.key.pos ||
                    (previous && candidate.key.pos <=
                        chunk.candidates[*previous].key.pos))
                    continue;
                previous = ci;
            }
            bool intervening_block = false;
            if (previous) {
                const hts_pos_t start = chunk.candidates[*previous].key.sort_pos();
                const hts_pos_t end = bridge.key.sort_pos();
                for (const CandidateVariant& candidate : chunk.candidates) {
                    const hts_pos_t pos = candidate.key.sort_pos();
                    if (start < pos && pos < end && candidate.phase_set > 0 &&
                        candidate.phase_set != phase_set) {
                        intervening_block = true;
                        break;
                    }
                }
            }
            const auto flip = exact_source && previous && !intervening_block ?
                local_run_boundary_flip(chunk, *previous, bridge_i,
                                        kMinMapq, kMaxWrongParity) :
                std::nullopt;
            shared_bridge = flip && !*flip;
        }
        if (!shared_bridge || original_graph_phase_sets == nullptr)
            return false;
        std::vector<size_t> graph_sites;
        hts_pos_t original_ps = 0;
        for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
            const CandidateVariant& candidate = chunk.candidates[ci];
            if (candidate.phase_set != phase_set || candidate.bam_injected ||
                candidate.hap_to_cons_alle[1] < 0 ||
                candidate.hap_to_cons_alle[2] < 0 ||
                candidate.hap_to_cons_alle[1] == candidate.hap_to_cons_alle[2])
                continue;
            if (ci >= gc.site_ids.size() || gc.site_ids[ci].empty()) {
                return false;
            }
            const auto original =
                original_graph_phase_sets->find(gc.site_ids[ci]);
            if (original == original_graph_phase_sets->end()) {
                // Recovery may have added an unphased graph indel to the
                // complete source. Its source allele must agree with the
                // current row before it can join the local graph path.
                const bool adopted = std::any_of(
                    gc.recovery_source_sites.begin(),
                    gc.recovery_source_sites.end(),
                    [&](const RecoverySourceSite& site) {
                        return site.candidate_index == ci &&
                            site.phase_set == phase_set && site.can_adopt &&
                            site.hap1_allele == candidate.hap_to_cons_alle[1] &&
                            site.hap2_allele == candidate.hap_to_cons_alle[2];
                    });
                if (!adopted) {
                    return false;
                }
            } else {
                if (original->second <= 0 ||
                    (original_ps > 0 && original_ps != original->second)) {
                    return false;
                }
                original_ps = original->second;
            }
            graph_sites.push_back(ci);
        }
        if (graph_sites.size() < 2 ||
            std::find(graph_sites.begin(), graph_sites.end(), bridge_i) ==
                graph_sites.end())
            return false;
        std::sort(graph_sites.begin(), graph_sites.end(),
                  [&](size_t a, size_t b) {
                      return chunk.candidates[a].key.sort_pos() <
                          chunk.candidates[b].key.sort_pos();
                  });
        for (size_t i = 1; i < graph_sites.size(); ++i) {
            if (chunk.candidates[graph_sites[i - 1]].key.sort_pos() >=
                chunk.candidates[graph_sites[i]].key.sort_pos())
                return false;
            const auto flip = local_run_boundary_flip(
                chunk, graph_sites[i - 1], graph_sites[i],
                kMinMapq, kMaxWrongParity);
            if (flip && !*flip) continue;
            if (flip && *flip) return false;
            // A weak direct edge cannot override the source path, but even
            // one callable reversed pair must veto that fallback.
            bool conflicting_pair = false;
            for (size_t ri = 0; ri < chunk.read_var_profile.size() &&
                                ri < chunk.reads.size(); ++ri) {
                const ReadVariantProfile& profile = chunk.read_var_profile[ri];
                if (chunk.reads[ri].is_skipped ||
                    chunk.reads[ri].mapq < kMinMapq ||
                    chunk.reads[ri].mapq == kUnknownQuality ||
                    profile.start_var_idx < 0 ||
                    graph_sites[i - 1] <
                        static_cast<size_t>(profile.start_var_idx) ||
                    graph_sites[i] <
                        static_cast<size_t>(profile.start_var_idx) ||
                    graph_sites[i - 1] >
                        static_cast<size_t>(profile.end_var_idx) ||
                    graph_sites[i] >
                        static_cast<size_t>(profile.end_var_idx))
                    continue;
                const size_t left_offset = graph_sites[i - 1] -
                    static_cast<size_t>(profile.start_var_idx);
                const size_t right_offset = graph_sites[i] -
                    static_cast<size_t>(profile.start_var_idx);
                if (std::max(left_offset, right_offset) >=
                    profile.alleles.size())
                    continue;
                const int left = profile.alleles[left_offset];
                const int right = profile.alleles[right_offset];
                if ((left == 0 || left == 1) && (right == 0 || right == 1) &&
                    (left == chunk.candidates[graph_sites[i - 1]]
                                 .hap_to_cons_alle[1]) !=
                    (right == chunk.candidates[graph_sites[i]]
                                  .hap_to_cons_alle[1])) {
                    conflicting_pair = true;
                    break;
                }
            }
            if (conflicting_pair) return false;
            // A complete source path can bridge an absent or one-sided
            // agreeing edge between exact shared SNPs in the original block.
            const auto shared_source = [&](size_t ci) {
                const CandidateVariant& candidate = chunk.candidates[ci];
                return std::any_of(gc.recovery_source_sites.begin(),
                                   gc.recovery_source_sites.end(),
                                   [&](const RecoverySourceSite& site) {
                    return site.candidate_index == ci &&
                        site.phase_set == phase_set && site.clean_shared_snp &&
                        site.can_adopt &&
                        site.hap1_allele == candidate.hap_to_cons_alle[1] &&
                        site.hap2_allele == candidate.hap_to_cons_alle[2];
                });
            };
            if (!shared_source(graph_sites[i - 1]) ||
                !shared_source(graph_sites[i]))
                return false;
        }
        return true;
    };
    const auto unique_graph_anchor = [&](hts_pos_t pos) -> std::optional<size_t> {
        std::optional<size_t> found;
        for (size_t ci = 0; ci < chunk.candidates.size() &&
                            ci < gc.site_ids.size(); ++ci) {
            if (chunk.candidates[ci].bam_injected || gc.site_ids[ci].empty())
                continue;
            const auto site = physical_substitution(ci);
            if (!site || site->pos != pos) continue;
            if (found) return std::nullopt;
            found = ci;
        }
        return found;
    };
    for (size_t seam_i = 0; seam_i < stitch_seams.size(); ++seam_i) {
        // A join can expose a closer seam. Copy this entry because appending
        // a new seam may grow the vector and invalidate its references.
        const RecoverySeam seam = stitch_seams[seam_i];
        record_physical_edge(pending_edge);
        pending_edge.reset();
        if (original_graph_phase_sets != nullptr &&
            stitched_graph_phase_sets != nullptr) {
            const auto left_i = unique_graph_anchor(seam.beg);
            const auto right_i = unique_graph_anchor(seam.end);
            if (left_i && right_i) {
                const std::string& left_id = gc.site_ids[*left_i];
                const std::string& right_id = gc.site_ids[*right_i];
                const auto left_original =
                    original_graph_phase_sets->find(left_id);
                const auto right_original =
                    original_graph_phase_sets->find(right_id);
                if (left_original != original_graph_phase_sets->end() &&
                    right_original != original_graph_phase_sets->end() &&
                    left_original->second > 0 &&
                    right_original->second > 0 &&
                    left_original->second != right_original->second &&
                    stitched_graph_phase_sets->count(left_id) != 0 &&
                    stitched_graph_phase_sets->count(right_id) != 0) {
                    pending_edge = PendingPhysicalEdge{
                        *left_i, *right_i,
                        chunk.candidates[*left_i].phase_set,
                        chunk.candidates[*right_i].phase_set};
                }
            }
        }
        std::array<std::set<hts_pos_t>, 2> current_ps;
        for (const RecoverySourceSite& site : gc.recovery_source_sites) {
            const size_t side = site.graph_phase_set == seam.left_phase_set ? 0 :
                site.graph_phase_set == seam.right_phase_set ? 1 : 2;
            if (side == 2 || site.candidate_index >= chunk.candidates.size())
                continue;
            const CandidateVariant& candidate =
                chunk.candidates[site.candidate_index];
            if (!candidate.bam_injected && physical_substitution(site.candidate_index))
                current_ps[side].insert(candidate.phase_set);
        }
        // A graph boundary need not also be present in the BAM source list.
        for (size_t side = 0; side < current_ps.size(); ++side) {
            if (!current_ps[side].empty()) continue;
            const hts_pos_t original_ps = side == 0 ?
                seam.left_phase_set : seam.right_phase_set;
            for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
                if (!chunk.candidates[ci].bam_injected &&
                    chunk.candidates[ci].phase_set == original_ps &&
                    physical_substitution(ci)) {
                    current_ps[side].insert(original_ps);
                    break;
                }
            }
        }
        if (std::find_if(gc.recovery_windows.begin(),
                         gc.recovery_windows.end(),
                         [&seam](const RecoverySeam& original) {
                             return original.beg == seam.beg &&
                                 original.end == seam.end &&
                                 original.left_phase_set ==
                                     seam.left_phase_set &&
                                 original.right_phase_set ==
                                     seam.right_phase_set;
                         }) == gc.recovery_windows.end()) {
            current_ps[0] = {seam.left_phase_set};
            current_ps[1] = {seam.right_phase_set};
            if (insertion_seams.count({seam.left_phase_set,
                                       seam.right_phase_set}) != 0) {
                current_ps[0].clear();
                current_ps[1].clear();
                for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
                    const auto site = physical_substitution(ci);
                    if (site && site->pos == seam.beg)
                        current_ps[0].insert(chunk.candidates[ci].phase_set);
                }
                for (const RecoverySourceSite& source : gc.recovery_source_sites) {
                    if (source.candidate_index >= chunk.candidates.size()) continue;
                    const CandidateVariant& candidate =
                        chunk.candidates[source.candidate_index];
                    if (candidate.key.type == VariantType::Insertion &&
                        candidate.key.sort_pos() == seam.end &&
                        candidate.phase_set > 0)
                        current_ps[1].insert(candidate.phase_set);
                }
            }
        }
        if (current_ps[0].size() != 1 || current_ps[1].size() != 1)
            continue;
        const hts_pos_t left_ps = *current_ps[0].begin();
        const hts_pos_t right_ps = *current_ps[1].begin();
        if (left_ps == right_ps) continue;
        const bool deletion_snp_retry = deletion_snp_seams.count(
            {seam.left_phase_set, seam.right_phase_set}) != 0;
        std::optional<size_t> left_i, right_i;
        std::optional<PhysicalSubstitution> left_site, right_site;
        for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
            const auto site = physical_substitution(ci);
            if (!site || site->pos < seam.beg || site->pos > seam.end)
                continue;
            const hts_pos_t ps = chunk.candidates[ci].phase_set;
            if (ps == left_ps &&
                (!left_site || site->pos > left_site->pos)) {
                left_i = ci;
                left_site = site;
            }
            if (ps == right_ps &&
                (!right_site || site->pos < right_site->pos)) {
                right_i = ci;
                right_site = site;
            }
        }
        if (!deletion_snp_retry && !left_i && right_i && right_site) {
            if (stitch_complementary_deletions_to_snp(
                    gc, context, tid, seam, left_ps, right_ps,
                    right_site->pos, right_site->ref, right_site->alt, *right_i,
                    original_graph_phase_sets, stitched_graph_phase_sets) ||
                stitch_msa_indel_to_snp(
                    gc, context, tid, seam, left_ps, right_ps,
                    right_site->pos, right_site->ref, right_site->alt, *right_i,
                    original_graph_phase_sets, stitched_graph_phase_sets))
                continue;
        }
        if (!deletion_snp_retry && left_i && left_site &&
            left_site->ref.size() == 1 && !right_i &&
            (stitch_snp_to_msa_insertion(
                 gc, context, tid, seam, left_ps, right_ps,
                 left_site->pos, left_site->ref[0], left_site->alt[0],
                 *left_i, original_graph_phase_sets,
                 stitched_graph_phase_sets) ||
             stitch_snp_to_msa_deletion(
                 gc, context, tid, seam, left_ps, right_ps,
                 left_site->pos, left_site->ref[0], left_site->alt[0],
                 *left_i, original_graph_phase_sets,
                 stitched_graph_phase_sets)))
            continue;
        if (deletion_snp_retry && !left_i) {
            for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
                const CandidateVariant& candidate = chunk.candidates[ci];
                if (candidate.phase_set != left_ps ||
                    candidate.counts.category != VariantCategory::CleanHetSnp)
                    continue;
                const auto site = physical_substitution(ci);
                if (!site || site->pos >= seam.beg ||
                    site->ref.size() != 1 || site->alt.size() != 1 ||
                    (left_site && site->pos <= left_site->pos))
                    continue;
                left_i = ci;
                left_site = site;
            }
        }
        // Recovery can expose a right insertion just before the first clean
        // SNP in that same block. If the insertion cannot be called reliably,
        // use the next SNP for the physical bridge instead of stopping at the
        // old seam coordinate. Do not jump over another phased block.
        if (left_i && !right_i && left_site &&
            left_site->ref.size() == 1 && left_site->alt.size() == 1 &&
            chunk.candidates[*left_i].counts.category ==
                VariantCategory::CleanHetSnp) {
            for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
                const CandidateVariant& candidate = chunk.candidates[ci];
                if (candidate.phase_set != right_ps ||
                    candidate.counts.category != VariantCategory::CleanHetSnp)
                    continue;
                const auto site = physical_substitution(ci);
                if (!site || site->pos <= seam.end ||
                    site->ref.size() != 1 || site->alt.size() != 1 ||
                    (right_site && site->pos >= right_site->pos))
                    continue;
                right_i = ci;
                right_site = site;
            }
            if (right_site) {
                for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
                    const CandidateVariant& candidate = chunk.candidates[ci];
                    if (candidate.phase_set <= 0 ||
                        candidate.phase_set == left_ps ||
                        candidate.phase_set == right_ps ||
                        (candidate.hap_to_cons_alle[1] == 1) ==
                            (candidate.hap_to_cons_alle[2] == 1))
                        continue;
                    hts_pos_t pos = candidate.key.sort_pos();
                    if (!candidate.bam_injected && ci < gc.site_meta.size()) {
                        const std::string* alt =
                            selected_graph_candidate_alt(gc, ci);
                        if (alt != nullptr && !gc.site_meta[ci].ref.empty())
                            pos = vcf_to_variant_key(
                                candidate.key.tid, gc.site_meta[ci].pos,
                                gc.site_meta[ci].ref, *alt).sort_pos();
                    }
                    if (pos > seam.end && pos < right_site->pos) {
                        right_i.reset();
                        right_site.reset();
                        break;
                    }
                }
            }
        }
        // Let the ordinary right-hand joins settle first. Extending this
        // seam early would change its right block before those path checks.
        if (right_i && right_site->pos > seam.end &&
            deferred_clean_snp_pairs.emplace(seam.left_phase_set,
                                             seam.right_phase_set).second) {
            stitch_seams.push_back(seam);
            continue;
        }
        if (!left_i || !right_i || left_site->pos >= right_site->pos)
            continue;
        bool noisy_bridge = false;
        const size_t clean_right_i = *right_i;
        const PhysicalSubstitution clean_right_site = *right_site;
        if (left_site->ref.size() == 1 && left_site->alt.size() == 1) {
            // A demoted SNP may still have directly callable high-quality BAM
            // bases. Prefer the nearest one when it lies before the clean SNP
            // and test its paired-read parity independently below.
            for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
                const CandidateVariant& candidate = chunk.candidates[ci];
                if (candidate.phase_set != right_ps ||
                    candidate.counts.category != VariantCategory::NoisyCandHet)
                    continue;
                const auto site = physical_substitution(ci, true);
                if (!site || site->ref.size() != 1 || site->alt.size() != 1 ||
                    site->pos <= left_site->pos ||
                    site->pos < seam.beg || site->pos > seam.end ||
                    site->pos >= right_site->pos)
                    continue;
                right_i = ci;
                right_site = site;
                noisy_bridge = true;
            }
        }
        if (deletion_snp_retry && noisy_bridge) continue;
        if (noisy_bridge) {
            // Keep the established clean bridge whenever even one high-quality
            // molecule calls both clean alleles. The noisy candidate is a
            // fallback for missing clean-pair coverage, not a replacement.
            std::unique_ptr<hts_itr_t, decltype(&hts_itr_destroy)> clean_iter(
                sam_itr_queryi(context.indexes.front().get(), tid,
                               left_site->pos - 1, clean_right_site.pos),
                &hts_itr_destroy);
            std::unique_ptr<bam1_t, AlignmentDeleter> clean_read(bam_init1());
            if (!clean_iter || !clean_read) continue;
            while (sam_itr_next(context.bams.front()->get(), clean_iter.get(),
                                clean_read.get()) >= 0) {
                const bam1_t* read = clean_read.get();
                if ((read->core.flag & (BAM_FUNMAP | BAM_FSECONDARY |
                                        BAM_FSUPPLEMENTARY | BAM_FDUP |
                                        BAM_FQCFAIL)) ||
                    read->core.qual < kMinMapq ||
                    read->core.qual == kUnknownQuality)
                    continue;
                int left_quality = 0, right_quality = 0;
                const int left_call = physical_substitution_call(
                    read, left_site->pos, left_site->ref, left_site->alt,
                    &left_quality);
                const int right_call = physical_substitution_call(
                    read, clean_right_site.pos, clean_right_site.ref,
                    clean_right_site.alt, &right_quality);
                if ((left_call != 0 && left_call != 2) ||
                    (right_call != 0 && right_call != 2) ||
                    left_quality < kMinBaseq || right_quality < kMinBaseq ||
                    left_quality == kUnknownQuality ||
                    right_quality == kUnknownQuality)
                    continue;
                right_i = clean_right_i;
                right_site = clean_right_site;
                noisy_bridge = false;
                break;
            }
        }
        // The clean pair had no callable molecule. Preserve the established
        // clean-site indel fallbacks before trying the demoted SNP; those
        // fallbacks must never receive the demoted allele as their anchor.
        if (noisy_bridge &&
            (stitch_graph_snp_to_complementary_deletions(
                 gc, context, tid, seam, left_ps, right_ps,
                 left_site->pos, *left_i, original_graph_phase_sets,
                 stitched_graph_phase_sets) ||
             stitch_msa_indel_to_snp(
                 gc, context, tid, seam, left_ps, right_ps,
                 clean_right_site.pos, clean_right_site.ref,
                 clean_right_site.alt, clean_right_i,
                 original_graph_phase_sets, stitched_graph_phase_sets)))
            continue;
        std::unique_ptr<hts_itr_t, decltype(&hts_itr_destroy)> iterator(
            sam_itr_queryi(context.indexes.front().get(), tid,
                           left_site->pos - 1, right_site->pos),
            &hts_itr_destroy);
        if (!iterator) continue;
        std::unique_ptr<bam1_t, AlignmentDeleter> alignment(bam_init1());
        if (!alignment) continue;
        std::unordered_set<std::string> seen;
        int paired = 0;
        int supporting = 0, opposing = 0;
        std::array<int, 2> left_alleles{};
        double log_odds = 0.0;
        while (sam_itr_next(context.bams.front()->get(), iterator.get(),
                            alignment.get()) >= 0) {
            const bam1_t* read = alignment.get();
            if ((read->core.flag & (BAM_FUNMAP | BAM_FSECONDARY |
                                    BAM_FSUPPLEMENTARY | BAM_FDUP | BAM_FQCFAIL)) ||
                read->core.qual < kMinMapq ||
                read->core.qual == kUnknownQuality)
                continue;
            const char* qname = bam_get_qname(read);
            if (seen.count(qname) != 0) continue;
            int left_quality = 0, right_quality = 0;
            const int left_call = physical_substitution_call(
                read, left_site->pos, left_site->ref, left_site->alt,
                &left_quality);
            const int right_call = physical_substitution_call(
                read, right_site->pos, right_site->ref, right_site->alt,
                &right_quality);
            if ((left_call != 0 && left_call != 2) ||
                (right_call != 0 && right_call != 2) ||
                left_quality < kMinBaseq || right_quality < kMinBaseq ||
                left_quality == kUnknownQuality ||
                right_quality == kUnknownQuality)
                continue;
            seen.insert(qname);
            const double p =
                std::pow(10.0, -static_cast<double>(left_quality) / 10.0) +
                std::pow(10.0, -static_cast<double>(right_quality) / 10.0) +
                2.0 * std::pow(10.0,
                    -static_cast<double>(read->core.qual) / 10.0);
            if (p <= 0.0 || p >= 0.5) continue;
            const bool left_hap1 = (left_call == 2) ==
                (chunk.candidates[*left_i].hap_to_cons_alle[1] == 1);
            const bool right_hap1 = (right_call == 2) ==
                (chunk.candidates[*right_i].hap_to_cons_alle[1] == 1);
            const double weight = std::log((1.0 - p) / p);
            const bool flip = left_hap1 != right_hap1;
            log_odds += flip ? weight : -weight;
            (flip ? supporting : opposing)++;
            ++left_alleles[left_call == 2 ? 1 : 0];
            ++paired;
        }
        const double threshold = std::log((1.0 - kMaxWrongParity) /
                                           kMaxWrongParity);
        // An exposed deletion seam may have only one callable clean-SNP
        // allele. Try the exact complementary ALT bridge before rejecting
        // that retry; any available SNP votes remain an orientation veto.
        if (!noisy_bridge &&
            (paired == 0 || (deletion_snp_retry &&
                (left_alleles[0] == 0 || left_alleles[1] == 0))) &&
            (supporting == 0 || opposing == 0) && left_i && left_site &&
            left_site->ref.size() == 1 &&
            stitch_graph_snp_to_complementary_deletions(
                gc, context, tid, seam, left_ps, right_ps,
                left_site->pos, *left_i, original_graph_phase_sets,
                stitched_graph_phase_sets, paired == 0 ? std::nullopt :
                    std::optional<bool>(supporting != 0)))
            continue;
        if (deletion_snp_retry &&
            (paired == 0 || left_alleles[0] == 0 || left_alleles[1] == 0 ||
             (supporting != 0 && opposing != 0)))
            continue;
        // A nearer right deletion may be callable when the clean SNP pair
        // is not. Its original-source path and current gauge are checked by
        // the helper before it can orient the entire right block.
        if (!noisy_bridge && paired == 0 && left_i && left_site &&
            left_site->ref.size() == 1 &&
            stitch_snp_to_msa_deletion(
                gc, context, tid, seam, left_ps, right_ps,
                left_site->pos, left_site->ref[0], left_site->alt[0],
                *left_i, original_graph_phase_sets,
                stitched_graph_phase_sets))
            continue;
        // Prefer the clean-SNP pair. Without callable pairs, a nearer left
        // MSA locus may bridge the gap even though the left SNP exists. Try
        // separate complementary rows before the single-indel helper, which
        // rejects overlapping alleles. Both retain their full path checks.
        if (!noisy_bridge && paired == 0 && right_i && right_site &&
            (stitch_complementary_deletions_to_snp(
                 gc, context, tid, seam, left_ps, right_ps,
                 right_site->pos, right_site->ref, right_site->alt, *right_i,
                 original_graph_phase_sets, stitched_graph_phase_sets) ||
             stitch_msa_indel_to_snp(
                 gc, context, tid, seam, left_ps, right_ps,
                 right_site->pos, right_site->ref, right_site->alt, *right_i,
                 original_graph_phase_sets, stitched_graph_phase_sets)))
            continue;
        // Independent reads from one haplotype can establish the diploid
        // parity; requiring observations from both haplotypes discards valid
        // low-coverage bridges despite a decisive likelihood.
        if (noisy_bridge &&
            ((supporting != 0 && opposing != 0) ||
             left_alleles[0] == 0 || left_alleles[1] == 0 ||
             std::ldexp(1.0, -paired) > kMaxWrongParity))
            continue;
        if (noisy_bridge) {
            // The demoted SNP must also agree with the first clean SNP in
            // its own block. A strong cross-seam vote can otherwise invert a
            // large block when this one noisy row has the wrong orientation.
            std::optional<size_t> clean_i;
            std::optional<PhysicalSubstitution> clean_site;
            for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
                if (chunk.candidates[ci].phase_set != right_ps) continue;
                const auto site = physical_substitution(ci);
                if (!site || site->ref.size() != 1 ||
                    site->alt.size() != 1 ||
                    site->pos <= right_site->pos ||
                    (clean_site && site->pos >= clean_site->pos))
                    continue;
                clean_i = ci;
                clean_site = site;
            }
            if (!clean_i) continue;
            std::unique_ptr<hts_itr_t, decltype(&hts_itr_destroy)> check_iter(
                sam_itr_queryi(context.indexes.front().get(), tid,
                               right_site->pos - 1, clean_site->pos),
                &hts_itr_destroy);
            if (!check_iter) continue;
            std::unordered_set<std::string> checked_reads;
            int same = 0, cross = 0;
            std::array<int, 2> noisy_alleles{};
            while (sam_itr_next(context.bams.front()->get(), check_iter.get(),
                                alignment.get()) >= 0) {
                const bam1_t* read = alignment.get();
                if ((read->core.flag & (BAM_FUNMAP | BAM_FSECONDARY |
                                        BAM_FSUPPLEMENTARY | BAM_FDUP |
                                        BAM_FQCFAIL)) ||
                    read->core.qual < kMinMapq ||
                    read->core.qual == kUnknownQuality)
                    continue;
                const char* qname = bam_get_qname(read);
                if (checked_reads.count(qname) != 0) continue;
                int noisy_quality = 0, clean_quality = 0;
                const int noisy_call = physical_substitution_call(
                    read, right_site->pos, right_site->ref, right_site->alt,
                    &noisy_quality);
                const int clean_call = physical_substitution_call(
                    read, clean_site->pos, clean_site->ref, clean_site->alt,
                    &clean_quality);
                if ((noisy_call != 0 && noisy_call != 2) ||
                    (clean_call != 0 && clean_call != 2) ||
                    noisy_quality < kMinBaseq ||
                    clean_quality < kMinBaseq ||
                    noisy_quality == kUnknownQuality ||
                    clean_quality == kUnknownQuality)
                    continue;
                checked_reads.insert(qname);
                const bool noisy_hap1 = (noisy_call == 2) ==
                    (chunk.candidates[*right_i].hap_to_cons_alle[1] == 1);
                const bool clean_hap1 = (clean_call == 2) ==
                    (chunk.candidates[*clean_i].hap_to_cons_alle[1] == 1);
                ++(noisy_hap1 == clean_hap1 ? same : cross);
                ++noisy_alleles[noisy_call == 2 ? 1 : 0];
            }
            if (noisy_alleles[0] == 0 || noisy_alleles[1] == 0 ||
                !source_graph_vote_supported(same, cross, 0,
                                             kMaxWrongParity))
                continue;
        }
        if (paired == 0 || std::abs(log_odds) < threshold) continue;
        std::pair<hts_pos_t, hts_pos_t> left_cut;
        bool left_path = graph_snp_path_supported(
            gc, left_ps, &left_cut, original_graph_phase_sets,
            stitched_graph_phase_sets, false, false, 0, true) ||
            bam_source_run_supported(gc, left_ps) ||
            ((noisy_bridge || deletion_snp_retry) &&
             complete_mixed_left_path(left_ps, *left_i));
        // A newly merged block may contain a one-haplotype GAF edge. Direct
        // Q30 BAM pairs can corroborate that edge without relaxing the path
        // check elsewhere: the prefix already passed and the helper checks
        // the full suffix. A dominant GAF reversal never supplies a cut here.
        if (!left_path && left_cut.first > 0) {
            left_path = physical_graph_snp_edge_supported(
                gc, context, tid, left_ps, left_cut,
                original_graph_phase_sets, stitched_graph_phase_sets, true)
                    .has_value();
        }
        if (!left_path) continue;
        std::pair<hts_pos_t, hts_pos_t> disconnected;
        const auto source_path = gc.recovery_source_path_supported.find(right_ps);
        const auto weak_cuts = gc.recovery_source_weak_cuts.find(right_ps);
        const auto quality_cuts = gc.recovery_source_quality_cuts.find(right_ps);
        const bool complete_right_source =
            source_path != gc.recovery_source_path_supported.end() &&
            source_path->second &&
            (weak_cuts == gc.recovery_source_weak_cuts.end() ||
             weak_cuts->second.empty()) &&
            (quality_cuts == gc.recovery_source_quality_cuts.end() ||
             quality_cuts->second.empty());
        const bool right_graph_path = graph_snp_path_supported(
            gc, right_ps, &disconnected, original_graph_phase_sets,
            stitched_graph_phase_sets);
        if (right_graph_path ||
            (complete_right_source && graph_snp_path_supported(
                gc, right_ps, nullptr, original_graph_phase_sets,
                stitched_graph_phase_sets, false, false, 0, false, true,
                noisy_bridge))) {
            merge_phase_sets_in_place(chunk, left_ps, right_ps,
                                      log_odds > 0.0);
            // The previous boundary may have hidden a nearer BAM source.
            // Retry only a clean-SNP left boundary inside the target: an
            // exposed indel-to-singleton edge can absorb a mixed block.
            for (const RecoverySeam& next : collect_phase_set_seams(gc)) {
                const bool inside_target = std::any_of(
                    gc.recovery_windows.begin(), gc.recovery_windows.end(),
                    [&next](const RecoverySeam& target) {
                        return target.beg <= next.beg &&
                            next.end <= target.end;
                    });
                bool clean_snp_left = false;
                for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
                    if (chunk.candidates[ci].phase_set != next.left_phase_set)
                        continue;
                    const auto site = physical_substitution(ci);
                    if (site && site->pos == next.beg &&
                        site->ref.size() == 1 && site->alt.size() == 1) {
                        clean_snp_left = true;
                        break;
                    }
                }
                if (inside_target && clean_snp_left &&
                    seen_pairs.emplace(next.left_phase_set,
                                       next.right_phase_set).second)
                    stitch_seams.push_back(next);
            }
        } else if (disconnected.first > 0) {
            join_graph_prefix_at_disconnected_cut(
                gc, context, tid, left_ps, right_ps, log_odds > 0.0,
                right_site->pos, disconnected);
        }
    }
    record_physical_edge(pending_edge);
    // A BAM-only insertion pair has no graph substitution to populate
    // current_ps above. Revisit adjacent phase-set seams directly after
    // ordinary physical SNP stitches have settled their current labels.
    for (const RecoverySeam& seam : collect_phase_set_seams(gc)) {
        const bool inside_target = std::any_of(
            gc.recovery_windows.begin(), gc.recovery_windows.end(),
            [&seam](const RecoverySeam& target) {
                return target.beg <= seam.beg && seam.end <= target.end;
            });
        if (inside_target &&
            stitch_repeat_insertion_pair(gc, context, tid, seam,
                                         seam.left_phase_set,
                                         seam.right_phase_set))
            continue;
        if (stitch_deletion_to_phased_right_reads(gc, context, tid, seam)) {
            // Moving the deletion into the right block exposes its previous
            // left neighbor as a new seam; the initial seam list is a snapshot.
            for (const RecoverySeam& exposed : collect_phase_set_seams(gc)) {
                if (exposed.end != seam.beg) continue;
                const bool targeted = std::any_of(
                    gc.recovery_windows.begin(), gc.recovery_windows.end(),
                    [&exposed](const RecoverySeam& target) {
                        return target.beg <= exposed.beg &&
                               exposed.end <= target.end;
                    });
                if (targeted)
                    stitch_complementary_left_deletions(
                        gc, context, tid, exposed,
                        original_graph_phase_sets,
                        stitched_graph_phase_sets);
            }
            continue;
        }
        if (inside_target && stitch_complementary_left_deletions(
                gc, context, tid, seam, original_graph_phase_sets,
                stitched_graph_phase_sets))
            continue;
        if (inside_target && stitch_bam_insertion_run_to_left_reads(
                gc, seam, original_graph_phase_sets,
                stitched_graph_phase_sets))
            continue;
        if (inside_target && stitch_orphan_snp_block(seam))
            continue;
        // Catalog deletion pairs already have graph observations and need no
        // targeted BAM solve before direct primary-read support can join them.
        stitch_graph_deletion_pair(gc, context, tid, seam,
                                   original_graph_phase_sets,
                                   stitched_graph_phase_sets);
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

struct SavedWeakCutSite {
    VariantKey key;
    int hap1 = -1;
    int hap2 = -1;
};

struct SavedWeakCutComponent {
    hts_pos_t source_phase_set = 0;
    SavedWeakCutSite graph_snp;
    std::vector<SavedWeakCutSite> private_sites;
};

// The first pass may have already projected a terminal graph SNP into one
// supported run of a weak-cut BAM source. Save only blocks with that single
// graph SNP and private BAM sites in its local source component.
static std::vector<SavedWeakCutComponent> save_weak_cut_components(
        const GraphChunkBuildResult& gc) {
    const PhasingChunk& chunk = gc.chunk;
    const auto oriented = [](const CandidateVariant& candidate) {
        return candidate.phase_set > 0 &&
            candidate.hap_to_cons_alle[1] >= 0 &&
            candidate.hap_to_cons_alle[2] >= 0 &&
            candidate.hap_to_cons_alle[1] != candidate.hap_to_cons_alle[2];
    };
    std::vector<SavedWeakCutComponent> saved;
    for (const auto& [source_ps, cuts] : gc.recovery_source_weak_cuts) {
        if (cuts.empty()) continue;
        const CandidateVariant* graph_snp = nullptr;
        bool valid = true;
        for (const CandidateVariant& candidate : chunk.candidates) {
            if (candidate.phase_set != source_ps) continue;
            if (!oriented(candidate)) continue;
            if (!candidate.bam_injected) {
                if (graph_snp != nullptr ||
                    candidate.key.type != VariantType::Snp) {
                    valid = false;
                    break;
                }
                graph_snp = &candidate;
            }
        }
        if (!valid || graph_snp == nullptr) continue;
        const size_t graph_component = static_cast<size_t>(std::lower_bound(
            cuts.begin(), cuts.end(), graph_snp->key.sort_pos()) - cuts.begin());
        SavedWeakCutComponent component;
        component.source_phase_set = source_ps;
        component.graph_snp = {graph_snp->key,
            graph_snp->hap_to_cons_alle[1], graph_snp->hap_to_cons_alle[2]};
        for (const CandidateVariant& candidate : chunk.candidates) {
            if (candidate.phase_set != source_ps || !candidate.bam_injected ||
                !oriented(candidate)) continue;
            const size_t run = static_cast<size_t>(std::lower_bound(
                cuts.begin(), cuts.end(), candidate.key.sort_pos()) - cuts.begin());
            if (run == graph_component)
                component.private_sites.push_back({candidate.key,
                    candidate.hap_to_cons_alle[1],
                    candidate.hap_to_cons_alle[2]});
        }
        if (!component.private_sites.empty()) saved.push_back(std::move(component));
    }
    return saved;
}

// An upstream retry must not tear apart a sequence-validated graph SNP and
// private sites in its original weak-cut-free source run. Move only that run
// into the newly supported upstream PS; sites beyond the cut stay independent.
static void restore_weak_cut_components(
        GraphChunkBuildResult& gc,
        const std::vector<SavedWeakCutComponent>& saved) {
    PhasingChunk& chunk = gc.chunk;
    const auto same_key = [](const VariantKey& a, const VariantKey& b) {
        return a.tid == b.tid && a.pos == b.pos && a.type == b.type &&
               a.ref_len == b.ref_len && a.alt == b.alt;
    };
    const auto find_site = [&](const VariantKey& key) -> std::optional<size_t> {
        std::optional<size_t> found;
        for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
            if (!same_key(chunk.candidates[ci].key, key)) continue;
            if (found) return std::nullopt;
            found = ci;
        }
        return found;
    };
    const auto orientation_flip = [](const CandidateVariant& candidate,
                                     const SavedWeakCutSite& before)
            -> std::optional<bool> {
        if (candidate.hap_to_cons_alle[1] == before.hap1 &&
            candidate.hap_to_cons_alle[2] == before.hap2) return false;
        if (candidate.hap_to_cons_alle[1] == before.hap2 &&
            candidate.hap_to_cons_alle[2] == before.hap1) return true;
        return std::nullopt;
    };
    const auto flip_hap = [](int hap) {
        return hap == 1 ? 2 : (hap == 2 ? 1 : hap);
    };
    for (const SavedWeakCutComponent& component : saved) {
        const auto graph_i = find_site(component.graph_snp.key);
        if (!graph_i) continue;
        const CandidateVariant& graph_snp = chunk.candidates[*graph_i];
        const hts_pos_t upstream_ps = graph_snp.phase_set;
        if (upstream_ps <= 0 || upstream_ps == component.source_phase_set)
            continue;
        const auto graph_flip = orientation_flip(graph_snp, component.graph_snp);
        if (!graph_flip) continue;
        std::vector<size_t> private_indices;
        std::optional<bool> move_flip;
        bool valid = true;
        for (const SavedWeakCutSite& before : component.private_sites) {
            const auto ci = find_site(before.key);
            if (!ci || chunk.candidates[*ci].phase_set !=
                           component.source_phase_set) {
                valid = false;
                break;
            }
            const auto private_flip =
                orientation_flip(chunk.candidates[*ci], before);
            if (!private_flip) { valid = false; break; }
            const bool flip = *graph_flip != *private_flip;
            if (move_flip && *move_flip != flip) { valid = false; break; }
            move_flip = flip;
            private_indices.push_back(*ci);
        }
        if (!valid || private_indices.empty() || !move_flip) continue;
        std::vector<bool> moved(chunk.candidates.size(), false);
        for (const size_t ci : private_indices) moved[ci] = true;
        for (const size_t ci : private_indices) {
            CandidateVariant& candidate = chunk.candidates[ci];
            if (*move_flip) {
                std::swap(candidate.hap_to_cons_alle[1],
                          candidate.hap_to_cons_alle[2]);
                std::swap(candidate.hap_to_alle_profile[1],
                          candidate.hap_to_alle_profile[2]);
                candidate.hap_alt = flip_hap(candidate.hap_alt);
                candidate.hap_ref = flip_hap(candidate.hap_ref);
            }
            candidate.phase_set = upstream_ps;
        }
        for (size_t ri = 0; ri < chunk.reads.size(); ++ri) {
            if (ri >= chunk.phase_sets.size() || ri >= chunk.haps.size() ||
                ri >= chunk.read_var_profile.size() ||
                chunk.phase_sets[ri] != component.source_phase_set)
                continue;
            const ReadVariantProfile& profile = chunk.read_var_profile[ri];
            if (profile.start_var_idx < 0) continue;
            bool sees_moved = false, sees_remaining = false;
            for (size_t offset = 0; offset < profile.alleles.size(); ++offset) {
                if (profile.alleles[offset] < 0) continue;
                const size_t ci = static_cast<size_t>(profile.start_var_idx) + offset;
                if (ci >= chunk.candidates.size()) break;
                if (moved[ci]) {
                    const CandidateVariant& candidate = chunk.candidates[ci];
                    const int allele = profile.alleles[offset];
                    if (allele == candidate.hap_to_cons_alle[1] ||
                        allele == candidate.hap_to_cons_alle[2])
                        sees_moved = true;
                } else if (chunk.candidates[ci].phase_set ==
                           component.source_phase_set) {
                    sees_remaining = true;
                }
            }
            if (!sees_moved || sees_remaining) continue;
            chunk.phase_sets[ri] = upstream_ps;
            if (*move_flip) chunk.haps[ri] = flip_hap(chunk.haps[ri]);
        }
    }
}

// A newly exposed seam can add alignment-verified SNPs after graph read labels
// were assigned. Refresh a read only when it calls one of those sites and all
// informative SNPs in its current block agree on the same haplotype.
static void refresh_new_seam_read_haps(GraphChunkBuildResult& gc) {
    constexpr int kMinReadMapq = 30;
    constexpr int kMinBaseQuality = 20;
    PhasingChunk& chunk = gc.chunk;
    const auto eligible_snp = [](const CandidateVariant& candidate) {
        return candidate.key.type == VariantType::Snp &&
            (candidate.counts.category == VariantCategory::CleanHetSnp ||
             (candidate.counts.category == VariantCategory::NoisyCandHet &&
              candidate.msa_verified && candidate.alignment_verified));
    };
    std::vector<bool> new_snp(chunk.candidates.size(), false);
    for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
        const CandidateVariant& candidate = chunk.candidates[ci];
        if (!candidate.bam_injected || !candidate.alignment_verified ||
            !eligible_snp(candidate) || candidate.phase_set <= 0 ||
            candidate.hap_to_cons_alle[1] < 0 ||
            candidate.hap_to_cons_alle[2] < 0 ||
            candidate.hap_to_cons_alle[1] == candidate.hap_to_cons_alle[2])
            continue;
        const hts_pos_t pos = candidate.key.sort_pos();
        new_snp[ci] = std::any_of(gc.recovery_windows.begin(),
            gc.recovery_windows.end(), [pos](const RecoverySeam& seam) {
                return seam.beg < pos && pos < seam.end;
            });
    }
    for (size_t ri = 0; ri < chunk.reads.size(); ++ri) {
        if (ri >= chunk.read_var_profile.size() ||
            ri >= chunk.phase_sets.size() || ri >= chunk.haps.size() ||
            chunk.reads[ri].is_skipped)
            continue;
        const ReadVariantProfile& profile = chunk.read_var_profile[ri];
        if (profile.start_var_idx < 0 ||
            profile.bam_mapq < kMinReadMapq || profile.bam_mapq == 255)
            continue;
        if (chunk.phase_sets[ri] <= 0 || chunk.haps[ri] <= 0) {
            hts_pos_t new_phase_set = 0;
            int new_hap = 0;
            bool conflict = false;
            for (int vi = profile.start_var_idx; vi <= profile.end_var_idx; ++vi) {
                if (vi < 0 || static_cast<size_t>(vi) >= new_snp.size() ||
                    !new_snp[vi]) continue;
                const CandidateVariant& candidate = chunk.candidates[vi];
                if (chunk.phase_sets[ri] > 0 &&
                    candidate.phase_set != chunk.phase_sets[ri]) continue;
                const size_t offset =
                    static_cast<size_t>(vi - profile.start_var_idx);
                if (offset >= profile.bam_alleles.size() ||
                    offset >= profile.bam_base_qualities.size() ||
                    profile.bam_base_qualities[offset] < kMinBaseQuality)
                    continue;
                const int allele = profile.bam_alleles[offset];
                const int hap = allele == candidate.hap_to_cons_alle[1] ? 1 :
                    allele == candidate.hap_to_cons_alle[2] ? 2 : 0;
                if (hap == 0) continue;
                if ((new_phase_set > 0 && new_phase_set != candidate.phase_set) ||
                    (new_hap > 0 && new_hap != hap)) {
                    conflict = true;
                    break;
                }
                new_phase_set = candidate.phase_set;
                new_hap = hap;
            }
            if (!conflict && new_phase_set > 0 && new_hap > 0) {
                chunk.phase_sets[ri] = new_phase_set;
                chunk.haps[ri] = new_hap;
            }
            continue;
        }
        int votes[3] = {0, 0, 0};
        bool saw_new_snp = false;
        for (int vi = profile.start_var_idx; vi <= profile.end_var_idx; ++vi) {
            if (vi < 0 || static_cast<size_t>(vi) >= chunk.candidates.size())
                continue;
            const CandidateVariant& candidate = chunk.candidates[vi];
            if (candidate.phase_set != chunk.phase_sets[ri] ||
                !eligible_snp(candidate) ||
                candidate.hap_to_cons_alle[1] < 0 ||
                candidate.hap_to_cons_alle[2] < 0 ||
                candidate.hap_to_cons_alle[1] ==
                    candidate.hap_to_cons_alle[2])
                continue;
            const size_t offset = static_cast<size_t>(vi - profile.start_var_idx);
            if (offset >= profile.alleles.size()) continue;
            int allele = profile.alleles[offset];
            if (new_snp[vi]) {
                if (offset >= profile.bam_alleles.size() ||
                    offset >= profile.bam_base_qualities.size() ||
                    profile.bam_base_qualities[offset] < kMinBaseQuality)
                    continue;
                allele = profile.bam_alleles[offset];
            }
            if (allele == candidate.hap_to_cons_alle[1]) {
                ++votes[1];
                if (new_snp[vi]) saw_new_snp = true;
            } else if (allele == candidate.hap_to_cons_alle[2]) {
                ++votes[2];
                if (new_snp[vi]) saw_new_snp = true;
            }
        }
        if (saw_new_snp && (votes[1] == 0) != (votes[2] == 0))
            chunk.haps[ri] = votes[1] > 0 ? 1 : 2;
    }
}

// In-chunk recovery for one chunk.
///
/// The BAM sub-solve phases injected sites and reads in an independent local
/// gauge. Once merged, decisive allele edges stitch left graph PS -> local PS
/// -> right graph PS. Existing graph blocks are never globally re-clustered.
static bool run_in_chunk_recovery_once(
        GraphChunkBuildResult& gc, const Options& opts, WorkerContext& ctx,
        const char* contig,
        const std::map<std::string, hts_pos_t>& original_graph_phase_sets,
        std::map<std::string, hts_pos_t>& stitched_graph_phase_sets,
        const std::vector<RecoverySeam>* completed_seams = nullptr) {
    if (!recover_phase_set_seams_in_place(gc, opts, ctx, contig,
                                          completed_seams)) return false;

    Options stitch_opts = opts;
    constexpr int kRecoveryLinkWindow = 128;
    stitch_opts.block_link_window = kRecoveryLinkWindow;
    // Admit one-read edge nominations; the recovery stitch applies its
    // independent block and allele checks before joining phase sets.
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
    if (stitched_graph_phase_sets.empty()) {
        for (size_t ci = 0; ci < gc.chunk.candidates.size() &&
                            ci < gc.site_ids.size(); ++ci)
            if (!gc.site_ids[ci].empty() &&
                gc.chunk.candidates[ci].phase_set > 0)
                stitched_graph_phase_sets.emplace(
                    gc.site_ids[ci], gc.chunk.candidates[ci].phase_set);
    }
    detach_bam_sites_across_weak_cuts(gc);
    // These source blocks already oriented both graph flanks through a local
    // validated path. A later one-sided run attachment would split that join.
    attach_bam_source_phase_sets(gc, locally_bridged_sources);
    attach_unanchored_bam_source_runs(gc, locally_bridged_sources);
    attach_singleton_bam_sources(gc);
    stitch_complete_bam_sources_at_repeat_cut(gc);
    stitch_bam_snp_runs_to_graph(gc);
    const int tid = contig == nullptr ? -1 :
        sam_hdr_name2tid(ctx.primary_header(), contig);
    if (tid >= 0 && !ctx.bams.empty() && !ctx.indexes.empty())
        attach_readless_insertion_source_blocks(
            gc, ctx, tid, &original_graph_phase_sets,
            &stitched_graph_phase_sets);
    stitch_physical_allele_seams(gc, ctx, contig,
                              &original_graph_phase_sets,
                              &stitched_graph_phase_sets);
    dump_recovery_phase_state(gc.chunk, stitch_opts, "recovery-final");
    return true;
}

static void run_in_chunk_recovery(GraphChunkBuildResult& gc,
                                  const Options& opts,
                                  WorkerContext& ctx,
                                  const char* contig) {
    // Candidate indices can change when BAM sites are injected. Graph site
    // IDs retain the original atomic block identity across both recovery
    // passes, so an already certified join can fill a GAF-only path gap.
    std::map<std::string, hts_pos_t> original_graph_phase_sets;
    for (size_t ci = 0; ci < gc.chunk.candidates.size() &&
                        ci < gc.site_ids.size(); ++ci)
        if (!gc.site_ids[ci].empty() &&
            gc.chunk.candidates[ci].phase_set > 0)
            original_graph_phase_sets.emplace(
                gc.site_ids[ci], gc.chunk.candidates[ci].phase_set);
    std::map<std::string, hts_pos_t> stitched_graph_phase_sets;
    if (!run_in_chunk_recovery_once(
            gc, opts, ctx, contig, original_graph_phase_sets,
            stitched_graph_phase_sets)) return;
    // A source attachment may split a graph block after seam discovery. Retry
    // only newly exposed, nonoverlapping intervals once; the previous source
    // decisions are already materialized in the chunk's candidate/read state.
    std::vector<RecoveryPhaseGauge> complete_gauges = gc.recovery_phase_gauges;
    auto complete_source_paths = gc.recovery_source_path_supported;
    const std::vector<RecoverySeam> completed = gc.recovery_windows;
    const auto weak_components = save_weak_cut_components(gc);
    if (run_in_chunk_recovery_once(
            gc, opts, ctx, contig, original_graph_phase_sets,
            stitched_graph_phase_sets, &completed)) {
        restore_weak_cut_components(gc, weak_components);
        refresh_new_seam_read_haps(gc);
        complete_gauges.insert(complete_gauges.begin(),
            gc.recovery_phase_gauges.begin(), gc.recovery_phase_gauges.end());
        for (const auto& path : gc.recovery_source_path_supported) {
            const auto inserted = complete_source_paths.emplace(path);
            // PS allocation is local to each solve. An older absorbed label
            // can be reused: preserve its evidence, but do not conflate the
            // two source identities into one reusable block certificate.
            if (!inserted.second) inserted.first->second = false;
        }
    }
    // Source attachment and retries may change both PS ownership and gauge.
    // Complete-block unions must observe that finalized state, and must not
    // change which source blocks the earlier transfer stages adopt.
    stitch_complete_recovery_phase_blocks(
        gc.chunk, collect_phase_set_seams(gc), complete_gauges, opts,
        complete_source_paths, complete_graph_source_paths(gc));
    // Resolve repeat-placement seams after both solves have finalized their
    // candidate indices and gauges. Core-only proof does not certify the
    // output-only read rescues, so defer the actual merges until those finish.
    const int tid = contig == nullptr ? -1 :
        sam_hdr_name2tid(ctx.primary_header(), contig);
    if (tid < 0) return;
    for (const RecoverySeam& seam : collect_phase_set_seams(gc)) {
        const bool targeted = std::any_of(
            gc.recovery_windows.begin(), gc.recovery_windows.end(),
            [&seam](const RecoverySeam& window) {
                return window.beg <= seam.beg && seam.end <= window.end;
            });
        if (!targeted) continue;
        const auto join = find_equivalent_source_insertion_join(
            gc, ctx, tid, seam, &original_graph_phase_sets,
            &stitched_graph_phase_sets);
        if (join) gc.equivalent_insertion_joins.push_back(*join);
    }
}

static void apply_equivalent_insertion_joins(
        std::vector<GraphChunkBuildResult>& graph_chunks) {
    for (GraphChunkBuildResult& gc : graph_chunks) {
        for (const auto& [left_i, right_i] : gc.equivalent_insertion_joins) {
            const CandidateVariant& left = gc.chunk.candidates[left_i];
            const CandidateVariant& right = gc.chunk.candidates[right_i];
            // Cross-chunk stitching may already have joined or flipped these
            // blocks. Recompute from current allele orientations and update
            // every chunk carrying the downstream PS, including its tail.
            if (!is_phase_set_anchor(left) || !is_phase_set_anchor(right)) continue;
            const hts_pos_t left_ps = left.phase_set;
            const hts_pos_t right_ps = right.phase_set;
            const bool flip = left.hap_to_cons_alle[1] != right.hap_to_cons_alle[1];
            for (GraphChunkBuildResult& block : graph_chunks)
                merge_phase_sets_in_place(block.chunk, left_ps, right_ps, flip);
        }
    }
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
                // Output reconstruction must retain the recovery proof used by
                // duplicate selection. Otherwise a verified noisy heterozygote
                // loses to a deeper unphased graph repeat at the same allele.
                cand.bam_injected = mcand.bam_injected;
                cand.msa_verified = mcand.msa_verified;
                cand.alignment_verified = mcand.alignment_verified;

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
                // Recovery already classified an MSA heterozygote from its
                // own reads and haplotype consensus. Its REF/ALT depth can be
                // very uneven; reapplying graph discovery thresholds here
                // turns a valid phased source row into LOW_AF at emission.
                const bool verified_recovery_het = mcand.bam_injected &&
                    mcand.counts.category == VariantCategory::NoisyCandHet &&
                    mcand.msa_verified && mcand.alignment_verified &&
                    mcand.phase_set > 0 && het_by_consensus;
                if (verified_recovery_het) {
                    cand.counts.category = mcand.counts.category;
                    cand.counts.candvarcate_initial =
                        mcand.counts.candvarcate_initial;
                    cand.lcd_var_i_to_cate = mcand.lcd_var_i_to_cate;
                } else if (alt_cov < opts.min_alt_depth || total_cov < opts.min_depth) {
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

    // The same physical variant can come from a catalog snarl and a verified
    // BAM recovery row. Retain its phased heterozygote when the other copy is
    // an unphased repeat call; coverage alone would discard the recovered
    // phase whenever the graph row had a few more observations.
    const auto phased_het = [](const CandidateVariant& candidate) {
        const auto category = candidate.counts.category;
        return candidate.phase_set > 0 &&
            candidate.hap_to_cons_alle[1] >= 0 &&
            candidate.hap_to_cons_alle[2] >= 0 &&
            candidate.hap_to_cons_alle[1] != candidate.hap_to_cons_alle[2] &&
            (category == VariantCategory::CleanHetSnp ||
             category == VariantCategory::CleanHetIndel ||
             (category == VariantCategory::NoisyCandHet &&
              candidate.msa_verified && candidate.alignment_verified));
    };
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
                const bool write_phased = phased_het(result[write]);
                const bool read_phased = phased_het(result[read]);
                if ((read_phased && !write_phased) ||
                    (read_phased == write_phased &&
                     result[read].counts.total_cov >
                         result[write].counts.total_cov)) {
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
                    if (thread_recovery_ctx != nullptr)
                        supplement_phased_snp_branches(
                            chunk_view, chunk_rows, graph_chunks[offset], opts);

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
        if (!graph_chunks[i].recovery_windows.empty())
            refresh_recovered_read_haps_from_bam_snps(graph_chunks[i].chunk);
        rescue_unphased_graph_reads(
            graph_chunks[i].chunk, graph_chunks[i].recovery_windows);
        apply_independent_bam_read_blocks(graph_chunks[i].chunk);
    }
    apply_equivalent_insertion_joins(graph_chunks);

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
                    if (thread_recovery_ctx != nullptr)
                        supplement_phased_snp_branches(
                            chunk_view, chunk_rows, graph_chunks[offset], opts);

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
        if (!graph_chunks[i].recovery_windows.empty())
            refresh_recovered_read_haps_from_bam_snps(graph_chunks[i].chunk);
        rescue_unphased_graph_reads(
            graph_chunks[i].chunk, graph_chunks[i].recovery_windows);
        apply_independent_bam_read_blocks(graph_chunks[i].chunk);
    }
    apply_equivalent_insertion_joins(graph_chunks);

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
