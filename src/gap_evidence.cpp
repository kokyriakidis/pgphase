#include "gap_evidence.hpp"

#include "collect_phase.hpp"
#include "collect_var.hpp"

#include <algorithm>
#include <numeric>
#include <fstream>
#include <filesystem>
#include <set>
#include <stdexcept>

namespace pgphase_collect {

bool gap_owns_variant(const PhaseGap& gap, const VariantKey& key) {
    if (key.tid != gap.tid) return false;
    const hts_pos_t first = key.sort_pos();
    const hts_pos_t last = key.ref_len > 0 ? key.pos + key.ref_len - 1 : first;
    return first > gap.left_end && last < gap.right_beg;
}

static ReadRecord copy_evidence_read(const ReadRecord& source) {
    ReadRecord read;
    read.tid = source.tid;
    read.input_index = source.input_index;
    read.beg = source.beg;
    read.end = source.end;
    read.reverse = source.reverse;
    read.nm = source.nm;
    read.mapq = source.mapq;
    read.qname = source.qname;
    if (source.alignment) read.alignment.reset(bam_dup1(source.alignment.get()));
    read.qual = source.qual;
    read.digars = source.digars;
    read.noisy_regions = source.noisy_regions;
    read.is_skipped = source.is_skipped;
    read.is_ont_palindrome = source.is_ont_palindrome;
    read.total_cand_events = source.total_cand_events;
    return read;
}

static GapAllele evidence_allele(int allele, int qi, const ReadRecord& read,
                                 const GapEvent& event, bool binary_source) {
    GapAllele call;
    if (read.alignment && qi >= 0 && qi < read.alignment->core.l_qseq) {
        call.query_index = qi;
        const int quality = bam_get_qual(read.alignment.get())[qi];
        if (quality != 255) call.base_quality = quality;
    }
    if (allele < 0) {
        call.status = allele == -2 ? GapAlleleStatus::LowQuality :
                      allele == -3 ? GapAlleleStatus::Conflicting : GapAlleleStatus::Missing;
        return call;
    }
    // GAF/BAM binary allele 1 spells key.alt, which need not be the first
    // alternate in an MSA dictionary. Never equate source-local integers.
    if (binary_source && allele == 1) {
        const auto found = std::find(event.alleles.begin(), event.alleles.end(), event.prototype.key.alt);
        allele = found == event.alleles.end() ? -1 : static_cast<int>(found - event.alleles.begin());
    }
    if (allele < 0 || static_cast<size_t>(allele) >= event.alleles.size()) {
        call.status = GapAlleleStatus::Conflicting;
        return call;
    }
    call.status = GapAlleleStatus::Observed;
    call.index = allele;
    return call;
}

GapEvidence::GapEvidence(PhasingChunk input, const PhaseGap& gap)
    : gap_(gap), region_(input.region), ref_beg_(input.ref_beg), ref_end_(input.ref_end),
      reference_(std::move(input.ref_seq)), low_complexity_(std::move(input.low_complexity_regions)),
      noisy_(std::move(input.noisy_regions)), reads_(std::move(input.reads)) {
    std::set<std::pair<int, std::string>> molecules;
    for (const auto& read : reads_)
        if (!molecules.emplace(read.input_index, read.qname).second)
            throw std::runtime_error("duplicate molecule in gap evidence: " + read.qname);
    std::vector<int> event_index(input.candidates.size(), -1);
    for (size_t vi = 0; vi < input.candidates.size(); ++vi) {
        const auto& candidate = input.candidates[vi];
        const bool owned = gap_owns_variant(gap, candidate.key);
        const bool anchor = candidate.graph_site &&
            (candidate.phase_set == gap.left_ps || candidate.phase_set == gap.right_ps) &&
            candidate.hap_to_cons_alle[1] >= 0 && candidate.hap_to_cons_alle[2] >= 0 &&
            candidate.hap_to_cons_alle[1] != candidate.hap_to_cons_alle[2];
        const bool overlaps = candidate.key.pos + candidate.key.ref_len > gap.left_end &&
                              candidate.key.sort_pos() < gap.right_beg;
        if (!owned && !anchor && !overlaps) continue;
        GapEvent event;
        event.beg = candidate.key.pos - 1;
        event.end = event.beg + candidate.key.ref_len;
        event.role = anchor ? GapEventRole::Anchor : owned ?
            (candidate.graph_site ? GapEventRole::Graph : GapEventRole::Private) : GapEventRole::Boundary;
        event.prototype = candidate;
        if (anchor) {
            event.anchor_ps = candidate.phase_set;
            event.anchor_alleles = {candidate.hap_to_cons_alle[1], candidate.hap_to_cons_alle[2]};
        }
        event.prototype.phase_set = -1;
        event.prototype.hap_alt = event.prototype.hap_ref = 0;
        event.prototype.hap_to_cons_alle = {-1, -1, -1};
        for (auto& profile : event.prototype.hap_to_alle_profile) profile.clear();
        event.prototype.gap_link_supported = false;
        const bool supported = (candidate.key.type == VariantType::Snp && candidate.key.ref_len == 1 && candidate.key.alt.size() == 1) ||
            (candidate.key.type == VariantType::Insertion && candidate.key.ref_len == 0) ||
            (candidate.key.type == VariantType::Deletion && candidate.key.alt.empty());
        if (event.beg < ref_beg_ - 1 || event.end > ref_end_ || event.end < event.beg) {
            event.role = GapEventRole::Unsupported;
            event.alleles = {""};
        } else {
            event.alleles.push_back(reference_.substr(static_cast<size_t>(event.beg - ref_beg_ + 1),
                                                      static_cast<size_t>(event.end - event.beg)));
            if (!candidate.msa_insertion_alts.empty())
                event.alleles.insert(event.alleles.end(), candidate.msa_insertion_alts.begin(), candidate.msa_insertion_alts.end());
            else event.alleles.push_back(candidate.key.alt);
            if (!supported) event.role = GapEventRole::Unsupported;
        }
        if (event.alleles.size() > 1) {
            std::sort(event.alleles.begin() + 1, event.alleles.end());
            event.alleles.erase(std::unique(event.alleles.begin() + 1, event.alleles.end()), event.alleles.end());
            const auto primary = std::find(event.alleles.begin() + 1, event.alleles.end(), candidate.key.alt);
            if (primary != event.alleles.end()) std::rotate(event.alleles.begin() + 1, primary, primary + 1);
        }
        event.msa_to_event.push_back(0);
        const std::vector<std::string> source_alts = candidate.msa_insertion_alts.empty()
            ? std::vector<std::string>{candidate.key.alt} : candidate.msa_insertion_alts;
        for (const auto& allele : source_alts) {
            const auto found = std::find(event.alleles.begin(), event.alleles.end(), allele);
            event.msa_to_event.push_back(found == event.alleles.end() ? -1 : static_cast<int>(found - event.alleles.begin()));
        }
        // Anchor consensus uses the same source-local dictionary as its
        // working observations. Translate it before replacing that dictionary.
        if (anchor) {
            for (auto& allele : event.anchor_alleles)
                allele = allele >= 0 && static_cast<size_t>(allele) < event.msa_to_event.size()
                    ? event.msa_to_event[allele] : -1;
        }
        if (!candidate.msa_insertion_alts.empty())
            event.prototype.msa_insertion_alts.assign(event.alleles.begin() + 1, event.alleles.end());
        event.id = std::to_string(gap.tid) + ":" + std::to_string(event.beg) + ":" + std::to_string(event.end);
        for (const auto& allele : event.alleles) event.id += ":" + allele;
        event_index[vi] = static_cast<int>(events_.size());
        events_.push_back(std::move(event));
    }
    for (size_t ri = 0; ri < input.read_var_profile.size(); ++ri) {
        const auto& profile = input.read_var_profile[ri];
        for (size_t pi = 0; pi < profile.alleles.size(); ++pi) {
            const int vi = profile.start_var_idx + static_cast<int>(pi);
            if (vi < 0 || static_cast<size_t>(vi) >= event_index.size() || event_index[vi] < 0) continue;
            GapObservation observation;
            observation.read = ri;
            observation.event = static_cast<size_t>(event_index[vi]);
            const auto& event = events_[observation.event];
            const int bam = pi < profile.bam_alleles.size() ? profile.bam_alleles[pi] : -1;
            const int qi = pi < profile.bam_qi.size() ? profile.bam_qi[pi] : -1;
            const int graph = pi < profile.graph_alleles.size() ? profile.graph_alleles[pi] : -1;
            observation.bam = evidence_allele(bam, qi, reads_[ri], event, true);
            observation.graph = evidence_allele(graph, -1, reads_[ri], event, true);
            if (event.prototype.msa_verified) {
                int allele = profile.alleles[pi];
                if (allele >= 0) allele = static_cast<size_t>(allele) < event.msa_to_event.size()
                    ? event.msa_to_event[allele] : -3;
                observation.msa = evidence_allele(allele, -1, reads_[ri], event, false);
            }
            if (observation.bam.status != GapAlleleStatus::Missing ||
                observation.graph.status != GapAlleleStatus::Missing ||
                observation.msa.status != GapAlleleStatus::Missing)
                observations_.push_back(observation);
        }
    }
}

static int legacy_allele(const GapAllele& call) {
    return call.status == GapAlleleStatus::Observed ? call.index :
           call.status == GapAlleleStatus::LowQuality ? -2 :
           call.status == GapAlleleStatus::Conflicting ? -3 : -1;
}

PhasingChunk GapEvidence::project(const Options& opts, bool bam_only) const {
    PhasingChunk chunk;
    chunk.region = region_;
    chunk.ref_beg = ref_beg_;
    chunk.ref_end = ref_end_;
    chunk.ref_seq = reference_;
    chunk.low_complexity_regions = low_complexity_;
    chunk.noisy_regions = noisy_;
    std::vector<int> indices(events_.size(), -1);
    for (size_t i = 0; i < events_.size(); ++i) {
        const auto& event = events_[i];
        if (event.role == GapEventRole::Boundary || event.role == GapEventRole::Unsupported) continue;
        indices[i] = static_cast<int>(chunk.candidates.size());
        chunk.candidates.push_back(event.prototype);
        auto& candidate = chunk.candidates.back();
        if (event.role == GapEventRole::Anchor) {
            candidate.phase_set = event.anchor_ps;
            candidate.hap_to_cons_alle = {-1, event.anchor_alleles[0], event.anchor_alleles[1]};
        }
        candidate.counts.total_cov = candidate.counts.ref_cov = candidate.counts.alt_cov = 0;
        candidate.counts.low_qual_cov = 0;
        candidate.counts.forward_ref = candidate.counts.reverse_ref = 0;
        candidate.counts.forward_alt = candidate.counts.reverse_alt = 0;
        candidate.counts.alle_covs.assign(event.alleles.size(), 0);
        candidate.counts.n_uniq_alles = static_cast<int>(event.alleles.size());
    }
    for (const auto& read : reads_) chunk.reads.push_back(copy_evidence_read(read));
    chunk.read_var_profile.resize(reads_.size());
    for (size_t ri = 0; ri < reads_.size(); ++ri) {
        auto& profile = chunk.read_var_profile[ri];
        profile.read_id = static_cast<int>(ri);
        profile.start_var_idx = chunk.candidates.empty() ? -1 : 0;
        profile.end_var_idx = static_cast<int>(chunk.candidates.size()) - 1;
        profile.alleles.assign(chunk.candidates.size(), -1);
        profile.alt_qi.assign(chunk.candidates.size(), -1);
        profile.graph_alleles.assign(chunk.candidates.size(), -1);
        profile.bam_alleles.assign(chunk.candidates.size(), -1);
        profile.bam_qi.assign(chunk.candidates.size(), -1);
    }
    for (const auto& observation : observations_) {
        const int vi = indices[observation.event];
        if (vi < 0) continue;
        const auto& event = events_[observation.event];
        auto& profile = chunk.read_var_profile[observation.read];
        const bool graph_site = event.role == GapEventRole::Anchor || event.role == GapEventRole::Graph;
        const auto& selected = graph_site && !bam_only ? observation.graph :
            event.prototype.msa_verified ? observation.msa :
            observation.bam.status != GapAlleleStatus::Missing || bam_only ? observation.bam : observation.graph;
        profile.alleles[vi] = legacy_allele(selected);
        profile.alt_qi[vi] = selected.query_index;
        profile.bam_alleles[vi] = legacy_allele(observation.bam);
        profile.bam_qi[vi] = observation.bam.query_index;
        profile.graph_alleles[vi] = legacy_allele(observation.graph);
        if (!bam_only && selected.status == GapAlleleStatus::Observed &&
            observation.graph.status == GapAlleleStatus::Observed && selected.index == observation.graph.index)
            profile.alt_qi[vi] = kGraphConfirmedAltQi;
    }
    if (bam_only) select_graph_gap_bam_reads(chunk, gap_, opts);
    for (size_t ri = 0; ri < chunk.reads.size(); ++ri) {
        if (chunk.reads[ri].is_skipped) continue;
        const auto& profile = chunk.read_var_profile[ri];
        for (size_t vi = 0; vi < chunk.candidates.size(); ++vi) {
            auto& counts = chunk.candidates[vi].counts;
            const int allele = profile.alleles[vi];
            if (allele == -2) ++counts.low_qual_cov;
            if (allele < 0) continue;
            ++counts.total_cov;
            ++counts.alle_covs[allele];
            const bool reverse = chunk.reads[ri].reverse;
            if (allele == 0) reverse ? ++counts.reverse_ref : ++counts.forward_ref;
            else if (allele == 1) reverse ? ++counts.reverse_alt : ++counts.forward_alt;
        }
    }
    for (auto& candidate : chunk.candidates) {
        auto& counts = candidate.counts;
        counts.ref_cov = counts.alle_covs[0];
        counts.alt_cov = counts.alle_covs.size() > 1 ? counts.alle_covs[1] : 0;
        counts.allele_fraction = counts.total_cov ? static_cast<double>(counts.alt_cov) / counts.total_cov : 0;
        if (candidate.lcd_var_i_to_cate & (kCandHetVarCate | kCandNonAnchorHet)) {
            const int supported = static_cast<int>(std::count_if(counts.alle_covs.begin(), counts.alle_covs.end(), [&](int depth) {
                return depth >= opts.min_alt_depth && counts.total_cov > 0 &&
                       static_cast<double>(depth) / counts.total_cov >= opts.min_af;
            }));
            if (counts.total_cov < opts.min_depth || supported < 2) {
                counts.category = VariantCategory::LowCoverage;
                candidate.lcd_var_i_to_cate = kLongcalldLowCovVar;
            }
        }
    }
    chunk.haps.assign(reads_.size(), 0);
    chunk.phase_sets.assign(reads_.size(), -1);
    chunk.ordered_read_ids.resize(reads_.size());
    std::iota(chunk.ordered_read_ids.begin(), chunk.ordered_read_ids.end(), 0);
    std::sort(chunk.ordered_read_ids.begin(), chunk.ordered_read_ids.end(), [&](int a, int b) {
        return std::tie(reads_[a].beg, reads_[a].input_index, reads_[a].qname) <
               std::tie(reads_[b].beg, reads_[b].input_index, reads_[b].qname);
    });
    chunk.read_var_cr.reset(cr_init());
    for (const auto& profile : chunk.read_var_profile)
        if (profile.start_var_idx >= 0)
            cr_add(chunk.read_var_cr.get(), "cr", profile.start_var_idx, profile.end_var_idx + 1, profile.read_id);
    cr_index(chunk.read_var_cr.get());
    return chunk;
}

void GapEvidence::write_audit(const std::string& prefix) const {
    const auto event_path = prefix + ".events.tsv";
    const auto observation_path = prefix + ".observations.tsv";
    if (std::filesystem::exists(event_path) || std::filesystem::exists(observation_path))
        throw std::runtime_error("gap source audit exists; use a new directory: " + prefix);
    std::ofstream events(event_path), observations(observation_path);
    if (!events || !observations) throw std::runtime_error("cannot write gap source audit: " + prefix);
    events << "EVENT\tID\tBEGIN_0\tEND_0\tROLE\tALLELES\n";
    const char* roles[] = {"anchor", "private", "graph", "boundary", "unsupported"};
    for (size_t i = 0; i < events_.size(); ++i) {
        const auto& event = events_[i];
        events << i << '\t' << event.id << '\t' << event.beg << '\t' << event.end << '\t'
               << roles[static_cast<int>(event.role)] << '\t';
        for (size_t ai = 0; ai < event.alleles.size(); ++ai)
            events << (ai ? "," : "") << event.alleles[ai];
        events << '\n';
    }
    observations << "INPUT\tREAD\tEVENT\tBAM_ALLELE\tBAM_STATUS\tBAM_QI\tBAM_BQ"
                    "\tGRAPH_ALLELE\tGRAPH_STATUS\tMSA_ALLELE\tMSA_STATUS\n";
    const char* statuses[] = {"missing", "observed", "low_quality", "conflicting"};
    for (const auto& observation : observations_) {
        const auto& read = reads_[observation.read];
        observations << read.input_index << '\t' << read.qname << '\t' << observation.event << '\t'
                     << observation.bam.index << '\t' << statuses[static_cast<int>(observation.bam.status)] << '\t'
                     << observation.bam.query_index << '\t' << observation.bam.base_quality << '\t'
                     << observation.graph.index << '\t' << statuses[static_cast<int>(observation.graph.status)] << '\t'
                     << observation.msa.index << '\t' << statuses[static_cast<int>(observation.msa.status)] << '\n';
    }
    events.close(); observations.close();
    if (!events || !observations) throw std::runtime_error("failed writing gap source audit: " + prefix);
}

const std::vector<GapEvent>& GapEvidence::events() const { return events_; }
const std::vector<GapObservation>& GapEvidence::observations() const { return observations_; }
const PhaseGap& GapEvidence::gap() const { return gap_; }

} // namespace pgphase_collect
