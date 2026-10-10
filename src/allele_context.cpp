#include "allele_context.hpp"
#include "phasing_types.hpp"

#include "edlib.h"

#include <algorithm>
#include <cctype>
#include <set>
#include <stdexcept>
#include <tuple>

namespace pgphase_collect {

// Bound diagnostic work; unsupported large/edge contexts remain unknown.
static constexpr hts_pos_t kAlleleContextFlank = 16;
static constexpr size_t kMaxAlleleContextLength = 4096;
static constexpr size_t kMaxAllelePathSites = 64;

static char context_base(char base) {
    base = static_cast<char>(std::toupper(static_cast<unsigned char>(base)));
    return base == 'A' || base == 'C' || base == 'G' || base == 'T' ? base : 'N';
}

static hts_pos_t edit_beg(const CandidateIdentityKey& edit) {
    return edit.pos + (edit.type == static_cast<int>(VariantType::Snp) ? 0 : 1);
}

ReadAlleleCatalog build_read_allele_catalog(
        const std::vector<MoleculeAlleleSequence>& observations,
        const std::string& excluded_molecule, size_t max_hypotheses) {
    ReadAlleleCatalog result;
    std::map<std::string, std::string> molecules;
    std::set<std::string> conflicts, unsupported;
    for (const auto& observation : observations) {
        if (!excluded_molecule.empty() && observation.molecule == excluded_molecule) continue;
        std::string sequence = observation.sequence;
        bool valid = !observation.molecule.empty() && !sequence.empty() && sequence.size() <= kMaxAlleleContextLength;
        for (char& base : sequence) {
            base = context_base(base);
            if (base == 'N') valid = false;
        }
        if (!valid) { unsupported.insert(observation.molecule); continue; }
        const auto [entry, inserted] = molecules.emplace(observation.molecule, sequence);
        if (!inserted && entry->second != sequence) conflicts.insert(observation.molecule);
    }
    std::map<std::string, std::vector<std::string>> supporters;
    for (const auto& [molecule, sequence] : molecules) {
        if (unsupported.count(molecule) || conflicts.count(molecule)) continue;
        ++result.eligible_molecules;
        supporters[sequence].push_back(molecule);
    }
    result.unsupported_molecules.assign(unsupported.begin(), unsupported.end());
    for (const auto& molecule : conflicts)
        if (!unsupported.count(molecule)) result.conflicting_molecules.push_back(molecule);
    for (auto& [sequence, names] : supporters) {
        if (names.size() < kMinReadAlleleSupport) continue;
        if (result.hypotheses.size() == max_hypotheses) {
            result.status = ReadAlleleCatalogStatus::Limited;
            result.hypotheses.clear();
            return result;
        }
        result.hypotheses.push_back({sequence, std::move(names)});
    }
    return result;
}

AlleleReferenceMap map_allele_reference(
        const std::string& reference, const std::string& allele, size_t max_cells) {
    AlleleReferenceMap result;
    if (reference.size() > kMaxAlleleContextLength || allele.size() > kMaxAlleleContextLength)
        return result;
    std::string ref = reference, alt = allele;
    for (auto* sequence : {&ref, &alt})
        for (char& base : *sequence) {
            base = context_base(base);
            if (base == 'N') return result;
        }
    const size_t rows = ref.size() + 1, cols = alt.size() + 1;
    result.reference = ref;
    result.allele = alt;
    if (rows * cols > max_cells) {
        result.status = AlleleMapStatus::Limited;
        return result;
    }
    // Backward costs plus a rolling forward row identify all optimal consuming
    // edges, without trusting one arbitrary traceback through a repeat.
    std::vector<uint16_t> backward(rows * cols);
    const auto back = [&](size_t i, size_t j) -> uint16_t& { return backward[i * cols + j]; };
    for (size_t i = rows; i-- > 0;)
        for (size_t j = cols; j-- > 0;) {
            if (i == ref.size()) back(i, j) = alt.size() - j;
            else if (j == alt.size()) back(i, j) = ref.size() - i;
            else back(i, j) = std::min({back(i + 1, j) + 1, back(i, j + 1) + 1,
                                      back(i + 1, j + 1) + (ref[i] != alt[j])});
        }
    result.status = AlleleMapStatus::Complete;
    result.distance = back(0, 0);
    std::vector<uint16_t> forward(cols), next(cols);
    for (size_t j = 0; j < cols; ++j) forward[j] = j;
    for (size_t i = 0; i < ref.size(); ++i) {
        size_t consuming_edges = 0, matched_alt = cols;
        // Every complete alignment consumes each reference base exactly once.
        // A unique optimal consuming edge must therefore occur on every path.
        for (size_t j = 0; j < cols; ++j) {
            if (forward[j] + 1 + back(i + 1, j) == result.distance) ++consuming_edges;
            if (j < alt.size() && forward[j] + (ref[i] != alt[j]) + back(i + 1, j + 1) == result.distance) {
                ++consuming_edges;
                if (ref[i] == alt[j]) matched_alt = j;
            }
        }
        if (consuming_edges == 1 && matched_alt < cols) {
            if (!result.matches.empty() && result.matches.back().ref_beg + result.matches.back().length == i &&
                result.matches.back().alt_beg + result.matches.back().length == matched_alt)
                ++result.matches.back().length;
            else result.matches.push_back({i, matched_alt, 1});
        }
        next[0] = i + 1;
        for (size_t j = 1; j < cols; ++j)
            next[j] = std::min({forward[j] + 1, next[j - 1] + 1,
                               forward[j - 1] + (ref[i] != alt[j - 1])});
        forward.swap(next);
    }
    return result;
}

std::optional<AlleleSequenceContext> build_allele_sequence_context(
        const AlleleContrastKey& contrast, const std::vector<AlleleSiteContext>& sites,
        const std::function<char(hts_pos_t)>& reference_base) {
    if (sites.empty() || contrast.alleles[0] == contrast.alleles[1]) return std::nullopt;
    hts_pos_t beg = sites.front().pos;
    hts_pos_t end = beg;
    for (const auto& site : sites) {
        if (site.pos < 1 || site.ref.empty() || site.ref.size() > kMaxAlleleContextLength)
            return std::nullopt;
        beg = std::min(beg, site.pos);
        end = std::max(end, site.pos + static_cast<hts_pos_t>(site.ref.size()) - 1);
    }
    for (const auto& allele : contrast.alleles) {
        if (!allele) continue;
        if (allele->ref_len < 0 || allele->alt.size() > kMaxAlleleContextLength)
            return std::nullopt;
        const hts_pos_t pos = edit_beg(*allele);
        beg = std::min(beg, pos);
        end = std::max(end, pos + std::max(1, allele->ref_len) - 1);
    }
    beg -= kAlleleContextFlank;
    end += kAlleleContextFlank;
    if (beg < 1 || end - beg + 1 > static_cast<hts_pos_t>(kMaxAlleleContextLength))
        return std::nullopt;
    std::string reference;
    for (hts_pos_t pos = beg; pos <= end; ++pos) {
        const char base = context_base(reference_base(pos));
        if (base == 'N') return std::nullopt;
        reference.push_back(base);
    }
    AlleleSequenceContext result{beg, end, {reference, reference}, {}};
    for (size_t i = 0; i < contrast.alleles.size(); ++i) {
        if (!contrast.alleles[i]) continue;
        const auto& edit = *contrast.alleles[i];
        result.alleles[i].replace(static_cast<size_t>(edit_beg(edit) - beg), edit.ref_len, edit.alt);
        if (result.alleles[i].size() > kMaxAlleleContextLength) return std::nullopt;
    }
    const auto add_other = [&](const std::string& sequence) {
        if (sequence != result.alleles[0] && sequence != result.alleles[1] &&
            std::find(result.other_alleles.begin(), result.other_alleles.end(), sequence) ==
                result.other_alleles.end()) result.other_alleles.push_back(sequence);
    };
    add_other(reference);
    for (const auto& site : sites) {
        const size_t offset = static_cast<size_t>(site.pos - beg);
        for (size_t i = 0; i < site.ref.size(); ++i)
            if (context_base(site.ref[i]) != reference[offset + i]) return std::nullopt;
        for (std::string alt : site.alts) {
            if (alt.empty() || alt.size() > kMaxAlleleContextLength) return std::nullopt;
            for (char& base : alt) {
                base = context_base(base);
                if (base == 'N') return std::nullopt;
            }
            std::string sequence = reference;
            sequence.replace(offset, site.ref.size(), alt);
            if (sequence.size() > kMaxAlleleContextLength) return std::nullopt;
            add_other(sequence);
        }
    }
    return result;
}

std::optional<AlleleReadSlice> extract_allele_read_slice(
        const bam1_t* alignment, const AlleleSequenceContext& context) {
    if (alignment == nullptr || (alignment->core.flag & BAM_FUNMAP) != 0 ||
        context.beg < 1 || context.end < context.beg || alignment->core.l_qseq <= 0)
        return std::nullopt;
    hts_pos_t ref_pos = alignment->core.pos + 1;
    int query_pos = 0;
    int query_beg = -1;
    int query_end = -1;
    const uint32_t* cigar = bam_get_cigar(alignment);
    for (uint32_t i = 0; i < alignment->core.n_cigar; ++i) {
        const int op = bam_cigar_op(cigar[i]);
        const int length = bam_cigar_oplen(cigar[i]);
        if (op == BAM_CREF_SKIP && ref_pos <= context.end && ref_pos + length > context.beg)
            return std::nullopt;
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            if (ref_pos <= context.beg && context.beg < ref_pos + length)
                query_beg = query_pos + static_cast<int>(context.beg - ref_pos);
            if (ref_pos <= context.end && context.end < ref_pos + length)
                query_end = query_pos + static_cast<int>(context.end - ref_pos) + 1;
        }
        if ((bam_cigar_type(op) & 1) != 0) query_pos += length;
        if ((bam_cigar_type(op) & 2) != 0) ref_pos += length;
    }
    if (query_beg < 0 || query_end <= query_beg || query_end > alignment->core.l_qseq ||
        static_cast<size_t>(query_end - query_beg) > kMaxAlleleContextLength) return std::nullopt;
    AlleleReadSlice result{query_beg, query_end, alignment->core.qual, {}, {}};
    const uint8_t* sequence = bam_get_seq(alignment);
    const uint8_t* qualities = bam_get_qual(alignment);
    for (int i = query_beg; i < query_end; ++i) {
        const char base = context_base(seq_nt16_str[bam_seqi(sequence, i)]);
        if (base == 'N') return std::nullopt;
        result.sequence.push_back(base);
        result.qualities.push_back(qualities[i]);
    }
    return result;
}

std::map<std::string, AlleleReadSlice> collect_allele_read_slices(
        const std::vector<const bam1_t*>& alignments, const AlleleSequenceContext& context,
        int min_mapq, bool include_filtered) {
    std::map<std::string, const bam1_t*> unique;
    for (const bam1_t* alignment : alignments) {
        if (alignment == nullptr || alignment->core.tid < 0 ||
            (alignment->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) != 0 ||
            (!include_filtered && (alignment->core.flag & (BAM_FQCFAIL | BAM_FDUP)) != 0) ||
            alignment->core.qual < min_mapq) continue;
        const auto [entry, inserted] = unique.emplace(bam_get_qname(alignment), alignment);
        if (!inserted) entry->second = nullptr;
    }
    std::map<std::string, AlleleReadSlice> result;
    for (const auto& [name, alignment] : unique) {
        if (alignment == nullptr || alignment->core.pos + 1 > context.beg ||
            bam_endpos(alignment) < context.end) continue;
        auto slice = extract_allele_read_slice(alignment, context);
        if (slice) result.emplace(name, std::move(*slice));
    }
    return result;
}

static int context_distance(const std::string& query, const std::string& allele) {
    const auto result = edlibAlign(query.data(), query.size(), allele.data(), allele.size(),
        edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, nullptr, 0));
    const int distance = result.editDistance;
    const int status = result.status;
    edlibFreeAlignResult(result);
    if (status != EDLIB_STATUS_OK || distance < 0)
        throw std::runtime_error("cannot score allele sequence context with edlib");
    return distance;
}

std::vector<std::string> full_allele_sequences(const AlleleSequenceContext& context) {
    std::vector<std::string> result = context.other_alleles;
    result.insert(result.end(), context.alleles.begin(), context.alleles.end());
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    return result;
}

bool PhysicalAlleleContextKey::operator<(const PhysicalAlleleContextKey& other) const {
    return std::tie(beg, end, alleles) < std::tie(other.beg, other.end, other.alleles);
}

std::map<PhysicalAlleleContextKey, std::vector<size_t>> group_allele_contexts(
        const std::vector<std::optional<AlleleSequenceContext>>& contexts) {
    std::map<PhysicalAlleleContextKey, std::vector<size_t>> groups;
    for (size_t i = 0; i < contexts.size(); ++i) {
        if (!contexts[i]) continue;
        const auto& context = *contexts[i];
        groups[{context.beg, context.end, full_allele_sequences(context)}].push_back(i);
    }
    return groups;
}

AlleleComposition compose_allele_sequence(
        hts_pos_t beg, hts_pos_t end, const std::vector<AlleleSequenceEdit>& edits,
        const std::function<char(hts_pos_t)>& reference_base) {
    if (beg < 1 || end < beg || end - beg + 1 > static_cast<hts_pos_t>(kMaxAlleleContextLength))
        return {};
    std::string reference;
    for (hts_pos_t pos = beg; pos <= end; ++pos) {
        const char base = context_base(reference_base(pos));
        if (base == 'N') return {};
        reference.push_back(base);
    }
    std::vector<CandidateIdentityKey> minimal;
    struct AmbiguousEdit {
        CandidateIdentityKey edit;
        hts_pos_t beg;
        hts_pos_t end;
    };
    std::vector<AmbiguousEdit> ambiguous;
    std::set<std::tuple<hts_pos_t, std::string, std::string>> anchored_edits;
    for (const auto& edit : edits) {
        if (edit.pos < beg || edit.pos > end + 1 ||
            edit.ref.size() > kMaxAlleleContextLength || edit.alt.size() > kMaxAlleleContextLength ||
            edit.pos + static_cast<hts_pos_t>(edit.ref.size()) > end + 1) return {};
        std::string ref = edit.ref;
        std::string alt = edit.alt;
        for (size_t i = 0; i < ref.size(); ++i) {
            ref[i] = context_base(ref[i]);
            if (ref[i] == 'N' || ref[i] != context_base(reference_base(edit.pos + i))) return {};
        }
        for (char& base : alt) {
            base = context_base(base);
            if (base == 'N') return {};
        }
        if (ref == alt) continue;
        if (!anchored_edits.emplace(edit.pos, ref, alt).second) continue;
        size_t anchored_prefix = 0;
        while (anchored_prefix < std::min(ref.size(), alt.size()) &&
            ref[anchored_prefix] == alt[anchored_prefix]) ++anchored_prefix;
        // Repeat shifts preserve one edit on REF but can change its meaning
        // beside another edit. Trim padding without changing the anchored path.
        while (!ref.empty() && !alt.empty() && ref.back() == alt.back()) {
            ref.pop_back();
            alt.pop_back();
        }
        size_t prefix = 0;
        while (prefix < std::min(ref.size(), alt.size()) && ref[prefix] == alt[prefix]) ++prefix;
        const hts_pos_t pos = edit.pos + prefix;
        ref.erase(0, prefix);
        alt.erase(0, prefix);
        const VariantType type = ref.size() == 1 && alt.size() == 1 ? VariantType::Snp
            : alt.size() > ref.size() ? VariantType::Insertion : VariantType::Deletion;
        minimal.push_back({type == VariantType::Snp ? pos : pos - 1,
            static_cast<int>(type), static_cast<int>(ref.size()), alt});
        if (edit.pos + anchored_prefix != static_cast<size_t>(pos))
            ambiguous.push_back({minimal.back(), pos,
                edit.pos + static_cast<hts_pos_t>(anchored_prefix + ref.size())});
    }
    std::sort(minimal.begin(), minimal.end(), [](const auto& left, const auto& right) {
        return std::make_tuple(edit_beg(left), left.ref_len, left.alt) <
               std::make_tuple(edit_beg(right), right.ref_len, right.alt);
    });
    for (const auto& region : ambiguous) {
        if (std::count(minimal.begin(), minimal.end(), region.edit) > 1)
            return {AlleleCompositionStatus::Overlap, {}};
        for (const auto& edit : minimal) {
            if (edit == region.edit) continue;
            const hts_pos_t pos = edit_beg(edit);
            if (pos <= region.end && pos + edit.ref_len >= region.beg)
                return {AlleleCompositionStatus::Overlap, {}};
        }
    }
    minimal.erase(std::unique(minimal.begin(), minimal.end()), minimal.end());
    std::map<std::string, CandidateIdentityKey> effects;
    for (const auto& edit : minimal) {
        std::string sequence = reference;
        sequence.replace(edit_beg(edit) - beg, edit.ref_len, edit.alt);
        const auto [found, inserted] = effects.emplace(sequence, edit);
        if (!inserted && !(found->second == edit))
            return {AlleleCompositionStatus::Overlap, {}};
    }
    std::string result;
    hts_pos_t cursor = beg;
    std::optional<hts_pos_t> insertion;
    for (const auto& edit : minimal) {
        const hts_pos_t pos = edit_beg(edit);
        // Boundary insertions can precede/follow a replacement; interior insertions
        // and different insertions at one boundary have no defined composition.
        if (pos < cursor || (edit.ref_len == 0 && insertion == pos))
            return {AlleleCompositionStatus::Overlap, {}};
        result.append(reference, cursor - beg, pos - cursor);
        result += edit.alt;
        if (edit.ref_len == 0) insertion = pos;
        cursor = pos + edit.ref_len;
        if (result.size() > kMaxAlleleContextLength) return {};
    }
    result.append(reference, cursor - beg, end + 1 - cursor);
    if (result.size() > kMaxAlleleContextLength) return {};
    return {AlleleCompositionStatus::Valid, std::move(result)};
}

AlleleComposition compose_allele_on_subpaths(
        hts_pos_t beg, const AlleleReferenceMap& parent,
        const std::vector<AlleleSequenceEdit>& neighbors) {
    if (beg < 1 || parent.reference.empty()) return {};
    const hts_pos_t end = beg + static_cast<hts_pos_t>(parent.reference.size()) - 1;
    const auto reference_base = [&](hts_pos_t pos) {
        return pos >= beg && pos <= end ? parent.reference[pos - beg] : 'N';
    };
    std::vector<AlleleSequenceEdit> original{{beg, parent.reference, parent.allele}};
    original.insert(original.end(), neighbors.begin(), neighbors.end());
    const auto anchored = compose_allele_sequence(beg, end, original, reference_base);
    if (anchored.status != AlleleCompositionStatus::Overlap) return anchored;
    if (parent.status != AlleleMapStatus::Complete) return {};
    const auto independent = compose_allele_sequence(beg, end, neighbors, reference_base);
    if (independent.status != AlleleCompositionStatus::Valid) return independent;
    std::vector<AlleleSequenceEdit> mapped;
    for (const auto& neighbor : neighbors) {
        std::string ref = neighbor.ref, alt = neighbor.alt;
        for (char& base : ref) base = context_base(base);
        for (char& base : alt) base = context_base(base);
        if (ref == alt) continue;
        const size_t offset = static_cast<size_t>(neighbor.pos - beg);
        const auto match = std::find_if(parent.matches.begin(), parent.matches.end(), [&](const auto& span) {
            if (offset < span.ref_beg || offset + ref.size() > span.ref_beg + span.length) return false;
            // Match edges fix bases, but do not order an extra insertion beside
            // an existing parent edit at the edge of a matched subpath.
            return ref.size() == alt.size() ||
                (offset > span.ref_beg && offset + ref.size() < span.ref_beg + span.length);
        });
        if (match == parent.matches.end()) return {AlleleCompositionStatus::Overlap, {}};
        mapped.push_back({static_cast<hts_pos_t>(match->alt_beg + offset - match->ref_beg + 1),
                          std::move(ref), std::move(alt)});
    }
    const auto allele_base = [&](hts_pos_t pos) {
        return pos >= 1 && pos <= static_cast<hts_pos_t>(parent.allele.size()) ? parent.allele[pos - 1] : 'N';
    };
    return compose_allele_sequence(1, parent.allele.size(), mapped, allele_base);
}

AllelePathCatalog build_allele_sequence_paths(
        hts_pos_t beg, hts_pos_t end, const std::vector<std::string>& parents,
        const std::vector<AlleleNeighborSite>& neighbors,
        const std::function<char(hts_pos_t)>& reference_base, size_t max_prefixes) {
    AllelePathCatalog result;
    const auto original = compose_allele_sequence(beg, end, {}, reference_base);
    if (original.status != AlleleCompositionStatus::Valid || parents.empty()) return result;
    auto sites = neighbors;
    std::sort(sites.begin(), sites.end(), [](const auto& left, const auto& right) {
        return left.candidate < right.candidate;
    });
    for (size_t i = 1; i < sites.size(); ++i)
        if (sites[i - 1].candidate == sites[i].candidate) return result;
    struct NeighborGroup {
        AlleleNeighborSite site;
        std::vector<AlleleNeighborSite> members;
    };
    std::vector<NeighborGroup> groups;
    std::map<std::tuple<hts_pos_t, std::string, std::vector<std::string>>, size_t> by_site;
    for (const auto& site : sites) {
        auto alts = site.alts;
        std::sort(alts.begin(), alts.end());
        alts.erase(std::unique(alts.begin(), alts.end()), alts.end());
        const auto [entry, inserted] = by_site.emplace(std::make_tuple(site.pos, site.ref, alts), groups.size());
        if (inserted) groups.push_back({{site.candidate, site.pos, site.ref, std::move(alts)}, {}});
        groups[entry->second].members.push_back(site);
    }
    result.status = AllelePathStatus::Complete;
    if (groups.size() > kMaxAllelePathSites) {
        result.status = AllelePathStatus::Limited;
        return result;
    }
    std::vector<AlleleSequenceEdit> edits;
    std::vector<std::pair<size_t, size_t>> selected;
    std::function<void(size_t, size_t)> visit = [&](size_t parent, size_t next) {
        if (result.status == AllelePathStatus::Limited) return;
        if (result.visited_prefixes == max_prefixes) {
            result.status = AllelePathStatus::Limited;
            result.paths.clear();
            return;
        }
        ++result.visited_prefixes;
        const auto composed = compose_allele_sequence(beg, end, edits, reference_base);
        if (composed.status == AlleleCompositionStatus::Overlap) {
            ++result.overlap_prefixes;
            return;
        }
        if (next == groups.size()) {
            if (composed.status == AlleleCompositionStatus::Valid) {
                auto provenance = selected;
                std::sort(provenance.begin(), provenance.end());
                result.paths.push_back({parent, std::move(provenance), composed.sequence});
            }
            else ++result.unsupported_paths;
            return;
        }
        // Omission is a no-edit hypothesis, not evidence for a literal REF call.
        visit(parent, next + 1);
        const auto& group = groups[next];
        const auto& site = group.site;
        for (size_t alt = 0; alt < site.alts.size(); ++alt) {
            if (result.status == AllelePathStatus::Limited) return;
            edits.push_back({site.pos, site.ref, site.alts[alt]});
            const size_t previous_selected = selected.size();
            for (const auto& member : group.members) {
                const auto allele = std::find(member.alts.begin(), member.alts.end(), site.alts[alt]);
                selected.emplace_back(member.candidate, static_cast<size_t>(std::distance(member.alts.begin(), allele)) + 1);
            }
            visit(parent, next + 1);
            selected.resize(previous_selected);
            edits.pop_back();
        }
    };
    for (size_t parent = 0; parent < parents.size(); ++parent) {
        if (result.status == AllelePathStatus::Limited) break;
        edits = {{beg, original.sequence, parents[parent]}};
        visit(parent, 0);
    }
    return result;
}

std::vector<int> score_allele_sequences(
        const std::vector<std::string>& alleles, const std::string& query) {
    std::vector<int> result;
    result.reserve(alleles.size());
    for (const auto& allele : alleles) result.push_back(context_distance(query, allele));
    return result;
}

AlleleContextScore score_allele_context(
        const AlleleSequenceContext& context, const std::string& query) {
    AlleleContextScore result{{context_distance(query, context.alleles[0]),
                              context_distance(query, context.alleles[1])}, -1, -1};
    for (const auto& allele : context.other_alleles) {
        const int distance = context_distance(query, allele);
        if (result.other_distance < 0 || distance < result.other_distance)
            result.other_distance = distance;
    }
    const int nearest = result.distances[0] < result.distances[1] ? 0 : 1;
    if (result.distances[nearest] < result.distances[1 - nearest] &&
        (result.other_distance < 0 || result.distances[nearest] < result.other_distance))
        result.nearest = nearest;
    return result;
}

} // namespace pgphase_collect
