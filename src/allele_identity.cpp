#include "allele_identity.hpp"
#include "phasing_types.hpp"

#include <algorithm>
#include <cctype>
#include <tuple>

namespace pgphase_collect {

bool CandidateIdentityKey::operator<(const CandidateIdentityKey& other) const {
    return std::tie(pos, type, ref_len, alt) <
           std::tie(other.pos, other.type, other.ref_len, other.alt);
}

bool CandidateIdentityKey::operator==(const CandidateIdentityKey& other) const {
    return pos == other.pos && type == other.type && ref_len == other.ref_len &&
           alt == other.alt;
}

bool AlleleContrastKey::operator<(const AlleleContrastKey& other) const {
    return alleles < other.alleles;
}

bool AlleleContrastKey::operator==(const AlleleContrastKey& other) const {
    return alleles == other.alleles;
}

void CandidateIdentityIndex::insert(const CandidateIdentityKey& key, size_t index) {
    const auto [entry, inserted] = entries_.emplace(key, index);
    if (!inserted && entry->second != std::optional<size_t>(index))
        entry->second.reset();
}

std::optional<size_t> CandidateIdentityIndex::find(const CandidateIdentityKey& key) const {
    const auto entry = entries_.find(key);
    return entry == entries_.end() ? std::nullopt : entry->second;
}

const std::map<CandidateIdentityKey, std::optional<size_t>>&
CandidateIdentityIndex::entries() const {
    return entries_;
}

void MoleculeAlleleEvidence::add(size_t candidate_index, int allele, int query_index) {
    if (allele != 0 && allele != 1 && allele != kConflictingBamAllele) return;
    for (const SourceAlleleObservation& previous : observations_) {
        if (previous.candidate_index == candidate_index && previous.allele == allele &&
            previous.query_index == query_index) return;
    }
    observations_.push_back({candidate_index, allele, query_index});
    if (allele == kConflictingBamAllele) {
        allele_ = kConflictingBamAllele;
        query_index_ = 0;
    } else if (allele_ == -1) {
        allele_ = allele;
        query_index_ = query_index;
    } else if (allele_ != allele) {
        allele_ = kConflictingBamAllele;
        query_index_ = 0;
    } else if (query_index_ != query_index) {
        // A repeat-shifted edit can occur at a different query coordinate.
        // The index is not a quality score and cannot be maximized or added.
        query_index_ = 0;
    }
}

int MoleculeAlleleEvidence::allele() const { return allele_; }
int MoleculeAlleleEvidence::query_index() const { return query_index_; }
const std::vector<SourceAlleleObservation>& MoleculeAlleleEvidence::observations() const {
    return observations_;
}

std::optional<CandidateIdentityMatch> match_candidate_identity(
        const CandidateIdentityIndex& raw_index,
        const CandidateIdentityIndex& sequence_index,
        const CandidateIdentityKey& key) {
    const auto raw = raw_index.entries().find(key);
    const auto sequence = sequence_index.entries().find(key);
    if (raw != raw_index.entries().end()) {
        if (!raw->second || (sequence != sequence_index.entries().end() &&
                             sequence->second != raw->second))
            return std::nullopt;
        return CandidateIdentityMatch{*raw->second, true};
    }
    if (sequence != sequence_index.entries().end() && sequence->second)
        return CandidateIdentityMatch{*sequence->second, false};
    return std::nullopt;
}


static char identity_base(char base) {
    base = static_cast<char>(std::toupper(static_cast<unsigned char>(base)));
    return base == 'A' || base == 'C' || base == 'G' || base == 'T' ? base : 'N';
}

std::optional<CandidateIdentityKey> normalize_candidate_identity(
        hts_pos_t pos, std::string ref, std::string alt,
        const std::function<char(hts_pos_t)>& reference_base) {
    if (pos < 1) return std::nullopt;
    for (size_t i = 0; i < ref.size(); ++i) {
        ref[i] = identity_base(ref[i]);
        if (ref[i] == 'N' || ref[i] != identity_base(reference_base(pos + i)))
            return std::nullopt;
    }
    for (char& base : alt) {
        base = identity_base(base);
        if (base == 'N') return std::nullopt;
    }
    if (ref == alt) return std::nullopt;
    // Remove common suffix first. Extending an empty allele on the left then
    // trimming again rotates repeat motifs, including padded replacements.
    for (;;) {
        while (!ref.empty() && !alt.empty() && ref.back() == alt.back()) {
            ref.pop_back();
            alt.pop_back();
        }
        if ((!ref.empty() && !alt.empty()) || pos == 1) break;
        const char previous = identity_base(reference_base(pos - 1));
        if (previous == 'N') return std::nullopt;
        --pos;
        ref.insert(ref.begin(), previous);
        alt.insert(alt.begin(), previous);
    }
    size_t prefix = 0;
    while (prefix < std::min(ref.size(), alt.size()) && ref[prefix] == alt[prefix])
        ++prefix;
    pos += prefix;
    ref.erase(0, prefix);
    alt.erase(0, prefix);
    const VariantType type = ref.size() == 1 && alt.size() == 1 ? VariantType::Snp
        : alt.size() > ref.size() ? VariantType::Insertion : VariantType::Deletion;
    return CandidateIdentityKey{type == VariantType::Snp ? pos : pos - 1,
                               static_cast<int>(type), static_cast<int>(ref.size()), alt};
}

std::optional<NormalizedAlleleContrast> normalize_allele_contrast(
        hts_pos_t pos, const std::string& ref, const std::string& first,
        const std::string& second, const std::function<char(hts_pos_t)>& reference_base) {
    if (pos < 1 || ref.empty()) return std::nullopt;
    std::string validated_ref = ref;
    for (size_t i = 0; i < ref.size(); ++i) {
        validated_ref[i] = identity_base(ref[i]);
        if (validated_ref[i] == 'N' ||
            validated_ref[i] != identity_base(reference_base(pos + i))) return std::nullopt;
    }
    AlleleContrastKey key;
    size_t index = 0;
    for (const std::string* allele : {&first, &second}) {
        std::string sequence = *allele;
        if (sequence.empty()) return std::nullopt;
        for (char& base : sequence) {
            base = identity_base(base);
            if (base == 'N') return std::nullopt;
        }
        if (sequence != validated_ref) {
            key.alleles[index] = normalize_candidate_identity(pos, validated_ref, sequence, reference_base);
            if (!key.alleles[index]) return std::nullopt;
        }
        ++index;
    }
    if (key.alleles[0] == key.alleles[1]) return std::nullopt;
    const bool reversed = key.alleles[1] < key.alleles[0];
    if (reversed) std::swap(key.alleles[0], key.alleles[1]);
    return NormalizedAlleleContrast{std::move(key), reversed};
}

std::optional<CandidateIdentityMatch> match_normalized_candidate_identity(
        const CandidateIdentityIndex& raw_index,
        const CandidateIdentityIndex& sequence_index,
        const CandidateIdentityIndex& normalized_index,
        const CandidateIdentityKey& key,
        const std::optional<CandidateIdentityKey>& normalized_key) {
    // Invalid reference context is never sufficient to certify a new alias.
    if (!normalized_key) return match_candidate_identity(raw_index, sequence_index, key);
    const auto normalized = normalized_index.entries().find(*normalized_key);
    const bool exact_present = raw_index.entries().count(key) != 0 ||
                               sequence_index.entries().count(key) != 0;
    const auto exact = match_candidate_identity(raw_index, sequence_index, key);
    if (normalized != normalized_index.entries().end()) {
        if (!normalized->second || (exact_present &&
            (!exact || exact->index != *normalized->second))) return std::nullopt;
        if (exact) return exact;
        return CandidateIdentityMatch{*normalized->second, false, true};
    }
    return exact;
}

bool supports_normalized_observation(int graph_allele,
        hts_pos_t read_beg, hts_pos_t read_end,
        hts_pos_t source_pos, size_t source_ref_len,
        hts_pos_t parent_pos, size_t parent_ref_len) {
    const hts_pos_t beg = std::max<hts_pos_t>(1, std::min(source_pos, parent_pos) - 1);
    const hts_pos_t end = std::max(source_pos + static_cast<hts_pos_t>(source_ref_len),
                                 parent_pos + static_cast<hts_pos_t>(parent_ref_len));
    // Identity equivalence does not certify a source REF class. Until complete
    // allele contrasts are retained, normalized calls only supplement reads
    // already observed independently by the graph, including disagreements.
    return (graph_allele == 0 || graph_allele == 1) && read_beg <= beg && read_end >= end;
}

} // namespace pgphase_collect
