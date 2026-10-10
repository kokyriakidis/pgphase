#ifndef PGPHASE_ALLELE_IDENTITY_HPP
#define PGPHASE_ALLELE_IDENTITY_HPP

#include <htslib/hts.h>

#include <array>
#include <cstddef>
#include <functional>
#include <map>
#include <optional>
#include <string>
#include <vector>

namespace pgphase_collect {

/// Exact identity within one contig; graph keys still require selected-ALT metadata.
struct CandidateIdentityKey {
    hts_pos_t pos;
    int type;
    int ref_len;
    std::string alt;
    bool operator<(const CandidateIdentityKey& other) const;
    bool operator==(const CandidateIdentityKey& other) const;
};

/// An unordered pair of complete sequence alleles relative to one reference.
/// A null member is literal reference, never an unknown observation.
struct AlleleContrastKey {
    std::array<std::optional<CandidateIdentityKey>, 2> alleles;
    bool operator<(const AlleleContrastKey& other) const;
    bool operator==(const AlleleContrastKey& other) const;
};

struct NormalizedAlleleContrast {
    AlleleContrastKey key;
    bool reversed;
};

/// Validate both full alleles; retain ALT/ALT and report their canonical order.
std::optional<NormalizedAlleleContrast> normalize_allele_contrast(
    hts_pos_t pos, const std::string& ref, const std::string& first,
    const std::string& second, const std::function<char(hts_pos_t)>& reference_base);

/// Multiple rows claiming one identity remain ambiguous regardless of insertion order.
class CandidateIdentityIndex {
public:
    void insert(const CandidateIdentityKey& key, size_t index);
    std::optional<size_t> find(const CandidateIdentityKey& key) const;
    const std::map<CandidateIdentityKey, std::optional<size_t>>& entries() const;
private:
    std::map<CandidateIdentityKey, std::optional<size_t>> entries_;
};

struct CandidateIdentityMatch {
    size_t index;
    bool is_raw;
    bool is_normalized = false;
};

struct SourceAlleleObservation {
    size_t candidate_index;
    int allele;
    int query_index;
};

/// One molecule and allele contrast: source descriptions never add voting weight.
/// Opposing calls remain conflicting; query coordinates survive only if all agree.
class MoleculeAlleleEvidence {
public:
    void add(size_t candidate_index, int allele, int query_index);
    int allele() const;
    int query_index() const;
    const std::vector<SourceAlleleObservation>& observations() const;
private:
    int allele_ = -1;
    int query_index_ = 0;
    std::vector<SourceAlleleObservation> observations_;
};

/// Raw and selected-sequence identities must agree when both contain the key.
std::optional<CandidateIdentityMatch> match_candidate_identity(
    const CandidateIdentityIndex& raw_index,
    const CandidateIdentityIndex& sequence_index,
    const CandidateIdentityKey& key);


/// Validate REF and left-align an edit using 1-based reference bases (N means unknown).
/// Empty REF/ALT supports physical insertions/deletions; no reference allele is returned.
std::optional<CandidateIdentityKey> normalize_candidate_identity(
    hts_pos_t pos, std::string ref, std::string alt,
    const std::function<char(hts_pos_t)>& reference_base);

/// Exact matches cannot bypass a conflicting or ambiguous normalized identity.
std::optional<CandidateIdentityMatch> match_normalized_candidate_identity(
    const CandidateIdentityIndex& raw_index,
    const CandidateIdentityIndex& sequence_index,
    const CandidateIdentityIndex& normalized_index,
    const CandidateIdentityKey& key,
    const std::optional<CandidateIdentityKey>& normalized_key);

/// A normalized alias requires an independent graph call and coverage of both edits.
bool supports_normalized_observation(int graph_allele,
                          hts_pos_t read_beg, hts_pos_t read_end,
                          hts_pos_t source_pos, size_t source_ref_len,
                          hts_pos_t parent_pos, size_t parent_ref_len);

} // namespace pgphase_collect

#endif
