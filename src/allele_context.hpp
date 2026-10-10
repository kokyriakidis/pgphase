#ifndef PGPHASE_ALLELE_CONTEXT_HPP
#define PGPHASE_ALLELE_CONTEXT_HPP

#include "allele_identity.hpp"

#include <htslib/sam.h>

#include <map>

namespace pgphase_collect {

struct AlleleSiteContext {
    hts_pos_t pos;
    std::string ref;
    std::vector<std::string> alts;
};

struct AlleleSequenceContext {
    hts_pos_t beg;
    hts_pos_t end;
    std::array<std::string, 2> alleles;
    std::vector<std::string> other_alleles;
};

struct PhysicalAlleleContextKey {
    hts_pos_t beg;
    hts_pos_t end;
    std::vector<std::string> alleles;

    bool operator<(const PhysicalAlleleContextKey& other) const;
};

struct AlleleSequenceEdit {
    hts_pos_t pos;
    std::string ref;
    std::string alt;
};

enum class AlleleCompositionStatus { Valid, Unsupported, Overlap };

struct AlleleComposition {
    AlleleCompositionStatus status = AlleleCompositionStatus::Unsupported;
    std::string sequence;
};

struct AlleleNeighborSite {
    size_t candidate;
    hts_pos_t pos;
    std::string ref;
    std::vector<std::string> alts;
};

struct AlleleSequencePath {
    size_t parent_allele;
    // Selected candidate and one-based ALT index; omission asserts no edit.
    std::vector<std::pair<size_t, size_t>> neighbors;
    std::string sequence;
};

enum class AllelePathStatus { Complete, Limited, Unsupported };

struct AllelePathCatalog {
    AllelePathStatus status = AllelePathStatus::Unsupported;
    size_t visited_prefixes = 0;
    size_t overlap_prefixes = 0;
    size_t unsupported_paths = 0;
    std::vector<AlleleSequencePath> paths;
};

constexpr size_t kMaxAllelePathPrefixes = 65536;

enum class AlleleMapStatus { Complete, Limited, Unsupported };

struct AlleleMatchedSubpath {
    // Zero-based half-open offsets in the complete REF and parent ALT strings.
    size_t ref_beg;
    size_t alt_beg;
    size_t length;
};

struct AlleleReferenceMap {
    AlleleMapStatus status = AlleleMapStatus::Unsupported;
    int distance = -1;
    std::vector<AlleleMatchedSubpath> matches;
    // Bind the map to its validated, case-normalized input sequences.
    std::string reference;
    std::string allele;
};

constexpr size_t kMaxAlleleMapCells = 1048576;

/// Return only exact matched subpaths shared by every optimal unit-edit alignment.
/// Complete can contain unresolved positions; a bounded-out map has no matches.
AlleleReferenceMap map_allele_reference(
    const std::string& reference, const std::string& allele,
    size_t max_cells = kMaxAlleleMapCells);

/// Preserve anchored compositions; resolve remaining overlaps only through
/// a map returned by map_allele_reference for this complete parent context.
/// Length-changing edits require matched guards on both sides of their raw REF.
AlleleComposition compose_allele_on_subpaths(
    hts_pos_t beg, const AlleleReferenceMap& parent,
    const std::vector<AlleleSequenceEdit>& neighbors);

struct AlleleReadSlice {
    int query_beg;
    int query_end;
    int mapq;
    std::string sequence;
    std::vector<uint8_t> qualities;
};

struct MoleculeAlleleSequence {
    std::string molecule;
    std::string sequence;
};

struct ReadAlleleHypothesis {
    std::string sequence;
    std::vector<std::string> molecules;
};

enum class ReadAlleleCatalogStatus { Complete, Limited };

struct ReadAlleleCatalog {
    ReadAlleleCatalogStatus status = ReadAlleleCatalogStatus::Complete;
    size_t eligible_molecules = 0;
    std::vector<std::string> conflicting_molecules;
    std::vector<std::string> unsupported_molecules;
    std::vector<ReadAlleleHypothesis> hypotheses;
};

constexpr size_t kMinReadAlleleSupport = 2;
constexpr size_t kMaxReadAlleleHypotheses = 64;

/// Discover complete-context sequences repeated by distinct physical molecules.
/// Conflicting/unsupported names abstain; source aliases do not add support.
/// Exclusion removes a molecule from discovery before grouping or validation.
/// A limited catalog publishes no partial hypotheses; support is not confidence.
ReadAlleleCatalog build_read_allele_catalog(
    const std::vector<MoleculeAlleleSequence>& observations,
    const std::string& excluded_molecule = {},
    size_t max_hypotheses = kMaxReadAlleleHypotheses);

struct AlleleContextScore {
    std::array<int, 2> distances;
    int other_distance;
    // Unique closest selected allele, or -1 for a tie/closer unselected allele.
    // This is an edit-distance ranking, not a calibrated phasing confidence.
    int nearest;
};

/// Pad complete parent sites and normalized edits into one reference window.
/// Preserve every unselected parent allele as a competing hypothesis.
std::optional<AlleleSequenceContext> build_allele_sequence_context(
    const AlleleContrastKey& contrast, const std::vector<AlleleSiteContext>& sites,
    const std::function<char(hts_pos_t)>& reference_base);

/// Require aligned outer flanks; retain internal I/D, reject reference skips.
/// BAM SEQ is already in reference orientation, including reverse alignments.
std::optional<AlleleReadSlice> extract_allele_read_slice(
    const bam1_t* alignment, const AlleleSequenceContext& context);

/// Collect complete-context sequences from original alignments on one contig.
/// Admission is independent of variant calls; duplicate eligible primary names
/// remain ambiguous even if only one of their alignments covers the context.
std::map<std::string, AlleleReadSlice> collect_allele_read_slices(
    const std::vector<const bam1_t*>& alignments, const AlleleSequenceContext& context,
    int min_mapq, bool include_filtered);

/// Score a saved query slice without BAM access or source-channel preferences.
AlleleContextScore score_allele_context(
    const AlleleSequenceContext& context, const std::string& query);

/// Canonical sequence order independent of selected pair or parent ALT order.
std::vector<std::string> full_allele_sequences(const AlleleSequenceContext& context);

/// Group descriptions on one contig by identical bounds and full hypotheses.
/// Selected contrasts are provenance, not additional physical evidence.
std::map<PhysicalAlleleContextKey, std::vector<size_t>> group_allele_contexts(
    const std::vector<std::optional<AlleleSequenceContext>>& contexts);

/// Compose reference-validated edits in one complete reference window.
/// Duplicate anchored edits count once; shifted aliases and unresolved overlaps
/// return no sequence. Independently left-aligning neighbors is unsafe.
/// REF/no-op entries impose no edit, rather than a genotype constraint.
AlleleComposition compose_allele_sequence(
    hts_pos_t beg, hts_pos_t end, const std::vector<AlleleSequenceEdit>& edits,
    const std::function<char(hts_pos_t)>& reference_base);

/// Enumerate supported parent/multiple-neighbor paths, one ALT per physical site.
/// Equal raw bounds/REF and full ALT tables share a choice; retain source indices.
/// A work-limited search returns no paths; its prefix is not a complete catalog.
/// Overlap pruning is monotone, while unsupported length prefixes may recover.
AllelePathCatalog build_allele_sequence_paths(
    hts_pos_t beg, hts_pos_t end, const std::vector<std::string>& parents,
    const std::vector<AlleleNeighborSite>& neighbors,
    const std::function<char(hts_pos_t)>& reference_base,
    size_t max_prefixes = kMaxAllelePathPrefixes);

/// Retain every competing allele's distance, rather than only its minimum.
std::vector<int> score_allele_sequences(
    const std::vector<std::string>& alleles, const std::string& query);

} // namespace pgphase_collect

#endif
