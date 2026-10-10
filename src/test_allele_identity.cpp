#include "allele_identity.hpp"
#include "phasing_types.hpp"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace pgphase_collect;

static bool check(bool condition, const std::string& message) {
    if (!condition) std::cerr << "FAIL: " << message << "\n";
    return condition;
}

static bool check_molecule_evidence() {
    bool ok = true;
    struct Fixture {
        std::vector<SourceAlleleObservation> calls;
        int allele;
        int query_index;
        size_t provenance_count;
    };
    const std::vector<Fixture> fixtures = {
        {{{2, 1, 42}, {7, 1, 42}, {2, 1, 42}}, 1, 42, 2},
        {{{2, 0, 42}, {7, 0, 45}, {9, 0, 42}}, 0, 0, 3},
        {{{2, 1, 42}, {7, 0, 45}, {9, 1, 42}}, kConflictingBamAllele, 0, 3},
        {{{2, 1, 42}, {2, 0, 42}, {2, 1, 42}}, kConflictingBamAllele, 0, 2},
        {{{2, -1, 42}, {7, 1, 45}, {9, -1, 42}}, 1, 45, 1},
        {{{2, -1, 42}, {7, -1, 45}, {9, -1, 42}}, -1, 0, 0},
        {{{2, 1, 42}, {7, 1, 0}, {9, 1, 42}}, 1, 0, 3},
        {{{2, kConflictingBamAllele, 0}, {7, 1, 42}, {9, -1, 42}}, kConflictingBamAllele, 0, 2},
    };
    size_t permutations = 0;
    for (const Fixture& fixture : fixtures) {
        std::vector<size_t> order{0, 1, 2};
        do {
            MoleculeAlleleEvidence evidence;
            for (const size_t i : order) {
                const auto& call = fixture.calls[i];
                evidence.add(call.candidate_index, call.allele, call.query_index);
            }
            ok &= check(evidence.allele() == fixture.allele &&
                        evidence.query_index() == fixture.query_index &&
                        evidence.observations().size() == fixture.provenance_count,
                        "molecule: order-independent single call, query index and provenance");
            for (int repeat = 0; repeat < 10; ++repeat) {
                for (const auto& call : fixture.calls)
                    evidence.add(call.candidate_index, call.allele, call.query_index);
            }
            ok &= check(evidence.allele() == fixture.allele &&
                        evidence.query_index() == fixture.query_index &&
                        evidence.observations().size() == fixture.provenance_count,
                        "molecule: repeated descriptions cannot add weight or resolve conflict");
            for (const auto& call : fixture.calls) {
                if (call.allele != 0 && call.allele != 1 && call.allele != kConflictingBamAllele) continue;
                const auto& provenance = evidence.observations();
                ok &= check(std::any_of(provenance.begin(), provenance.end(), [&](const auto& obs) {
                    return obs.candidate_index == call.candidate_index &&
                           obs.allele == call.allele && obs.query_index == call.query_index;
                }), "molecule: original source calls survive reduction");
            }
            ++permutations;
        } while (std::next_permutation(order.begin(), order.end()));
    }
    // Cohorts are per molecule and complete allele contrast, not per position.
    const CandidateIdentityKey first{10, static_cast<int>(VariantType::Insertion), 0, "A"};
    const CandidateIdentityKey second{10, static_cast<int>(VariantType::Insertion), 0, "AA"};
    std::map<CandidateIdentityKey, MoleculeAlleleEvidence> loci;
    loci[first].add(2, 1, 42);
    loci[second].add(7, 0, 45);
    MoleculeAlleleEvidence other_read;
    other_read.add(2, 0, 51);
    ok &= check(loci[first].allele() == 1 && loci[second].allele() == 0 &&
                other_read.allele() == 0, "molecule: different alleles and reads remain independent");
    const std::string reference = "CAAATG";
    const auto base = [&](hts_pos_t pos) {
        return pos >= 1 && pos <= static_cast<hts_pos_t>(reference.size())
            ? reference[static_cast<size_t>(pos - 1)] : 'N';
    };
    const auto canonical = normalize_candidate_identity(1, "C", "CA", base);
    CandidateIdentityIndex raw, selected, normalized;
    normalized.insert(*canonical, 11);
    const CandidateIdentityKey source1{2, static_cast<int>(VariantType::Insertion), 0, "A"};
    const CandidateIdentityKey source2{3, static_cast<int>(VariantType::Insertion), 0, "A"};
    const auto parent1 = match_normalized_candidate_identity(raw, selected, normalized,
        source1, normalize_candidate_identity(3, "", "A", base));
    const auto parent2 = match_normalized_candidate_identity(raw, selected, normalized,
        source2, normalize_candidate_identity(4, "", "A", base));
    ok &= check(parent1 && parent2 && parent1->index == parent2->index,
                "molecule: repeat-shifted source descriptions share one unique destination");
    if (parent1 && parent2) {
        std::map<size_t, MoleculeAlleleEvidence> destinations;
        destinations[parent1->index].add(2, 1, 42);
        destinations[parent2->index].add(7, 1, 43);
        ok &= check(destinations.size() == 1 && destinations.begin()->second.allele() == 1 &&
                    destinations.begin()->second.query_index() == 0 &&
                    destinations.begin()->second.observations().size() == 2,
                    "molecule: normalized duplicate rows supply one vote with both sources retained");
        destinations[parent2->index].add(9, 0, 44);
        ok &= check(destinations.begin()->second.allele() == kConflictingBamAllele,
                    "molecule: opposing normalized calls abstain at their shared destination");
    }
    normalized.insert(*canonical, 12);
    ok &= check(!match_normalized_candidate_identity(raw, selected, normalized,
        source1, canonical), "molecule: source merging cannot bypass an ambiguous graph destination");
    std::cout << "Molecule permutations: " << permutations << "\n";
    return ok;
}

static bool check_normalization() {
    bool ok = true;
    const std::string reference = "CAAATATATGACCT";
    const auto base = [&](hts_pos_t pos) {
        return pos >= 1 && pos <= static_cast<hts_pos_t>(reference.size())
            ? reference[static_cast<size_t>(pos - 1)] : 'N';
    };
    const auto normalize = [&](hts_pos_t pos, const std::string& ref, const std::string& alt) {
        return normalize_candidate_identity(pos, ref, alt, base);
    };
    const auto insertion = normalize(1, "C", "CA");
    ok &= check(insertion && insertion == normalize(3, "A", "AA") &&
                insertion == normalize(4, "", "a"), "normalize: homopolymer shifts and case agree");
    ok &= check(normalize(4, "A", "ATA") == normalize(8, "A", "ATA") &&
                normalize(4, "A", "ATA") == normalize(7, "T", "TAT"),
                "normalize: tandem motif rotation agrees");
    ok &= check(normalize(1, "CA", "C") == normalize(3, "AA", "A"),
                "normalize: shifted deletions agree");
    ok &= check(normalize(10, "GACC", "GTTC") == normalize(11, "AC", "TT"),
                "normalize: compound replacements retain both changed bases");
    ok &= check(!(normalize(11, "AC", "TT") == normalize(11, "AC", "GG")),
                "normalize: same position and length do not collapse different sequences");
    ok &= check(normalize(10, "GAC", "GTC") == normalize(11, "A", "T"),
                "normalize: padded substitutions agree");
    ok &= check(!normalize(10, "AAC", "ATC") && !normalize(15, "A", "T") &&
                !normalize(0, "C", "T"), "normalize: mismatching and missing REF abstain");
    ok &= check(!normalize(1, "C", "<INS>") && !normalize(1, "C", "*") &&
                !normalize(1, "C", "N") && !normalize(1, "C", "C"),
                "normalize: symbolic, unknown and unchanged alleles abstain");
    const auto incomplete = [&](hts_pos_t pos) { return pos >= 3 ? base(pos) : 'N'; };
    ok &= check(!normalize_candidate_identity(3, "A", "AA", incomplete),
                "normalize: a repeat at an unknown left boundary cannot claim canonical identity");
    ok &= check(normalize(1, "", "C") == normalize(1, "C", "CC"),
                "normalize: contig start supports a right-anchored insertion");
    ok &= check(normalize(1, "C", "") == normalize(1, "CA", "A"),
                "normalize: contig start supports a right-anchored deletion");

    // Exhaustively compare resulting haplotypes, independently of the trim and
    // rotation algorithm. Every equivalent edit must receive the same identity.
    size_t tested = 0;
    for (int bits = 0; bits < 32; ++bits) {
        std::string seq = "G";
        for (int bit = 0; bit < 5; ++bit) seq += (bits & (1 << bit)) ? 'A' : 'C';
        seq += 'T';
        const auto ref_base = [&](hts_pos_t pos) {
            return pos >= 1 && pos <= static_cast<hts_pos_t>(seq.size())
                ? seq[static_cast<size_t>(pos - 1)] : 'N';
        };
        std::map<std::string, CandidateIdentityKey> haplotypes;
        for (size_t offset = 1; offset < seq.size(); ++offset) {
            for (size_t len = 0; len <= 2 && offset + len <= seq.size(); ++len) {
                for (const std::string alt : {"", "A", "C", "AA", "AC", "CA", "CC"}) {
                    const std::string ref = seq.substr(offset, len);
                    if (ref == alt) continue;
                    const auto id = normalize_candidate_identity(offset + 1, ref, alt, ref_base);
                    ok &= check(id.has_value(), "normalize: complete small-reference edit is valid");
                    if (!id) continue;
                    const hts_pos_t edit_pos = id->pos +
                        (id->type != static_cast<int>(VariantType::Snp));
                    const std::string haplotype = seq.substr(0, offset) + alt + seq.substr(offset + len);
                    const std::string normalized_haplotype = seq.substr(0, edit_pos - 1) + id->alt +
                        seq.substr(edit_pos - 1 + id->ref_len);
                    ok &= check(haplotype == normalized_haplotype,
                                "normalize: canonical edit preserves the full local haplotype");
                    const auto [prior, inserted] = haplotypes.emplace(haplotype, *id);
                    ok &= check(inserted || prior->second == *id,
                                "normalize: equivalent complete haplotypes have one identity");
                    const auto again = normalize_candidate_identity(edit_pos,
                        seq.substr(edit_pos - 1, id->ref_len), id->alt, ref_base);
                    ok &= check(again == id, "normalize: canonical identity is idempotent");
                    ++tested;
                }
            }
        }
    }
    const CandidateIdentityKey source{3, static_cast<int>(VariantType::Insertion), 0, "A"};
    CandidateIdentityIndex raw, sequence, normalized;
    normalized.insert(*insertion, 7);
    auto match = match_normalized_candidate_identity(raw, sequence, normalized, source, insertion);
    ok &= check(match && match->index == 7 && match->is_normalized && !match->is_raw,
                "normalize: unique shifted identity is a sequence match");
    sequence.insert(source, 8);
    ok &= check(!match_normalized_candidate_identity(raw, sequence, normalized, source, insertion),
                "normalize: exact and canonical destinations must agree");
    sequence = CandidateIdentityIndex{};
    raw.insert(source, 7);
    match = match_normalized_candidate_identity(raw, sequence, normalized, source, insertion);
    ok &= check(match && match->is_raw && !match->is_normalized,
                "normalize: agreeing indexes preserve an exact observation");
    normalized.insert(*insertion, 9);
    ok &= check(!match_normalized_candidate_identity(raw, sequence, normalized, source, insertion),
                "normalize: exact match cannot bypass shifted duplicate ambiguity");
    ok &= check(match_normalized_candidate_identity(raw, sequence, normalized, source, std::nullopt).has_value(),
                "normalize: absent canonical context retains the established exact rule");
    ok &= check(supports_normalized_observation(0, 1, 8, 4, 0, 2, 0) &&
                !supports_normalized_observation(0, 3, 8, 4, 0, 2, 0) &&
                !supports_normalized_observation(0, 1, 6, 4, 3, 2, 0),
                "normalize: short reads cannot borrow a shifted call across an uncovered flank");
    ok &= check(supports_normalized_observation(1, 1, 8, 4, 0, 2, 0) &&
                !supports_normalized_observation(-1, 1, 8, 4, 0, 2, 0) &&
                !supports_normalized_observation(2, 1, 8, 4, 0, 2, 0),
                "normalize: aliases retain independent graph calls but cannot certify a missing REF class");
    std::cout << "Normalized " << tested << " complete-haplotype fixtures\n";
    return ok;
}

static bool check_allele_contrasts() {
    bool ok = true;
    size_t tested = 0;
    for (int bits = 0; bits < 32; ++bits) {
        std::string seq = "G";
        for (int bit = 0; bit < 5; ++bit) seq += (bits & (1 << bit)) ? 'A' : 'C';
        seq += 'T';
        const auto base = [&](hts_pos_t pos) {
            return pos >= 1 && pos <= static_cast<hts_pos_t>(seq.size())
                ? seq[static_cast<size_t>(pos - 1)] : 'N';
        };
        std::map<std::array<std::string, 2>, AlleleContrastKey> pairs;
        for (size_t offset = 1; offset < seq.size(); ++offset) {
            for (size_t len = 1; len <= 2 && offset + len <= seq.size(); ++len) {
                const std::string ref = seq.substr(offset, len);
                std::vector<std::string> alleles = {ref, "A", "C", "AA", "AC", "CA", "CC"};
                std::sort(alleles.begin(), alleles.end());
                alleles.erase(std::unique(alleles.begin(), alleles.end()), alleles.end());
                for (size_t a = 0; a < alleles.size(); ++a) {
                    for (size_t b = a + 1; b < alleles.size(); ++b) {
                        const auto contrast = normalize_allele_contrast(offset + 1, ref, alleles[a], alleles[b], base);
                        ok &= check(contrast.has_value(), "contrast: complete DNA pair is valid");
                        if (!contrast) continue;
                        const auto reverse = normalize_allele_contrast(offset + 1, ref, alleles[b], alleles[a], base);
                        ok &= check(reverse && reverse->key == contrast->key &&
                            reverse->reversed != contrast->reversed, "contrast: source order changes mapping, not identity");
                        std::array<std::string, 2> haplotypes = {
                            seq.substr(0, offset) + alleles[a] + seq.substr(offset + len),
                            seq.substr(0, offset) + alleles[b] + seq.substr(offset + len)};
                        const auto complete = [&](const std::optional<CandidateIdentityKey>& allele) {
                            if (!allele) return seq;
                            const hts_pos_t pos = allele->pos + (allele->type != static_cast<int>(VariantType::Snp));
                            return seq.substr(0, pos - 1) + allele->alt + seq.substr(pos - 1 + allele->ref_len);
                        };
                        ok &= check(complete(contrast->key.alleles[0]) == haplotypes[contrast->reversed ? 1 : 0] &&
                            complete(contrast->key.alleles[1]) == haplotypes[contrast->reversed ? 0 : 1],
                            "contrast: both normalized alleles reconstruct their independent full haplotypes");
                        std::sort(haplotypes.begin(), haplotypes.end());
                        const auto [prior, inserted] = pairs.emplace(haplotypes, contrast->key);
                        ok &= check(inserted || prior->second == contrast->key,
                            "contrast: equivalent full haplotype pairs share one identity");
                        ++tested;
                    }
                }
            }
        }
    }
    const std::string ref = "CAAATG";
    const auto base = [&](hts_pos_t pos) { return pos >= 1 && pos <= 6 ? ref[pos - 1] : 'N'; };
    const auto alt_pair = normalize_allele_contrast(2, "A", "AA", "AAA", base);
    const auto shifted = normalize_allele_contrast(3, "A", "AAA", "AA", base);
    const auto ref_pair = normalize_allele_contrast(2, "A", "A", "AA", base);
    ok &= check(alt_pair && shifted && ref_pair && alt_pair->key == shifted->key &&
        !(alt_pair->key == ref_pair->key) && alt_pair->key.alleles[0] && alt_pair->key.alleles[1],
        "contrast: ALT/ALT is retained and cannot collapse into REF/selected-ALT");
    ok &= check(!normalize_allele_contrast(2, "C", "AA", "AAA", base) &&
        !normalize_allele_contrast(2, "A", "AA", "AA", base) &&
        !normalize_allele_contrast(2, "A", "N", "AA", base) &&
        !normalize_allele_contrast(2, "A", "*", "AA", base) &&
        !normalize_allele_contrast(2, "A", "", "AA", base),
        "contrast: invalid REF, identical, unknown, symbolic and empty allele pairs abstain");
    std::cout << "Complete allele-pair fixtures: " << tested << "\n";
    return ok;
}

// Test-only batch entry point lets the external normalization oracle exercise
// the production helper without starting a phasing pipeline.
int main(int argc, char** argv) {
    if (argc == 2 && std::string(argv[1]) == "--contrasts") {
        const bool ok = check_allele_contrasts();
        if (ok) std::cout << "ALL PASS\n";
        return ok ? 0 : 1;
    }
    if (argc == 3 && std::string(argv[1]) == "--normalize") {
        std::ifstream fasta(argv[2]);
        if (!fasta) return 2;
        std::string reference, line;
        while (std::getline(fasta, line)) {
            if (!line.empty() && line[0] != '>') reference += line;
        }
        const auto base = [&](hts_pos_t pos) {
            return pos >= 1 && pos <= static_cast<hts_pos_t>(reference.size())
                ? reference[static_cast<size_t>(pos - 1)] : 'N';
        };
        while (std::getline(std::cin, line)) {
            std::istringstream row(line);
            hts_pos_t pos;
            std::string ref, alt;
            if (!(row >> pos >> ref >> alt)) return 2;
            const auto id = normalize_candidate_identity(pos, ref == "." ? "" : ref,
                alt == "." ? "" : alt, base);
            if (!id) std::cout << ".\n";
            else std::cout << id->pos << '\t' << id->type << '\t' << id->ref_len << '\t'
                           << (id->alt.empty() ? "." : id->alt) << '\n';
        }
        return 0;
    }
    bool ok = true;
    ok &= check_allele_contrasts();
    ok &= check_molecule_evidence();
    ok &= check_normalization();
    {
        const CandidateIdentityKey allele{100, static_cast<int>(VariantType::Insertion), 0, "A"};
        const CandidateIdentityKey other{100, static_cast<int>(VariantType::Insertion), 0, "AA"};
        const CandidateIdentityKey compound{100, static_cast<int>(VariantType::Deletion), 2, "G"};
        const CandidateIdentityKey other_compound{100, static_cast<int>(VariantType::Deletion), 2, "T"};
        for (const bool reverse : {false, true}) {
            CandidateIdentityIndex index;
            index.insert(allele, reverse ? 1 : 0);
            ok &= check(index.find(allele) == std::optional<size_t>(reverse ? 1 : 0),
                        "identity: a unique candidate remains available");
            index.insert(allele, reverse ? 1 : 0);
            ok &= check(index.find(allele).has_value(),
                        "identity: repeating the same row is idempotent");
            index.insert(allele, reverse ? 0 : 1);
            index.insert(allele, reverse ? 1 : 0);
            ok &= check(!index.find(allele),
                        "identity: duplicate rows abstain permanently in either order");
            index.insert(other, 2);
            index.insert(compound, 3);
            index.insert(other_compound, 4);
            ok &= check(index.find(other) == std::optional<size_t>(2) &&
                        index.find(compound) == std::optional<size_t>(3) &&
                        index.find(other_compound) == std::optional<size_t>(4),
                        "identity: distinct insertion and replacement alleles remain separate");
            ok &= check(!index.find({101, allele.type, allele.ref_len, allele.alt}),
                        "identity: a missing candidate is not matched by proximity");
        }
        CandidateIdentityIndex raw, sequence;
        const CandidateIdentityKey graph_key{100, allele.type, 1, ">1>4"};
        raw.insert(graph_key, 0);
        raw.insert(graph_key, 1);
        sequence.insert(allele, 0);
        sequence.insert(other, 1);
        const auto first = match_candidate_identity(raw, sequence, allele);
        const auto second = match_candidate_identity(raw, sequence, other);
        ok &= check(first && second && first->index == 0 && second->index == 1 &&
                    !first->is_raw && !second->is_raw,
                    "identity: shared graph key cannot erase distinct selected alternatives");
        ok &= check(!match_candidate_identity(raw, sequence, graph_key),
                    "identity: ambiguous topology key cannot select the first graph row");
        sequence.insert(allele, 2);
        ok &= check(!match_candidate_identity(raw, sequence, allele),
                    "identity: duplicate sequence descriptions cannot select the first row");
        raw.insert(allele, 0);
        ok &= check(!match_candidate_identity(raw, sequence, allele),
                    "identity: unique raw entry cannot bypass ambiguous sequence identity");
        raw.insert(other, 1);
        const auto agreeing = match_candidate_identity(raw, sequence, other);
        ok &= check(agreeing && agreeing->index == 1 && agreeing->is_raw,
                    "identity: agreeing raw and sequence indexes retain a raw match");
        CandidateIdentityIndex conflicting;
        conflicting.insert(other, 3);
        ok &= check(!match_candidate_identity(raw, conflicting, other),
                    "identity: disagreeing unique indexes cannot borrow another allele gauge");
    }
    if (ok) std::cout << "ALL PASS\n";
    return ok ? 0 : 1;
}
