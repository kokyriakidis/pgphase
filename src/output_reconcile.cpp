/**
 * @file output_reconcile.cpp
 * Post-phasing record reconciliation ("output" trace stage).
 */

#include "output_reconcile.hpp"
#include "trace.hpp"

#include <algorithm>
#include <cctype>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace pgphase_collect {

// One haplotype cannot carry two different alleles at one position, and we
// were emitting exactly that: 106 positions on chr20 where a SNP and the
// insertion that CONTAINS it were both called on the same haplotype --
// T>G together with T>GAG at 288,018, T>A with T>AC at 890,261. Applying
// both is incoherent; one of the pair has to go.
//
// The reads decide which, and they say KEEP THE LONGER. Classifying every
// read at five of these loci against both forms: at 288,018, 40 reads carry
// the reference and 24 carry GAG, with none carrying a bare G; at 890,261,
// 32 carry AC against 25 reference and none carry a bare A; at 764,862 and
// 1,474,295 the short form appears on 6 of 54 and 14 of 79 reads, all on the
// same parent as the long form, i.e. the same event aligned two ways. The
// alt-carrying reads carry the whole insertion, so dropping it and keeping
// the SNP would report an allele no read has and lose the inserted bases.
//
// Allele depth cannot arbitrate here -- both records score the same reads,
// 38,21 against 29,18 -- and neither can the alignment arm, which emits the
// SNP alone at these loci (at 890,261 a T>A plus a separate T>C at 890,262);
// its representation splits the event across records rather than showing
// that the bare SNP is what the haplotype carries. So the rule is
// containment resolved towards the complete allele: when one ALT is a prefix
// of the other on a haplotype already claimed at that position, the
// contained (shorter) record is dropped.
/// Make two co-located records complementary instead of discarding one.
///
/// longcallD cannot express a haplotype carrying two alleles at one position:
/// it writes one ALT per candidate and takes the genotype from that
/// candidate's own haplotype alleles, so its two rows at a locus are
/// complementary by construction, and it scores ZERO such positions on chr20.
/// Ours scored 412 with the multiallelic merge off -- two nested MSA
/// insertions, A>ATC and A>ATCTC at 882,277, both called 1|0 with 11 and 15
/// alternate reads against 2 reference.
///
/// Both alleles have real support there, so the locus is heterozygous with two
/// non-reference alleles and the right answer is one allele per haplotype --
/// which is what upstream emits. Dropping the contained form instead costs 798
/// records on chr20 and takes two-row positions to 1,542 against upstream's
/// 2,401.
///
/// So: when two records at one position claim the SAME haplotype and both
/// alleles carry alternate support, move the less supported one to the other
/// haplotype. Returns the number of records re-oriented. Anything without
/// support on both sides is left to drop_conflicting_haplotype_alleles.
// Two equal-length replacements starting at one position, the shorter a
// prefix of the longer: the longer allele makes the shorter one's change too,
// so the two are one haplotype's event described twice, not two alleles. A
// catalog snarl spanning several SNPs gives exactly this -- at chr12:1,069,108
// the 280 bp record (three SNPs) and the SNP T>C at its first base, both on
// the C-carrying haplotype, which the rule below then split across the two
// haplotypes: the SNP was written opposite to its own reads' tags. 38 of the
// 64 single-site orientation flips on chr12's arms were this. Both records
// stay: the SNP is the minimal form (the one the alignment arm and truth sets
// use), and the long record carries the rest of the snarl. Dropping the SNP
// instead cost 101 truth heterozygotes on chr12, 49 on chr18 and 74 on chr20;
// the same pair already sits side by side wherever the SNP is not the snarl's
// first base (1,069,183 in the record above).
static bool nested_substitutions(const CandidateVariant& a, const CandidateVariant& b) {
    const auto substitution = [](const CandidateVariant& c) {
        return c.key.type == VariantType::Snp &&
               c.key.ref_len == static_cast<int>(c.key.alt.size());
    };
    if (!substitution(a) || !substitution(b) || a.key.pos != b.key.pos) return false;
    const bool a_short = a.key.alt.size() <= b.key.alt.size();
    const std::string& shorter = a_short ? a.key.alt : b.key.alt;
    const std::string& longer = a_short ? b.key.alt : a.key.alt;
    return longer.compare(0, shorter.size(), shorter) == 0;
}

// Records by the POS they are written at (vcf_record_pos): the two rules below compare alleles that start at one
// reference base. Grouped by key.sort_pos() instead, a replacement (written unpadded at key.pos) sat with a SNP on
// the base before it, and the SNP, its ALT a prefix of the replacement's, was dropped as the contained form: chr20:
// 6,670,080 T>G beside T>GAGGAGGA at 6,670,081, chr12:18,697,237 G>C beside G>CT at 238, each the MSA's haplotype.
static std::map<std::pair<int, hts_pos_t>, std::vector<size_t>> records_by_vcf_pos(const CandidateTable& result) {
    std::map<std::pair<int, hts_pos_t>, std::vector<size_t>> by_pos;
    for (size_t i = 0; i < result.size(); ++i) by_pos[{result[i].key.tid, vcf_record_pos(result[i])}].push_back(i);
    return by_pos;
}

// The MSA's consensus-to-reference alignment can describe a replacement as a deletion and an insertion side by side
// (chr18:830,100, TTTG deleted then AAAAAA inserted after it, the catalog's TTTTG>TAAAAAA; chr20:30,110,228). Written
// apart, the insertion's padding base is the deletion's last deleted base: the two records claim one reference base on
// one haplotype, so no reader can apply both, and with the insertion before the deletion the two share a POS and the
// conflict rule below dropped one. Merged only where both are surely on one haplotype: phased in one phase set with
// the same genotype, or both 1/1.
size_t merge_adjacent_deletion_insertion(CandidateTable& result) {
    const auto same_haplotype = [](const CandidateVariant& a, const CandidateVariant& b) {
        const auto& ha = a.hap_to_cons_alle;
        const auto& hb = b.hap_to_cons_alle;
        if (ha[1] != hb[1] || ha[2] != hb[2] || (ha[1] <= 0 && ha[2] <= 0)) return false;
        if (ha[1] > 0 && ha[2] > 0) return true;
        return a.phase_set > 0 && a.phase_set == b.phase_set && !a.vcf_unphased && !b.vcf_unphased;
    };
    const auto plain_deletion = [](const CandidateVariant& c) {
        return c.key.type == VariantType::Deletion && c.key.ref_len > 0 && c.key.alt.empty() &&
               c.msa_insertion_alts.empty();
    };
    const auto plain_insertion = [](const CandidateVariant& c) {
        return c.key.type == VariantType::Insertion && c.key.ref_len == 0 && !c.key.alt.empty() &&
               c.msa_insertion_alts.empty();
    };
    std::map<std::pair<int, hts_pos_t>, std::vector<size_t>> insertions;
    for (size_t i = 0; i < result.size(); ++i)
        if (plain_insertion(result[i])) insertions[{result[i].key.tid, result[i].key.pos}].push_back(i);
    std::vector<char> drop(result.size(), 0);
    size_t merged = 0;
    for (size_t d = 0; d < result.size(); ++d) {
        CandidateVariant& del = result[d];
        if (!plain_deletion(del)) continue;
        // An insertion right after the deleted bases, else right before them.
        for (const hts_pos_t at : {del.key.pos + del.key.ref_len, del.key.pos}) {
            const auto it = insertions.find({del.key.tid, at});
            if (it == insertions.end()) continue;
            size_t found = result.size();
            for (const size_t i : it->second)
                if (!drop[i] && same_haplotype(del, result[i])) {
                    if (found != result.size()) { found = result.size() + 1; break; }  // two candidates: leave them
                    found = i;
                }
            if (found >= result.size()) continue;
            const CandidateVariant& ins = result[found];
            pgphase_trace::candidate("output", "merge", ins, "deletion_insertion_one_haplotype",
                                     "into=" + pgphase_trace::key_text(del.key));
            del.key.alt = ins.key.alt;
            if (del.key.alt.size() > static_cast<size_t>(del.key.ref_len)) del.key.type = VariantType::Insertion;
            drop[found] = 1;
            ++merged;
            break;
        }
    }
    if (merged > 0) {
        size_t write = 0;
        for (size_t i = 0; i < result.size(); ++i)
            if (!drop[i]) {
                if (write != i) result[write] = std::move(result[i]);
                ++write;
            }
        result.resize(write);
    }
    return merged;
}

// A catalog walk -- a deletion or replacement longer than kWalkMinLength, the snarl's path spelled whole -- and a
// record inside it on one haplotype of a phase set describe that haplotype twice; where the inner record changes a base
// the walk deletes or keeps as reference, the two cannot both hold. Against the truth (HG002, chr12/17/18/19/20/21) the
// walk was wrong in 88 of 94 such pairs (the inner record right in 48): the graph's reads are forced onto a walk where
// the sample's haplotype is none of the snarl's (chr19:842,288, a 219-base deletion over two true 120-base insertions).
// The walk is dropped: over chr12/17/18/19/20/21, 11 walks, 16 contradictory overlaps fewer; truth hets unchanged.
static constexpr int kWalkMinLength = 10;

size_t drop_contradicted_walks(CandidateTable& result, const std::function<char(int, hts_pos_t)>& ref_base) {
    const auto phased_hap = [](const CandidateVariant& c, int h) {
        return c.phase_set > 0 && !c.vcf_unphased && c.hap_to_cons_alle.size() > 2 && c.hap_to_cons_alle[h] > 0;
    };
    std::vector<size_t> order(result.size());
    for (size_t i = 0; i < order.size(); ++i) order[i] = i;
    std::sort(order.begin(), order.end(),
              [&](size_t a, size_t b) { return vcf_record_pos(result[a]) < vcf_record_pos(result[b]); });
    std::vector<hts_pos_t> starts;
    for (const size_t i : order) starts.push_back(vcf_record_pos(result[i]));
    std::vector<char> drop(result.size(), 0);
    size_t dropped = 0;
    for (size_t w = 0; w < result.size(); ++w) {
        const CandidateVariant& walk = result[w];
        if (!walk.msa_insertion_alts.empty() || walk.bam_injected || walk.key.type == VariantType::Snp) continue;
        // The bases the walk replaces (from key.pos) and what it writes there (key.alt; its padding base aside).
        const hts_pos_t beg = walk.key.pos;
        const int ref_len = walk.key.ref_len;
        const std::string& alt = walk.key.alt;
        if (ref_len <= kWalkMinLength && static_cast<int>(alt.size()) <= kWalkMinLength) continue;
        if (ref_len < 1) continue;
        std::string ref;
        for (int i = 0; i < ref_len; ++i) ref.push_back(ref_base(walk.key.tid, beg + i));
        // Where an inner record's change contradicts the walk: a base the walk deletes (a pure deletion), or one it
        // keeps -- its common prefix and suffix with REF. Inside a replacement's changed core the two may agree.
        size_t prefix = 0, suffix = 0;
        while (prefix < ref.size() && prefix < alt.size() && ref[prefix] == alt[prefix]) ++prefix;
        while (suffix + prefix < ref.size() && suffix + prefix < alt.size() &&
               ref[ref.size() - 1 - suffix] == alt[alt.size() - 1 - suffix]) ++suffix;
        const auto contradicts_at = [&](size_t off) {
            if (alt.empty()) return true;
            if (ref.size() == alt.size()) return ref[off] == alt[off];
            return off < prefix || off >= ref.size() - suffix;
        };
        for (size_t k = std::lower_bound(starts.begin(), starts.end(), beg) - starts.begin();
             k < order.size() && starts[k] < beg + ref_len; ++k) {
            const size_t i = order[k];
            if (i == w) continue;
            const CandidateVariant& inner = result[i];
            if (inner.phase_set != walk.phase_set) continue;
            bool same_hap = false;
            for (int h = 1; h <= 2; ++h) same_hap = same_hap || (phased_hap(walk, h) && phased_hap(inner, h));
            if (!same_hap || !contradicts_at(static_cast<size_t>(starts[k] - beg))) continue;
            drop[w] = 1;
            break;
        }
        if (drop[w]) {
            ++dropped;
            pgphase_trace::candidate("output", "drop", walk, "walk_contradicted", pgphase_trace::genotype_text(walk));
        }
    }
    if (dropped > 0) {
        size_t write = 0;
        for (size_t i = 0; i < result.size(); ++i)
            if (!drop[i]) {
                if (write != i) result[write] = std::move(result[i]);
                ++write;
            }
        result.resize(write);
    }
    return dropped;
}

size_t make_colocated_alleles_complementary(CandidateTable& result, int min_alt_support) {
    const auto by_pos = records_by_vcf_pos(result);
    // A BAM-discovered candidate carries its depths in ref_cov/alt_cov and
    // leaves alle_covs empty -- 662 such candidates per megabase on chr20 --
    // so reading alle_covs[1] alone returned 0 and the min_alt_support test
    // below rejected every pair on the alignment arm. This rule was dead there.
    const auto alt_depth = [](const CandidateVariant& c) {
        if (c.counts.alle_covs.size() > 1) return c.counts.alle_covs[1];
        return c.counts.alt_cov;
    };
    size_t flipped = 0;
    for (const auto& all : by_pos) {
        // Records with no ALT on either haplotype (an uncalled catalog copy) write nothing: they do not make a pair a
        // crowd (chr20:43,522,037, two uncalled repeat rows beside the MSA's insertion and deletion).
        std::vector<size_t> group;
        for (const size_t i : all.second) {
            const auto& h = result[i].hap_to_cons_alle;
            if (h.size() < 3 || h[1] > 0 || h[2] > 0) group.push_back(i);
        }
        if (group.size() != 2) continue;
        CandidateVariant& a = result[group[0]];
        CandidateVariant& b = result[group[1]];
        if (a.hap_to_cons_alle.size() < 3 || b.hap_to_cons_alle.size() < 3) continue;
        if (a.key.alt == b.key.alt) continue;         // same allele twice: not this rule
        if (nested_substitutions(a, b)) continue;     // one event described twice
        if (alt_depth(a) < min_alt_support || alt_depth(b) < min_alt_support) continue;
        // Two records the phaser put in one phase set were each phased from
        // their own reads, so where both claim one haplotype the reads say so:
        // they are one allele described twice (a catalog walk and the
        // alignment's call, two stutter lengths of one repeat insertion), not
        // two alleles. Moving one to the other haplotype wrote it on the
        // reference haplotype -- at chr12:13,569,119 +12 and +14 bp AT, both on
        // haplotype 2 at learned error 0.01, the +14 (the truth's allele) came
        // out 1|0: 5 of 36 arm flips over chr12/18/20. The pair is left to
        // drop_conflicting_haplotype_alleles, which keeps one description.
        if (a.phase_set > 0 && a.phase_set == b.phase_set) continue;
        // Records in two different phase sets: their haplotype numbers are
        // unrelated, so "both claim haplotype 1" says nothing. Moving one
        // inverted it against its own reads -- chr20:7,184,086 +4 bp AC
        // (phase set 7,184,047), carried by 38 paternal reads on haplotype 2,
        // was moved to haplotype 1 because a +2 bp record in phase set
        // 7,184,085 claimed haplotype 2 too. Arm switches -1 on chr20 and
        // chr12, flips -1 on chr12 and chr17; reads unchanged.
        if (a.phase_set > 0 && b.phase_set > 0 && a.phase_set != b.phase_set) continue;
        int hap = 0;
        for (int h = 1; h <= 2; ++h)
            if (a.hap_to_cons_alle[h] == 1 && b.hap_to_cons_alle[h] == 1) hap = h;
        if (hap == 0) {
            // The other way a co-located pair fails to be complementary: one
            // candidate claims BOTH haplotypes and the other claims NEITHER, so
            // the first is emitted 1|1 and the second dropped for carrying no
            // ALT. longcallD splits these -- at chr20:3,863,176 it writes
            // C>AAAAAAAAA 1|0 beside C>AAAAAAAAAAA 0|1 where we wrote 1|1 and
            // nothing. Give each allele the haplotype whose own profile prefers
            // it, and only when both profiles are present and disagree, so the
            // orientation comes from read evidence rather than from the order
            // the candidates happen to sit in.
            const bool a_both = a.hap_to_cons_alle[1] == 1 && a.hap_to_cons_alle[2] == 1;
            const bool b_none = b.hap_to_cons_alle[1] != 1 && b.hap_to_cons_alle[2] != 1;
            const bool b_both = b.hap_to_cons_alle[1] == 1 && b.hap_to_cons_alle[2] == 1;
            const bool a_none = a.hap_to_cons_alle[1] != 1 && a.hap_to_cons_alle[2] != 1;
            if (!((a_both && b_none) || (b_both && a_none))) continue;
            CandidateVariant& both = a_both ? a : b;
            CandidateVariant& none = a_both ? b : a;
            const auto alt_prof = [](const CandidateVariant& c, int h) -> int {
                if (c.hap_to_alle_profile.size() < 3) return -1;
                const auto& row = c.hap_to_alle_profile[static_cast<size_t>(h)];
                return row.size() > 1 ? row[1] : -1;
            };
            const int n1 = alt_prof(none, 1), n2 = alt_prof(none, 2);
            if (n1 < 0 || n2 < 0 || n1 == n2) continue;
            const int none_hap = n1 > n2 ? 1 : 2;       // where the dropped allele is stronger
            const int both_hap = none_hap == 1 ? 2 : 1;
            none.hap_to_cons_alle[none_hap] = 1;
            none.hap_to_cons_alle[both_hap] = 0;
            both.hap_to_cons_alle[both_hap] = 1;
            both.hap_to_cons_alle[none_hap] = 0;
            // Keep hap_alt/hap_ref in step with the consensus, the way the
            // rule below this one does with its swap: collect_phase.cpp:1132-1136
            // encodes both-alt as 3/0, hap1-alt as 1/2 and hap2-alt as 2/1.
            const auto sync_hap_fields = [](CandidateVariant& c) {
                const bool h1 = c.hap_to_cons_alle[1] != 0;
                const bool h2 = c.hap_to_cons_alle[2] != 0;
                c.hap_alt = (h1 && h2) ? 3 : (h1 ? 1 : (h2 ? 2 : 0));
                c.hap_ref = (h1 && h2) ? 0 : (h1 ? 2 : (h2 ? 1 : 3));
            };
            sync_hap_fields(none);
            sync_hap_fields(both);
            pgphase_trace::candidate("output", "regenotype", none, "colocated_split_both_none",
                                     pgphase_trace::genotype_text(none));
            pgphase_trace::candidate("output", "regenotype", both, "colocated_split_both_none",
                                     pgphase_trace::genotype_text(both));
            ++flipped;
            continue;
        }
        // A record written unphased says nothing of its haplotype; the phased one's comes from its reads. Move the
        // unphased one, whatever their depths.
        const auto phased = [](const CandidateVariant& c) { return c.phase_set > 0 && !c.vcf_unphased; };
        CandidateVariant& weaker = phased(a) != phased(b) ? (phased(a) ? b : a)
                                 : alt_depth(a) <= alt_depth(b) ? a : b;
        const int other = hap == 1 ? 2 : 1;
        pgphase_trace::candidate("output", "regenotype", weaker, "colocated_moved_to_other_hap",
                                 pgphase_trace::genotype_text(weaker) + " hap=" + std::to_string(hap) + " to_hap=" +
                                 std::to_string(other));
        weaker.hap_to_cons_alle[hap] = 0;
        weaker.hap_to_cons_alle[other] = 1;
        std::swap(weaker.hap_alt, weaker.hap_ref);
        ++flipped;
    }
    return flipped;
}

void drop_conflicting_haplotype_alleles(CandidateTable& result) {
    const auto by_pos = records_by_vcf_pos(result);
    std::vector<bool> drop(result.size(), false);
    for (const auto& group : by_pos) {
        if (group.second.size() < 2) continue;
        for (int hap = 1; hap <= 2; ++hap) {
            std::vector<size_t> carriers;
            for (size_t i : group.second) {
                const auto& h = result[i].hap_to_cons_alle;
                if (h.size() > static_cast<size_t>(hap) && h[hap] == 1) carriers.push_back(i);
            }
            // A substitution nested in a longer carrier is the same event
            // (nested_substitutions) and is no conflict: it stays.
            const std::vector<size_t> all_carriers = carriers;
            carriers.erase(std::remove_if(carriers.begin(), carriers.end(), [&](size_t i) {
                return std::any_of(all_carriers.begin(), all_carriers.end(), [&](size_t j) {
                    return result[i].key.alt.size() < result[j].key.alt.size() &&
                           nested_substitutions(result[i], result[j]);
                });
            }), carriers.end());
            if (carriers.size() < 2) continue;
            std::sort(carriers.begin(), carriers.end(), [&](size_t a, size_t b) {
                return result[a].key.alt.size() < result[b].key.alt.size();
            });
            // carriers is sorted shortest ALT first; a shorter allele that the
            // longest one contains is the contained form and is dropped.
            const std::string& longest = result[carriers.back()].key.alt;
            bool resolved = false;
            // A deletion's ALT is empty, a prefix of every allele: containment holds only between alleles of one type
            // (chr20:43,522,037, the MSA's 3 bp deletion dropped as "contained" in the insertion beside it).
            const VariantType longest_type = result[carriers.back()].key.type;
            for (size_t k = 0; k + 1 < carriers.size(); ++k) {
                const std::string& other = result[carriers[k]].key.alt;
                if (!other.empty() && result[carriers[k]].key.type == longest_type && other.size() < longest.size() &&
                    longest.compare(0, other.size(), other) == 0) {
                    if (!drop[carriers[k]])
                        pgphase_trace::candidate("output", "drop", result[carriers[k]], "conflict_contained",
                                                 pgphase_trace::genotype_text(result[carriers[k]]) + " hap=" +
                                                 std::to_string(hap) + " in=" +
                                                 pgphase_trace::key_text(result[carriers.back()].key));
                    drop[carriers[k]] = true;
                    resolved = true;
                }
            }
            if (resolved) continue;
            // Neither allele contains the other -- an insertion anchored on the
            // reference base together with a SNP changing that same base, which
            // is equally incoherent on one haplotype and cannot be settled by
            // containment. Here the two records do NOT score the same reads, so
            // allele depth separates them (at 1,907,006 the insertion carries 10
            // against the SNP's 5); keep the better supported one, breaking a
            // tie on the shorter allele so the choice is deterministic.
            auto alt_depth = [](const CandidateVariant& c) {
                return c.counts.alle_covs.size() > 1 ? c.counts.alle_covs[1] : 0;
            };
            size_t best = carriers.front();
            for (size_t i : carriers) {
                if (alt_depth(result[i]) > alt_depth(result[best]) ||
                    (alt_depth(result[i]) == alt_depth(result[best]) &&
                     result[i].key.alt.size() < result[best].key.alt.size()))
                    best = i;
            }
            for (size_t i : carriers) {
                if (i == best) continue;
                if (!drop[i])
                    pgphase_trace::candidate("output", "drop", result[i], "conflict_weaker",
                                             pgphase_trace::genotype_text(result[i]) + " hap=" + std::to_string(hap) +
                                             " kept=" + pgphase_trace::key_text(result[best].key));
                drop[i] = true;
            }
        }
    }
    CandidateTable kept;
    kept.reserve(result.size());
    for (size_t i = 0; i < result.size(); ++i)
        if (!drop[i]) kept.push_back(std::move(result[i]));
    result = std::move(kept);
}


} // namespace pgphase_collect
