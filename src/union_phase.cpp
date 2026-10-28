// Union gap phasing (--bam): the catalog's sites are phased
// together with the alignment's sample-specific heterozygotes -- clean calls
// and MSA-verified noisy calls the catalog lacks -- by EM over the whole
// read x site matrix, plus local haplotype windows over noisy loci.
// See docs/IMPLEMENTATION.md, "Union gap phasing".

#include "union_internal.hpp"
#include "noise_filter.hpp"
#include "trace.hpp"

namespace pgphase_collect {

// The learned error up to which a catalog repeat indel its allele window phases is written (union EM pass 2).
constexpr double kWindowWriteMaxError = 0.3;
// The median edit distance of each side's reads to its allele's window sequence for that write.
constexpr int kWindowWriteMaxMedianEdits = 1;

std::vector<char> bridge_weak_sites(const GraphChunkBuildResult& graph_chunk, const PhasingChunk& bam) {
    const PhasingChunk& chunk = graph_chunk.chunk;
    std::vector<char> weak(chunk.candidates.size(), 0);
    const bool have_ref = !bam.ref_seq.empty();
    const ReferenceView ref{bam};
    const bool have_meta = graph_chunk.site_meta.size() == chunk.candidates.size() &&
                           graph_chunk.site_allele_orig_idx.size() == chunk.candidates.size();
    for (size_t ci = 0; ci < chunk.candidates.size(); ++ci) {
        const CandidateVariant& c = chunk.candidates[ci];
        if ((c.lcd_var_i_to_cate & kCandGermlineVarCate) == 0) continue;
        // A call from inside a noisy region too long for MSA is unverified:
        // at chr18:20.78 Mb such calls are the real hets of a 98 kb region,
        // at chr20:49.65 Mb 695 of them in 24 kb are a divergent copy's
        // mismatches that flipped a block of 110 reads.
        if (c.unresolved_noisy || c.low_mapq_only) {
            weak[ci] = 1;
            continue;
        }
        if (c.key.type == VariantType::Snp) continue;
        // One homopolymer base between the two alleles.
        if (!c.msa_insertion_alts.empty()) {
            for (size_t a = 0; a < c.msa_insertion_alts.size(); ++a)
                for (size_t b = a + 1; b < c.msa_insertion_alts.size(); ++b)
                    if (std::abs(static_cast<long>(c.msa_insertion_alts[a].size()) -
                                 static_cast<long>(c.msa_insertion_alts[b].size())) == 1) weak[ci] = 1;
            continue;
        }
        if (c.bam_injected && have_ref) {
            const hts_pos_t p0 = c.key.pos - 1;
            if (c.key.type == VariantType::Insertion && c.key.ref_len == 0 && c.key.alt.size() == 1) {
                const char b = static_cast<char>(std::toupper(static_cast<unsigned char>(c.key.alt[0])));
                if (ref.base(p0 - 1) == b || ref.base(p0) == b) weak[ci] = 1;
            } else if (c.key.type == VariantType::Deletion && c.key.ref_len == 1) {
                const char b = ref.base(p0);
                if (ref.base(p0 - 1) == b || ref.base(p0 + 1) == b) weak[ci] = 1;
            }
            continue;
        }
        if (have_meta && !c.bam_injected) {
            const GraphSiteMeta& meta = graph_chunk.site_meta[ci];
            const std::vector<int>& orig = graph_chunk.site_allele_orig_idx[ci];
            // The row's two alleles (orig[0] is REF unless it is a 1/2 row).
            const auto allele_len = [&](int a) -> long {
                return a == 0 ? static_cast<long>(meta.ref.size())
                              : static_cast<long>(meta.alts[static_cast<size_t>(a) - 1].size());
            };
            if (orig.size() < 2 || orig[0] < 0 || orig[1] < 1 || static_cast<size_t>(orig[0]) > meta.alts.size() ||
                static_cast<size_t>(orig[1]) > meta.alts.size()) continue;
            if (std::abs(allele_len(orig[1]) - allele_len(orig[0])) == 1) weak[ci] = 1;
        }
    }
    return weak;
}

void phase_chunk_with_alignment_sites(GraphChunkBuildResult& graph_chunk, const PhasingChunk* bam) {
    PhasingChunk& chunk = graph_chunk.chunk;
    // Reads the graph chunk dropped for lack of an informative catalog allele
    // can still carry a variant only the windows see (a two-allele indel).
    if (bam != nullptr) {
        const std::vector<char> added = add_alignment_only_reads(chunk, *bam);
        // They carry no catalog calls; give them calls by realignment.
        fill_missing_observations(graph_chunk, *bam, added);
    }
    // Pileup calls at repeat indels are unreliable; the alignment's indel sites
    // are re-called per read by realignment to their two exact alleles.
    ObservationErrors errors;
    if (bam != nullptr) realign_indel_observations(chunk, *bam, &errors);
    // The graph's call at a catalog SNP, unbacked by the alignment's base call,
    // is the noisier one (see confirm_catalog_snp_calls).
    if (bam != nullptr) confirm_catalog_snp_calls(graph_chunk, *bam, &errors);
    // The first solve labels reads where the clean stage could not (dense
    // noisy regions); those labels seed the windows the second solve adds.
    const std::vector<char> weak = bam != nullptr ? bridge_weak_sites(graph_chunk, *bam) : std::vector<char>{};
    const std::vector<char>* weak_ptr = bam != nullptr ? &weak : nullptr;
    phase_chunk_by_global_em(chunk, nullptr, weak_ptr, &errors);
    std::vector<LocusWindowSite> loci;
    if (bam != nullptr) {
        loci = build_locus_window_sites(*bam, chunk);
        std::vector<LocusWindowSite> alleles = build_allele_window_sites(*bam, graph_chunk, loci);
        loci.insert(loci.end(), std::make_move_iterator(alleles.begin()), std::make_move_iterator(alleles.end()));
        // Length windows over repeat loci no other site covers, as island-only
        // sites: a last resort for reads the solve leaves unlabelled. As voting
        // sites they cost reads (-99 correct / +22 discordant over
        // chr12/18/20): most split their reads near 50/50 by parent.
        std::vector<LocusWindowSite> lengths = build_length_window_sites(*bam, graph_chunk, loci);
        for (LocusWindowSite& w : lengths) w.island_only = true;
        loci.insert(loci.end(), std::make_move_iterator(lengths.begin()), std::make_move_iterator(lengths.end()));
    }
    phase_chunk_by_global_em(chunk, &loci, weak_ptr, &errors);
    // A catalog repeat indel is never an EM site, so never written; its allele window, the reads realigned to the
    // row's two exact alleles, is. A window the solve phases reliably gives the row its genotype and block.
    for (const LocusWindowSite& w : loci) {
        if (w.catalog_row < 0 || static_cast<size_t>(w.catalog_row) >= chunk.candidates.size()) continue;
        if (w.phase < 0 || w.block <= 0 || w.error >= kWindowWriteMaxError) continue;
        // Each side's reads carry its allele, not merely one nearer to it than the other: at chr20:20,393,257 the
        // reads inserting TGGATGGA sat 4 edits from the catalog's TGGA, closer than to REF, and TGGA was written.
        if (w.median_edits[0] > kWindowWriteMaxMedianEdits || w.median_edits[1] > kWindowWriteMaxMedianEdits) continue;
        CandidateVariant& c = chunk.candidates[static_cast<size_t>(w.catalog_row)];
        if (c.counts.category != VariantCategory::RepeatHetIndel) continue;
        // A simple insertion or deletion only: a catalog walk (a long replacement) the window phases as a whole
        // contradicted the SNPs written inside it (chr20:20,415,921, 41 bases over a SNP on its haplotype).
        if (static_cast<size_t>(w.catalog_row) >= graph_chunk.site_meta.size() ||
            static_cast<size_t>(w.catalog_row) >= graph_chunk.site_allele_orig_idx.size()) continue;
        const GraphSiteMeta& meta = graph_chunk.site_meta[static_cast<size_t>(w.catalog_row)];
        const std::vector<int>& orig = graph_chunk.site_allele_orig_idx[static_cast<size_t>(w.catalog_row)];
        if (orig.size() < 2 || orig[0] != 0 || orig[1] < 1 || static_cast<size_t>(orig[1]) > meta.alts.size()) continue;
        hts_pos_t vpos = meta.pos;
        std::string vref = meta.ref, valt = meta.alts[static_cast<size_t>(orig[1]) - 1];
        trim_to_minimal_vcf(vpos, vref, valt);
        if (vref.size() == valt.size() || std::min(vref.size(), valt.size()) != 1 || vref[0] != valt[0]) continue;
        c.hap_to_cons_alle[1] = w.phase;
        c.hap_to_cons_alle[2] = 1 - w.phase;
        c.phase_set = w.block;
        c.vcf_unphased = false;
        c.window_phased = true;
        pgphase_trace::candidate("em", "window", c, "repeat_indel_window_phased",
                                 "err=" + std::to_string(w.error) + " block=" + std::to_string(w.block));
    }
}

} // namespace pgphase_collect
