// Union gap phasing (--bam): the catalog's sites are phased
// together with the alignment's sample-specific heterozygotes -- clean calls
// and MSA-verified noisy calls the catalog lacks -- by EM over the whole
// read x site matrix, plus local haplotype windows over noisy loci.
// See docs/IMPLEMENTATION.md, "Union gap phasing".

#include "union_internal.hpp"

namespace pgphase_collect {

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
            if (orig.size() < 2 || orig[1] < 1 || static_cast<size_t>(orig[1]) > meta.alts.size()) continue;
            const long sel = static_cast<long>(meta.alts[static_cast<size_t>(orig[1]) - 1].size());
            if (std::abs(sel - static_cast<long>(meta.ref.size())) == 1 && !meta.non_selected_alt_class) weak[ci] = 1;
            if (meta.non_selected_alt_class)
                for (size_t a = 0; a < meta.alts.size(); ++a)
                    if (static_cast<int>(a) + 1 != orig[1] && std::abs(sel - static_cast<long>(meta.alts[a].size())) == 1)
                        weak[ci] = 1;
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
    if (bam != nullptr) realign_indel_observations(chunk, *bam);
    // The first solve labels reads where the clean stage could not (dense
    // noisy regions); those labels seed the windows the second solve adds.
    const std::vector<char> weak = bam != nullptr ? bridge_weak_sites(graph_chunk, *bam) : std::vector<char>{};
    const std::vector<char>* weak_ptr = bam != nullptr ? &weak : nullptr;
    phase_chunk_by_global_em(chunk, nullptr, weak_ptr);
    std::vector<LocusWindowSite> loci;
    if (bam != nullptr) {
        loci = build_locus_window_sites(*bam, chunk);
        std::vector<LocusWindowSite> alleles = build_allele_window_sites(*bam, graph_chunk, loci);
        loci.insert(loci.end(), std::make_move_iterator(alleles.begin()), std::make_move_iterator(alleles.end()));
    }
    phase_chunk_by_global_em(chunk, &loci, weak_ptr);
    // Last resort: reads still unlabelled take their haplotype from MSA-verified
    // indels whose phase the labelled reads already establish.
    if (bam != nullptr) {
        // Verified indels first (exact alleles), then the labelled reads' own
        // haplotype markers for reads still undecided.
        label_reads_from_verified_indels(graph_chunk, *bam);
        label_reads_from_haplotype_consensus(graph_chunk, *bam);
    }
}

} // namespace pgphase_collect
