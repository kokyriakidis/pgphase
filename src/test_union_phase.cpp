// Unit tests for union gap phasing: site injection, global EM, local windows.

#include "union_phase.hpp"

#include "collect_phase.hpp"
#include "collect_phase_noisy.hpp"

#include <cstdint>
#include <iostream>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include <htslib/sam.h>

using namespace pgphase_collect;

static int failures = 0;

static void check(bool cond, const std::string& msg) {
    if (!cond) {
        std::cerr << "FAIL: " << msg << "\n";
        ++failures;
    }
}

static CandidateVariant snp(hts_pos_t pos, VariantCategory category, uint32_t cate) {
    CandidateVariant c;
    c.key.tid = 0;
    c.key.pos = pos;
    c.key.type = VariantType::Snp;
    c.key.ref_len = 1;
    c.key.alt = "T";
    c.counts.category = category;
    c.lcd_var_i_to_cate = cate;
    return c;
}

static CandidateVariant clean_snp(hts_pos_t pos) {
    return snp(pos, VariantCategory::CleanHetSnp, kCandCleanHetSnp);
}

static ReadRecord read_named(const std::string& name, int mapq) {
    ReadRecord r;
    r.tid = 0;
    r.qname = name;
    r.mapq = mapq;
    return r;
}

static ReadVariantProfile profile(int read_id, int start, std::vector<int> alleles) {
    ReadVariantProfile p;
    p.read_id = read_id;
    p.start_var_idx = start;
    p.end_var_idx = start + static_cast<int>(alleles.size()) - 1;
    p.alt_qi.assign(alleles.size(), 0);
    p.alleles = std::move(alleles);
    return p;
}

static void finish_chunk(PhasingChunk& chunk) {
    chunk.region.tid = 0;
    chunk.haps.assign(chunk.reads.size(), 0);
    chunk.phase_sets.assign(chunk.reads.size(), kUnphasedReadPhaseSet);
}

// Haplotype A carries allele (k % 2) at site k, haplotype B the other one.
static int hap_allele(int hap, size_t site) { return static_cast<int>((site + static_cast<size_t>(hap)) % 2); }

// Reads tiled over sites [first, last], each covering `span` sites, haplotypes alternating.
static void add_tiled_reads(PhasingChunk& chunk, std::vector<int>& truth, size_t first, size_t last,
                            size_t span, int mapq = 60) {
    for (size_t s = first; s + span - 1 <= last; ++s) {
        for (int hap = 0; hap < 2; ++hap) {
            for (int copy = 0; copy < 2; ++copy) {
                std::vector<int> alleles;
                for (size_t k = s; k < s + span; ++k) alleles.push_back(hap_allele(hap, k));
                const int id = static_cast<int>(chunk.reads.size());
                chunk.reads.push_back(read_named("r" + std::to_string(id), mapq));
                chunk.read_var_profile.push_back(profile(id, static_cast<int>(s), alleles));
                truth.push_back(hap);
            }
        }
    }
}

// Every labelled read agrees with the truth under one global orientation.
static bool reads_consistent(const PhasingChunk& chunk, const std::vector<int>& truth, size_t* labelled) {
    int orientation = 0;
    *labelled = 0;
    for (size_t r = 0; r < truth.size(); ++r) {
        if (chunk.haps[r] == 0) continue;
        ++*labelled;
        const int o = (chunk.haps[r] == 1) == (truth[r] == 0) ? 1 : -1;
        if (orientation == 0) orientation = o;
        if (o != orientation) return false;
    }
    return true;
}

static void test_em_phases_one_block() {
    PhasingChunk chunk;
    for (int k = 0; k < 12; ++k) chunk.candidates.push_back(clean_snp(100 + 100 * k));
    std::vector<int> truth;
    add_tiled_reads(chunk, truth, 0, 11, 4);
    finish_chunk(chunk);
    const size_t blocks = phase_chunk_by_global_em(chunk);
    check(blocks == 1, "EM: tiled reads give one block, got " + std::to_string(blocks));
    std::set<hts_pos_t> sets;
    int orientation = 0;
    bool sites_ok = true;
    for (size_t k = 0; k < chunk.candidates.size(); ++k) {
        const CandidateVariant& c = chunk.candidates[k];
        sets.insert(c.phase_set);
        const int o = c.hap_to_cons_alle[1] == hap_allele(0, k) ? 1 : -1;
        if (orientation == 0) orientation = o;
        sites_ok = sites_ok && o == orientation && c.hap_to_cons_alle[2] == 1 - c.hap_to_cons_alle[1];
    }
    check(sets.size() == 1 && *sets.begin() > 0, "EM: all sites share one phase set");
    check(sites_ok, "EM: site phases consistent with the truth");
    size_t labelled = 0;
    check(reads_consistent(chunk, truth, &labelled), "EM: read labels consistent with the truth");
    check(labelled == truth.size(), "EM: every read labelled");
}

static void test_em_cuts_unlinked_blocks() {
    PhasingChunk chunk;
    for (int k = 0; k < 10; ++k) chunk.candidates.push_back(clean_snp(100 + 100 * k));
    std::vector<int> truth;
    add_tiled_reads(chunk, truth, 0, 4, 3);
    add_tiled_reads(chunk, truth, 5, 9, 3);
    finish_chunk(chunk);
    const size_t blocks = phase_chunk_by_global_em(chunk);
    check(blocks == 2, "EM: no read spans sites 4-5, expected 2 blocks, got " + std::to_string(blocks));
    check(chunk.candidates[4].phase_set != chunk.candidates[5].phase_set, "EM: block cut between 4 and 5");
}

static void test_em_leaves_noise_unphased() {
    PhasingChunk chunk;
    for (int k = 0; k < 8; ++k) chunk.candidates.push_back(clean_snp(100 + 100 * k));
    chunk.candidates.push_back(clean_snp(850));  // index 8: alleles unrelated to haplotype
    std::vector<int> truth;
    // 32 reads: a sites learned error leans on its classs pooled error until
    // its own reads outweigh the prior (kEmClassPriorReads), so the property
    // under test -- a noise site at ordinary depth is rejected -- needs depth.
    for (int rep = 0; rep < 8; ++rep) add_tiled_reads(chunk, truth, 0, 7, 8);
    for (size_t r = 0; r < chunk.read_var_profile.size(); ++r) {
        ReadVariantProfile& p = chunk.read_var_profile[r];
        p.alleles.push_back(static_cast<int>((r / 2) % 2) ^ static_cast<int>(r % 2) ^ truth[r]);
        p.alt_qi.push_back(0);
        p.end_var_idx = 8;
    }
    // The noisy column above is a function of read index only; make it balanced.
    finish_chunk(chunk);
    phase_chunk_by_global_em(chunk);
    check(chunk.candidates[8].hap_to_cons_alle[1] == -1, "EM: unreliable site left unphased");
    check(chunk.candidates[0].hap_to_cons_alle[1] >= 0, "EM: clean sites stay phased");
}

static void test_em_recovers_from_switched_start() {
    PhasingChunk chunk;
    for (int k = 0; k < 12; ++k) chunk.candidates.push_back(clean_snp(100 + 100 * k));
    std::vector<int> truth;
    add_tiled_reads(chunk, truth, 0, 11, 4);
    finish_chunk(chunk);
    // Start from a switch: reads starting right of site 6 begin flipped.
    for (size_t r = 0; r < truth.size(); ++r) {
        const bool right = chunk.read_var_profile[r].start_var_idx > 6;
        chunk.haps[r] = (truth[r] == 0) != right ? 1 : 2;
    }
    const size_t blocks = phase_chunk_by_global_em(chunk);
    size_t labelled = 0;
    check(blocks == 1, "EM: switched start still gives one block");
    check(reads_consistent(chunk, truth, &labelled) && labelled == truth.size(),
          "EM: switched start resolved");
}

static void test_em_labels_low_mapq_reads() {
    PhasingChunk chunk;
    for (int k = 0; k < 8; ++k) chunk.candidates.push_back(clean_snp(100 + 100 * k));
    std::vector<int> truth;
    add_tiled_reads(chunk, truth, 0, 7, 4);
    add_tiled_reads(chunk, truth, 2, 5, 4, 5);  // MAPQ 5: labelled, not learned from
    finish_chunk(chunk);
    phase_chunk_by_global_em(chunk);
    size_t labelled = 0;
    check(reads_consistent(chunk, truth, &labelled) && labelled == truth.size(),
          "EM: low-MAPQ reads labelled consistently");
}

static void test_em_window_links_blocks() {
    // Sites 0-3 and 4-7; reads on the left also see a window locus, reads on
    // the right see it too. Only the window connects the two sides.
    std::vector<int> truth;
    const auto build = [&truth](PhasingChunk& chunk) {
        truth.clear();
        for (int k = 0; k < 8; ++k) chunk.candidates.push_back(clean_snp(100 + 100 * k));
        add_tiled_reads(chunk, truth, 0, 3, 4);
        add_tiled_reads(chunk, truth, 4, 7, 4);
        finish_chunk(chunk);
    };
    PhasingChunk without;
    build(without);
    check(phase_chunk_by_global_em(without) == 2, "EM: two blocks without the window");
    PhasingChunk chunk;
    build(chunk);
    std::vector<LocusWindowSite> loci(1);
    loci[0].pos = 450;
    for (size_t r = 0; r < truth.size(); ++r)
        loci[0].observations.emplace_back(r, truth[r], 0.01f);
    check(phase_chunk_by_global_em(chunk, &loci) == 1, "EM: the window joins the blocks");
    size_t labelled = 0;
    check(reads_consistent(chunk, truth, &labelled) && labelled == truth.size(),
          "EM: window join is consistent with the truth");
}

static void test_inject_selects_private_sites() {
    std::string ref;
    for (int i = 0; i < 500; ++i) ref += "ACGT"[(i * 7 + i / 3) % 4];
    ref.replace(399, 7, "CAAAAAG");  // 1-based 400: C, 401-405: A run, 406: G
    for (size_t i : {99, 199, 299, 349}) ref[i] = 'A';  // SNP REF differs from ALT T

    GraphChunkBuildResult gc;
    PhasingChunk& graph = gc.chunk;
    graph.ref_beg = 1;
    graph.ref_end = 500;
    graph.candidates.push_back(clean_snp(100));
    for (int r = 0; r < 6; ++r) {
        graph.reads.push_back(read_named("r" + std::to_string(r), 60));
        graph.read_var_profile.push_back(profile(r, 0, {r % 2}));
    }
    gc.site_ids = {"s100"};
    gc.site_meta.resize(1);
    gc.site_allele_orig_idx = {{0, 1}};
    finish_chunk(graph);

    PhasingChunk bam;
    bam.ref_beg = 1;
    bam.ref_end = 500;
    bam.ref_seq = ref;
    bam.candidates.push_back(clean_snp(100));                                    // 0: graph has it
    bam.candidates.push_back(clean_snp(200));                                    // 1: injected
    // 2: ALT only at low MAPQ, outside the het categories (a clean het read so
    // is injected: test_inject_low_mapq_clean_het)
    bam.candidates.push_back(snp(300, VariantCategory::NoisyCandHom, 0));
    bam.candidates.push_back(snp(350, VariantCategory::NoisyCandHet, kCandNoisyCandHet));  // 3: not verified
    for (hts_pos_t pos : {403, 406}) {                                           // 4, 5: one allele
        CandidateVariant ins;
        ins.key.tid = 0;
        ins.key.pos = pos;
        ins.key.type = VariantType::Insertion;
        ins.key.ref_len = 0;
        ins.key.alt = "A";
        ins.counts.category = VariantCategory::NoisyCandHet;
        ins.lcd_var_i_to_cate = kCandNoisyCandHet;
        ins.msa_verified = true;
        bam.candidates.push_back(ins);
    }
    for (int r = 0; r < 6; ++r) {
        const int a = r % 2;
        bam.reads.push_back(read_named("r" + std::to_string(r), 60));
        bam.read_var_profile.push_back(profile(r, 0, {a, a, 0, a, a, a}));
    }
    bam.reads.push_back(read_named("x6", 5));
    bam.read_var_profile.push_back(profile(6, 2, {1}));
    bam.reads.push_back(read_named("x7", 5));
    bam.read_var_profile.push_back(profile(7, 2, {1}));
    bam.reads.push_back(read_named("x8", 60));
    bam.read_var_profile.push_back(profile(8, 1, {1}));

    const size_t added = inject_alignment_sites(gc, bam, "chr1");
    // SNP 100 overlaps the voting graph row: it joins as a hidden complement
    // site (no VCF record) right after that row.
    check(added == 3, "inject: three sites added, got " + std::to_string(added));
    check(graph.candidates.size() == 4, "inject: four rows after merge");
    if (graph.candidates.size() != 4) return;
    check(!graph.candidates[0].bam_injected && graph.candidates[0].key.pos == 100, "inject: graph row first");
    check(graph.candidates[1].bam_injected && graph.candidates[1].key.pos == 100 && gc.site_meta[1].ref.empty(),
          "inject: overlapping call is a hidden complement site");
    check(graph.candidates[2].bam_injected && graph.candidates[2].key.pos == 200, "inject: SNP 200 injected");
    check(graph.candidates[3].bam_injected && graph.candidates[3].key.type == VariantType::Insertion,
          "inject: one homopolymer insertion injected");
    check(gc.site_meta.size() == 4 && gc.site_meta[2].pos == 200 && gc.site_meta[2].ref == ref.substr(199, 1),
          "inject: VCF metadata synthesised");
    check(graph.reads.size() == 7 && graph.reads[6].qname == "x8",
          "inject: only the alignment-only read with a kept call is added");
    const auto allele_at = [&](size_t read, size_t site) {
        const ReadVariantProfile& p = graph.read_var_profile[read];
        if (p.start_var_idx < 0 || static_cast<int>(site) < p.start_var_idx ||
            static_cast<int>(site) > p.end_var_idx) return -2;
        return p.alleles[site - static_cast<size_t>(p.start_var_idx)];
    };
    check(allele_at(0, 0) == 0 && allele_at(1, 0) == 1, "inject: graph observations kept");
    check(allele_at(0, 2) == 0 && allele_at(1, 2) == 1, "inject: alignment observations transferred");
    check(allele_at(5, 3) == 1, "inject: alias placement folds into the canonical row");
    check(allele_at(6, 2) == 1, "inject: alignment-only read carries its call");
    check(graph.haps.size() == 7 && graph.phase_sets.size() == 7, "inject: read arrays resized");
}

// ── Local windows ───────────────────────────────────────────────────────────

static std::unique_ptr<bam1_t, AlignmentDeleter> make_alignment(const std::string& name, const std::string& seq,
                                                                std::vector<uint32_t> cigar) {
    std::unique_ptr<bam1_t, AlignmentDeleter> b(bam_init1());
    const std::string qual(seq.size(), static_cast<char>(30));
    bam_set1(b.get(), name.size(), name.c_str(), 0, 0, 0, 60, cigar.size(), cigar.data(), -1, -1, 0,
             seq.size(), seq.c_str(), qual.c_str(), 0);
    return b;
}

static void window_case(bool heterozygous, bool seeded = true) {
    std::string ref;
    uint32_t state = 12345;
    for (int i = 0; i < 600; ++i) {
        state = state * 1103515245u + 12345u;
        ref += "ACGT"[(state >> 16) % 4];
    }
    const std::string ins = "GATTACAG";

    PhasingChunk bam;
    bam.ref_beg = 1;
    bam.ref_end = 600;
    bam.ref_seq = ref;
    CandidateVariant call;
    call.key.tid = 0;
    call.key.pos = 300;
    call.key.type = VariantType::Insertion;
    call.key.alt = ins;
    call.counts.category = VariantCategory::NoisyCandHet;
    bam.candidates.push_back(call);

    // Graph chunk: three phased sites label the reads (hap 1 carries ALT).
    PhasingChunk graph;
    for (int k = 0; k < 3; ++k) {
        CandidateVariant c = clean_snp(50 + 10 * k);
        c.phase_set = seeded ? 50 : 0;  // unseeded: no phased site labels the reads
        c.hap_to_cons_alle = {-1, 1, 0};
        graph.candidates.push_back(c);
    }
    for (int r = 0; r < 8; ++r) {
        const int hap = r % 2 == 0 ? 1 : 2;
        const std::string name = "w" + std::to_string(r);
        graph.reads.push_back(read_named(name, 60));
        const int a = hap == 1 ? 1 : 0;
        graph.read_var_profile.push_back(profile(r, 0, {a, a, a}));
        ReadRecord rec = read_named(name, 60);
        const bool carries = heterozygous ? hap == 1 : true;
        if (carries) {
            rec.alignment = make_alignment(name, ref.substr(0, 299) + ins + ref.substr(299),
                {bam_cigar_gen(299, BAM_CMATCH), bam_cigar_gen(ins.size(), BAM_CINS), bam_cigar_gen(301, BAM_CMATCH)});
        } else {
            rec.alignment = make_alignment(name, ref, {bam_cigar_gen(600, BAM_CMATCH)});
        }
        bam.reads.push_back(std::move(rec));
    }
    finish_chunk(graph);

    const std::vector<LocusWindowSite> loci = build_locus_window_sites(bam, graph);
    if (!heterozygous) {
        check(loci.empty(), "windows: identical haplotypes give no window");
        return;
    }
    check(loci.size() == 1, "windows: one heterozygous window, got " + std::to_string(loci.size()));
    if (loci.size() != 1) return;
    check(loci[0].observations.size() == 8, "windows: every spanning read assigned");
    // Seeded windows put haplotype 1 on side 0; an unseeded split may use
    // either orientation but must still separate the haplotypes.
    const int even_side = std::get<1>(loci[0].observations.front()) ^ (std::get<0>(loci[0].observations.front()) % 2);
    bool sides_ok = !seeded || even_side == 0;
    for (const auto& [read, side, error] : loci[0].observations)
        sides_ok = sides_ok && side == (static_cast<int>(read % 2) ^ even_side) && error < 0.01f;
    check(sides_ok, std::string("windows: reads assigned to their own haplotype with small error") +
                        (seeded ? "" : " (unseeded split)"));
}


static void test_realign_recalls_indel() {
    // An injected 8-base insertion that half the reads have no call at (the MSA
    // did not cover them) and half carry an MSA call: the re-call fills the
    // missing ones from each read's sequence and keeps the MSA's own calls.
    std::string ref;
    uint32_t state = 777;
    for (int i = 0; i < 600; ++i) {
        state = state * 1103515245u + 12345u;
        ref += "ACGT"[(state >> 16) % 4];
    }
    const std::string ins = "GATTACAG";
    PhasingChunk bam;
    bam.ref_beg = 1;
    bam.ref_end = 600;
    bam.ref_seq = ref;
    PhasingChunk chunk;
    CandidateVariant c;
    c.key.tid = 0;
    c.key.pos = 300;
    c.key.type = VariantType::Insertion;
    c.key.alt = ins;
    c.counts.category = VariantCategory::NoisyCandHet;
    c.lcd_var_i_to_cate = kCandNoisyCandHet;
    c.bam_injected = true;
    chunk.candidates.push_back(c);
    for (int r = 0; r < 8; ++r) {
        const std::string name = "q" + std::to_string(r);
        const bool carries = r % 2 == 0;
        chunk.reads.push_back(read_named(name, 60));
        // Reads 0-3 lack a call; reads 4-7 carry the MSA's (correct) call.
        chunk.read_var_profile.push_back(profile(r, 0, {r < 4 ? -1 : (carries ? 1 : 0)}));
        ReadRecord rec = read_named(name, 60);
        if (carries)
            rec.alignment = make_alignment(name, ref.substr(0, 299) + ins + ref.substr(299),
                {bam_cigar_gen(299, BAM_CMATCH), bam_cigar_gen(ins.size(), BAM_CINS), bam_cigar_gen(301, BAM_CMATCH)});
        else
            rec.alignment = make_alignment(name, ref, {bam_cigar_gen(600, BAM_CMATCH)});
        bam.reads.push_back(std::move(rec));
    }
    finish_chunk(chunk);
    ObservationErrors errors;
    const size_t changed = realign_indel_observations(chunk, bam, &errors);
    check(changed == 4, "realign: the four missing calls filled, got " + std::to_string(changed));
    bool ok = true;
    for (int r = 0; r < 8; ++r) ok = ok && chunk.read_var_profile[static_cast<size_t>(r)].alleles[0] == (r % 2 == 0 ? 1 : 0);
    check(ok, "realign: each read takes the allele its sequence carries");
    // Every graded call carries an error below one half; exact segments are confident.
    bool graded = !errors.empty();
    for (const auto& [key, e] : errors) graded = graded && e > 0.0 && e < 0.5;
    check(graded, "realign: calls carry a per-read error in (0, 0.5)");
}

static void test_realign_recalls_merged_deletion() {
    // Two deletions of different length at one anchor merged into one row
    // (merge_msa_colocated_deletions): REF spans the longer one, allele 1
    // deletes all of it, allele 2 keeps the bases the shorter one leaves.
    // Reads the MSA did not call are re-called against all three.
    std::string ref;
    uint32_t state = 4242;
    for (int i = 0; i < 600; ++i) {
        state = state * 1103515245u + 12345u;
        ref += "ACGT"[(state >> 16) % 4];
    }
    PhasingChunk bam;
    bam.ref_beg = 1;
    bam.ref_end = 600;
    bam.ref_seq = ref;
    PhasingChunk chunk;
    CandidateVariant c;
    c.key.tid = 0;
    c.key.pos = 300;  // first deleted base, 1-based
    c.key.type = VariantType::Deletion;
    c.key.ref_len = 6;
    c.msa_insertion_alts = {std::string(), ref.substr(301, 4)};
    c.counts.category = VariantCategory::NoisyCandHet;
    c.lcd_var_i_to_cate = kCandNoisyCandHet;
    c.bam_injected = true;
    chunk.candidates.push_back(c);
    for (int r = 0; r < 9; ++r) {
        const std::string name = "d" + std::to_string(r);
        chunk.reads.push_back(read_named(name, 60));
        chunk.read_var_profile.push_back(profile(r, 0, {-1}));
        ReadRecord rec = read_named(name, 60);
        const int del = r % 3 == 0 ? 0 : r % 3 == 1 ? 6 : 2;  // allele 0, 1, 2
        if (del == 0)
            rec.alignment = make_alignment(name, ref, {bam_cigar_gen(600, BAM_CMATCH)});
        else
            rec.alignment = make_alignment(name, ref.substr(0, 299) + ref.substr(299 + del),
                {bam_cigar_gen(299, BAM_CMATCH), bam_cigar_gen(del, BAM_CDEL), bam_cigar_gen(301 - del, BAM_CMATCH)});
        bam.reads.push_back(std::move(rec));
    }
    finish_chunk(chunk);
    const size_t changed = realign_indel_observations(chunk, bam, nullptr);
    check(changed == 9, "realign: every read of a merged deletion row called, got " + std::to_string(changed));
    bool ok = true;
    for (int r = 0; r < 9; ++r) ok = ok && chunk.read_var_profile[static_cast<size_t>(r)].alleles[0] == r % 3;
    check(ok, "realign: each read takes its own deletion allele");
}

static void test_bridge_weak_sites() {
    // A one-base insertion into a long homopolymer is a weak bridge; a SNP and
    // a one-base insertion outside any run are not.
    std::string ref;
    uint32_t state = 99;
    for (int i = 0; i < 300; ++i) {
        state = state * 1103515245u + 12345u;
        ref += "ACGT"[(state >> 16) % 4];
    }
    ref.replace(99, 10, "AAAAAAAAAA");  // 1-based 100-109: A10
    ref[59] = 'C'; ref[60] = 'G'; ref[58] = 'T';  // 1-based 59-61: TCG (no run)
    PhasingChunk bam;
    bam.ref_beg = 1;
    bam.ref_end = 300;
    bam.ref_seq = ref;
    GraphChunkBuildResult gc;
    PhasingChunk& chunk = gc.chunk;
    const auto injected_insertion = [](hts_pos_t pos, const std::string& alt) {
        CandidateVariant c;
        c.key.tid = 0;
        c.key.pos = pos;
        c.key.type = VariantType::Insertion;
        c.key.alt = alt;
        c.lcd_var_i_to_cate = kCandNoisyCandHet;
        c.bam_injected = true;
        return c;
    };
    chunk.candidates.push_back(injected_insertion(105, "A"));  // into A10
    chunk.candidates.push_back(injected_insertion(61, "A"));   // between C and G
    chunk.candidates.push_back(clean_snp(200));
    finish_chunk(chunk);
    const std::vector<char> weak = bridge_weak_sites(gc, bam);
    check(weak.size() == 3 && weak[0] == 1, "weak: homopolymer length site is weak");
    check(weak.size() == 3 && weak[1] == 0, "weak: insertion outside a run is not weak");
    check(weak.size() == 3 && weak[2] == 0, "weak: SNP is not weak");
}

static void test_em_join_needs_strong_reads() {
    // Sites 0-3 and 4-7, each side tiled; only `crossers` reads span sites 3
    // and 4, one call per side. Two such reads clear the ordinary cut but are
    // too few to make the join strong: it needs the weak-bridge margin.
    const auto blocks_with = [](int crossers) {
        PhasingChunk chunk;
        std::vector<int> truth;
        for (int k = 0; k < 8; ++k) chunk.candidates.push_back(clean_snp(100 + 100 * k));
        for (int rep = 0; rep < 2; ++rep) {
            add_tiled_reads(chunk, truth, 0, 3, 4);
            add_tiled_reads(chunk, truth, 4, 7, 4);
        }
        for (int c = 0; c < crossers; ++c) {
            const int hap = c % 2;
            const int id = static_cast<int>(chunk.reads.size());
            chunk.reads.push_back(read_named("x" + std::to_string(c), 60));
            chunk.read_var_profile.push_back(profile(id, 3, {hap_allele(hap, 3), hap_allele(hap, 4)}));
            truth.push_back(hap);
        }
        finish_chunk(chunk);
        const std::vector<char> no_weak(chunk.candidates.size(), 0);
        return phase_chunk_by_global_em(chunk, nullptr, &no_weak);
    };
    check(blocks_with(2) == 2, "EM: two crossing reads do not make a strong join");
    check(blocks_with(4) == 1, "EM: four crossing reads join the blocks");
}

static void test_msa_absence_of_alt() {
    // The haplotype without this ALT may carry another event at the site; a
    // read matching its consensus exactly carries absence of the ALT (0).
    VariantKey del;
    del.tid = 0;
    del.pos = 100;
    del.type = VariantType::Deletion;
    del.ref_len = 12;
    const std::string repeat = "AAAGAAAGAAAG";
    const std::array<std::string, 2> two_deletions = {"AAAAGAAAG", ""};  // 3 bp and 12 bp
    check(msa_site_allele_from_sequences(del, repeat, two_deletions, "") == 1, "msa: 12 bp deletion read is ALT");
    check(msa_site_allele_from_sequences(del, repeat, two_deletions, "AAAAGAAAG") == 0,
          "msa: 3 bp deletion read is absence of the 12 bp ALT");
    check(msa_site_allele_from_sequences(del, repeat, two_deletions, "AAAGAAAG") == -1,
          "msa: a third sequence stays uncalled");
    VariantKey snv;
    snv.tid = 0;
    snv.pos = 200;
    snv.type = VariantType::Snp;
    snv.ref_len = 1;
    snv.alt = "A";
    const std::array<std::string, 2> insertion_before = {"AC", "A"};  // other haplotype: +A before the base
    check(msa_site_allele_from_sequences(snv, "C", insertion_before, "A") == 1, "msa: SNP read is ALT");
    check(msa_site_allele_from_sequences(snv, "C", insertion_before, "AC") == 0,
          "msa: read with the other haplotype's insertion is absence of the SNP");
    check(msa_site_allele_from_sequences(snv, "C", insertion_before, "GC") == -1, "msa: unrelated read uncalled");
}

static void test_confirm_catalog_snp_calls() {
    // Catalog SNP at 100 (REF A, ALT T). Read 0's graph call disagrees with its
    // alignment call; read 1 agrees; read 2 has no alignment call there.
    GraphChunkBuildResult gc;
    PhasingChunk& graph = gc.chunk;
    graph.candidates.push_back(clean_snp(100));
    for (int r = 0; r < 3; ++r) {
        graph.reads.push_back(read_named("r" + std::to_string(r), 60));
        graph.read_var_profile.push_back(profile(r, 0, {1}));
    }
    GraphSiteMeta meta;
    meta.pos = 100;
    meta.ref = "A";
    meta.alts = {"T"};
    gc.site_meta = {meta};
    gc.site_allele_orig_idx = {{0, 1}};
    finish_chunk(graph);
    PhasingChunk bam;
    bam.candidates.push_back(clean_snp(100));
    for (int r = 0; r < 3; ++r) bam.reads.push_back(read_named("r" + std::to_string(r), 60));
    bam.read_var_profile.push_back(profile(0, 0, {0}));   // alignment: REF
    bam.read_var_profile.push_back(profile(1, 0, {1}));   // alignment: ALT
    ReadVariantProfile none;
    none.read_id = 2;
    bam.read_var_profile.push_back(none);                 // no alignment call
    ObservationErrors errors;
    const size_t changed = confirm_catalog_snp_calls(gc, bam, &errors);
    check(changed == 1 && graph.read_var_profile[0].alleles[0] == 0, "confirm: disagreeing graph call takes the alignment's");
    check(graph.read_var_profile[1].alleles[0] == 1 && errors.count(observation_key(1, 0)) == 0,
          "confirm: agreeing call untouched");
    check(graph.read_var_profile[2].alleles[0] == 1 && errors.count(observation_key(2, 0)) == 1 &&
              errors[observation_key(2, 0)] >= 0.2,
          "confirm: graph-only call kept with a per-read error");
}

static void test_inject_island_only_sites() {
    // An MSA-verified het whose ALT only low-MAPQ reads carry fails the
    // confident-read rule: it enters as an island-only site (no germline
    // bits, no VCF record). One at an indel a regular site already carries
    // (same left-aligned allele) is dropped, never made canonical.
    std::string ref;
    for (int i = 0; i < 500; ++i) ref += "ACGT"[(i * 7 + i / 3) % 4];
    ref.replace(399, 7, "CAAAAAG");  // 1-based 400: C, 401-405: A run, 406: G
    for (size_t i : {99, 199, 249}) ref[i] = 'A';

    GraphChunkBuildResult gc;
    PhasingChunk& graph = gc.chunk;
    graph.ref_beg = 1;
    graph.ref_end = 500;
    graph.candidates.push_back(clean_snp(100));
    for (int r = 0; r < 6; ++r) {
        graph.reads.push_back(read_named("r" + std::to_string(r), 60));
        graph.read_var_profile.push_back(profile(r, 0, {r % 2}));
    }
    gc.site_ids = {"s100"};
    gc.site_meta.resize(1);
    gc.site_allele_orig_idx = {{0, 1}};
    finish_chunk(graph);

    const auto msa_ins = [](hts_pos_t pos) {
        CandidateVariant ins;
        ins.key.tid = 0;
        ins.key.pos = pos;
        ins.key.type = VariantType::Insertion;
        ins.key.ref_len = 0;
        ins.key.alt = "A";
        ins.counts.category = VariantCategory::NoisyCandHet;
        ins.lcd_var_i_to_cate = kCandNoisyCandHet;
        ins.msa_verified = true;
        return ins;
    };
    PhasingChunk bam;
    bam.ref_beg = 1;
    bam.ref_end = 500;
    bam.ref_seq = ref;
    bam.candidates.push_back(clean_snp(200));                                    // 0: injected
    CandidateVariant weak = snp(250, VariantCategory::NoisyCandHet, kCandNoisyCandHet);
    weak.msa_verified = true;
    bam.candidates.push_back(weak);                                              // 1: island-only
    bam.candidates.push_back(msa_ins(403));                                      // 2: island-only, aliases 3
    bam.candidates.push_back(msa_ins(406));                                      // 3: injected
    for (int r = 0; r < 6; ++r) {
        const int a = r % 2;
        bam.reads.push_back(read_named("r" + std::to_string(r), 60));
        bam.read_var_profile.push_back(profile(r, 0, {a, 0, 0, a}));
    }
    for (int r = 6; r < 8; ++r) {
        bam.reads.push_back(read_named("x" + std::to_string(r), 5));
        bam.read_var_profile.push_back(profile(r, 1, {1, 1}));
    }
    inject_alignment_sites(gc, bam, "chr1");
    int at250 = -1, n_ins = 0;
    for (size_t ci = 0; ci < graph.candidates.size(); ++ci) {
        const CandidateVariant& c = graph.candidates[ci];
        if (c.key.pos == 250) at250 = static_cast<int>(ci);
        if (c.key.type == VariantType::Insertion) {
            ++n_ins;
            check(!c.island_only && (c.lcd_var_i_to_cate & kCandGermlineVarCate) != 0,
                  "inject: the regular insertion stays canonical");
        }
    }
    check(n_ins == 1, "inject: one insertion row, got " + std::to_string(n_ins));
    check(at250 >= 0, "inject: the low-confidence MSA het enters");
    if (at250 < 0) return;
    const CandidateVariant& c = graph.candidates[static_cast<size_t>(at250)];
    check(c.island_only && c.bam_injected, "inject: it is island-only");
    check((c.lcd_var_i_to_cate & kCandGermlineVarCate) == 0, "inject: it carries no germline bits");
    check(gc.site_meta[static_cast<size_t>(at250)].ref.empty(), "inject: it writes no VCF record");
}

static void test_inject_low_mapq_clean_het() {
    // A clean het read only at low MAPQ is injected as an ordinary site (the
    // confident-read rule binds MSA-verified hets only); an MSA het read the
    // same way stays island-only.
    std::string ref;
    for (int i = 0; i < 500; ++i) ref += "ACGT"[(i * 7 + i / 3) % 4];
    for (size_t i : {99, 199, 299, 349}) ref[i] = 'A';
    GraphChunkBuildResult gc;
    PhasingChunk& graph = gc.chunk;
    graph.ref_beg = 1;
    graph.ref_end = 500;
    graph.candidates.push_back(clean_snp(100));
    for (int r = 0; r < 6; ++r) {
        graph.reads.push_back(read_named("r" + std::to_string(r), 60));
        graph.read_var_profile.push_back(profile(r, 0, {r % 2}));
    }
    gc.site_ids = {"s100"};
    gc.site_meta.resize(1);
    gc.site_allele_orig_idx = {{0, 1}};
    finish_chunk(graph);
    PhasingChunk bam;
    bam.ref_beg = 1;
    bam.ref_end = 500;
    bam.ref_seq = ref;
    bam.candidates.push_back(clean_snp(200));                                    // 0: confident reads
    bam.candidates.push_back(clean_snp(300));                                    // 1: low MAPQ only
    CandidateVariant weak = snp(350, VariantCategory::NoisyCandHet, kCandNoisyCandHet);
    weak.msa_verified = true;
    bam.candidates.push_back(weak);                                              // 2: low MAPQ only
    for (int r = 0; r < 6; ++r) {
        bam.reads.push_back(read_named("r" + std::to_string(r), 60));
        bam.read_var_profile.push_back(profile(r, 0, {r % 2}));
    }
    for (int r = 6; r < 12; ++r) {
        bam.reads.push_back(read_named("x" + std::to_string(r), 5));
        bam.read_var_profile.push_back(profile(r, 1, {r % 2, r % 2}));
    }
    inject_alignment_sites(gc, bam, "chr1");
    int at300 = -1, at350 = -1;
    for (size_t ci = 0; ci < graph.candidates.size(); ++ci) {
        if (graph.candidates[ci].key.pos == 300) at300 = static_cast<int>(ci);
        if (graph.candidates[ci].key.pos == 350) at350 = static_cast<int>(ci);
    }
    check(at300 >= 0, "low-MAPQ clean het: injected");
    if (at300 >= 0) {
        const CandidateVariant& c = graph.candidates[static_cast<size_t>(at300)];
        check(!c.island_only && (c.lcd_var_i_to_cate & kCandGermlineVarCate) != 0 && c.low_mapq_only,
              "low-MAPQ clean het: a germline site flagged low_mapq_only");
    }
    for (const CandidateVariant& c : graph.candidates)
        if (c.key.pos == 200) check(!c.low_mapq_only, "confidently read clean het: not flagged");
    check(at350 >= 0 && graph.candidates[static_cast<size_t>(at350)].island_only,
          "low-MAPQ MSA het: island-only");
}

static void test_em_island_blocks() {
    // Sites 0-5 clean SNPs, tiled; site 6 a repeat het indel and site 7 an
    // island-only MSA het, neither an EM site, read by labelled reads that
    // agree with their haplotypes. Two reads each call only site 6 or only
    // site 7: the solve cannot label them, the islands do, by their own calls.
    PhasingChunk chunk;
    for (int k = 0; k < 6; ++k) chunk.candidates.push_back(clean_snp(100 + 100 * k));
    CandidateVariant rhi = snp(700, VariantCategory::RepeatHetIndel, 0);
    chunk.candidates.push_back(rhi);
    CandidateVariant only = snp(800, VariantCategory::NoisyCandHet, 0);
    only.island_only = true;
    chunk.candidates.push_back(only);
    std::vector<int> truth;
    add_tiled_reads(chunk, truth, 0, 5, 3);
    add_tiled_reads(chunk, truth, 4, 6, 3);
    add_tiled_reads(chunk, truth, 5, 7, 3);
    const size_t first_lone = chunk.reads.size();
    for (int site : {6, 7}) {
        for (int hap = 0; hap < 2; ++hap) {
            const int id = static_cast<int>(chunk.reads.size());
            chunk.reads.push_back(read_named("lone" + std::to_string(id), 60));
            chunk.read_var_profile.push_back(profile(id, site, {hap_allele(hap, static_cast<size_t>(site))}));
            truth.push_back(hap);
        }
    }
    finish_chunk(chunk);
    const std::vector<LocusWindowSite> no_loci;
    phase_chunk_by_global_em(chunk, &no_loci, nullptr, nullptr);
    for (size_t r = first_lone; r < chunk.reads.size(); ++r)
        check(chunk.haps[r] == 0, "islands: the solve leaves lone read " + std::to_string(r) + " unlabelled");
    check(chunk.island_haps.size() == chunk.reads.size(), "islands: island labels recorded");
    if (chunk.island_haps.size() != chunk.reads.size()) return;
    for (size_t pair = 0; pair < 2; ++pair) {
        const size_t a = first_lone + 2 * pair, b = a + 1;
        check(chunk.island_haps[a] != 0 && chunk.island_haps[b] != 0 && chunk.island_haps[a] != chunk.island_haps[b],
              "islands: lone reads at site " + std::to_string(6 + pair) + " take opposite haplotypes");
        check(chunk.island_phase_sets[a] == chunk.island_phase_sets[b] &&
              chunk.island_phase_sets[a] == chunk.candidates[6 + pair].key.sort_pos(),
              "islands: one phase set per island, at its site");
    }
    // Applied after the rescue: a read the rescue labelled keeps its label.
    chunk.gap_haps.assign(chunk.reads.size(), 0);
    chunk.gap_phase_sets.assign(chunk.reads.size(), 0);
    chunk.gap_haps[first_lone + 3] = 1;
    chunk.gap_phase_sets[first_lone + 3] = 100 + kGapFillPsOffset;
    const size_t applied = apply_island_labels(chunk);
    check(applied == 3, "islands: three lone reads labelled, got " + std::to_string(applied));
    check(chunk.gap_haps[first_lone] == chunk.island_haps[first_lone] &&
          chunk.gap_phase_sets[first_lone] == chunk.island_phase_sets[first_lone] + kGapFillPsOffset &&
          chunk.gap_from_bam_observation[first_lone],
          "islands: label goes through the gap channel, fill only");
    check(chunk.gap_haps[first_lone + 3] == 1 && chunk.gap_phase_sets[first_lone + 3] == 100 + kGapFillPsOffset,
          "islands: an earlier gap label is kept");
    for (size_t r = 0; r < first_lone; ++r)
        check(chunk.gap_haps[r] == 0, "islands: solve-labelled read " + std::to_string(r) + " untouched");
    check(chunk.candidates[6].phase_set == chunk.island_phase_sets[first_lone] &&
          chunk.candidates[6].hap_to_cons_alle[1] != chunk.candidates[6].hap_to_cons_alle[2],
          "islands: the island site is phased in its own block");
}

static void test_em_deferred_labels() {
    // A read whose one call carries a high per-read error stays below the
    // label posterior: its calls agree, so it takes its block's label last,
    // through the island channel. A read whose calls disagree does not.
    PhasingChunk chunk;
    for (int k = 0; k < 6; ++k) chunk.candidates.push_back(clean_snp(100 + 100 * k));
    std::vector<int> truth;
    add_tiled_reads(chunk, truth, 0, 5, 3);
    const size_t weak = chunk.reads.size();
    chunk.reads.push_back(read_named("weak", 60));
    chunk.read_var_profile.push_back(profile(static_cast<int>(weak), 2, {hap_allele(0, 2)}));
    truth.push_back(0);
    const size_t torn = chunk.reads.size();
    chunk.reads.push_back(read_named("torn", 60));
    chunk.read_var_profile.push_back(profile(static_cast<int>(torn), 2, {hap_allele(0, 2), hap_allele(1, 3)}));
    truth.push_back(0);
    finish_chunk(chunk);
    ObservationErrors errors;
    errors[observation_key(weak, 2)] = 0.25;
    errors[observation_key(torn, 2)] = 0.25;
    errors[observation_key(torn, 3)] = 0.25;
    const std::vector<LocusWindowSite> no_loci;
    phase_chunk_by_global_em(chunk, &no_loci, nullptr, &errors);
    check(chunk.haps[weak] == 0 && chunk.haps[torn] == 0, "deferred: both reads stay below the label posterior");
    check(chunk.island_haps.size() == chunk.reads.size(), "deferred: labels recorded");
    if (chunk.island_haps.size() != chunk.reads.size()) return;
    // Its label agrees with a solve-labelled read on the same haplotype.
    size_t mate = 0;
    while (mate < weak && !(truth[mate] == 0 && chunk.haps[mate] != 0)) ++mate;
    check(mate < weak && chunk.island_haps[weak] == chunk.haps[mate] &&
          chunk.island_phase_sets[weak] == chunk.phase_sets[mate],
          "deferred: the agreeing read takes its block's label");
    check(chunk.island_haps[torn] == 0, "deferred: a read whose calls disagree takes none");
    check(apply_island_labels(chunk) == 1 && chunk.gap_haps[weak] == chunk.island_haps[weak],
          "deferred: applied through the gap channel");
}

static void test_em_echo_sites() {
    // Sites 2 (300) and 3 (310) split every read the same way: one event
    // written twice. A read with calls at sites 1-4 carries the event on the
    // other haplotype; counted twice it ties the read's two other calls and
    // leaves the read unlabelled, counted once it is outvoted. Site 3 is
    // still phased.
    PhasingChunk chunk;
    for (hts_pos_t pos : {100, 200, 300, 310, 400, 500, 600}) chunk.candidates.push_back(clean_snp(pos));
    std::vector<int> truth;
    add_tiled_reads(chunk, truth, 0, 6, 5);
    add_tiled_reads(chunk, truth, 0, 6, 5);
    const size_t odd = chunk.reads.size();
    chunk.reads.push_back(read_named("odd", 60));
    chunk.read_var_profile.push_back(
        profile(static_cast<int>(odd), 1, {hap_allele(0, 1), hap_allele(1, 2), hap_allele(1, 3), hap_allele(0, 4)}));
    truth.push_back(0);
    finish_chunk(chunk);
    const std::vector<LocusWindowSite> no_loci;
    phase_chunk_by_global_em(chunk, &no_loci, nullptr, nullptr);
    size_t mate = 0;
    while (mate < odd && !(truth[mate] == 0 && chunk.haps[mate] != 0 && chunk.phase_sets[mate] == chunk.phase_sets[0]))
        ++mate;
    check(mate < odd, "echo: a haplotype-0 read is labelled");
    check(chunk.haps[odd] != 0 && chunk.haps[odd] == chunk.haps[mate] && chunk.phase_sets[odd] == chunk.phase_sets[mate],
          "echo: the read takes its other calls' haplotype");
    check(chunk.candidates[3].phase_set == chunk.candidates[2].phase_set &&
          chunk.candidates[3].hap_to_cons_alle[1] != chunk.candidates[3].hap_to_cons_alle[2],
          "echo: the echo site is phased with its block");
}

static void test_em_unverified_fallback(bool low_mapq) {
    // Sites 4 and 5 are unverified calls (a noisy region too long for MSA) or,
    // with low_mapq, clean hets no confidently mapped read carries -- treated
    // the same. They sit between clean sites, learned from tiled reads that
    // cross them. A read whose only calls are there gets no solve or deferred
    // label, then takes its block's label from them, last and fill-only.
    PhasingChunk chunk;
    for (int k = 0; k < 10; ++k) {
        chunk.candidates.push_back(clean_snp(100 + 100 * k));
        if (k == 4 || k == 5) (low_mapq ? chunk.candidates.back().low_mapq_only : chunk.candidates.back().unresolved_noisy) = true;
    }
    std::vector<int> truth;
    add_tiled_reads(chunk, truth, 0, 9, 5);
    add_tiled_reads(chunk, truth, 0, 9, 5);
    const size_t lone = chunk.reads.size();
    chunk.reads.push_back(read_named("lone", 60));
    chunk.read_var_profile.push_back(profile(static_cast<int>(lone), 4, {hap_allele(0, 4), hap_allele(0, 5)}));
    truth.push_back(0);
    finish_chunk(chunk);
    const std::vector<LocusWindowSite> no_loci;
    phase_chunk_by_global_em(chunk, &no_loci, nullptr, nullptr);
    check(chunk.haps[lone] == 0, "unverified: the solve leaves the read unlabelled");
    size_t mate = 0;
    while (mate < lone && !(truth[mate] == 0 && chunk.haps[mate] != 0 && chunk.phase_sets[mate] == chunk.candidates[4].phase_set))
        ++mate;
    check(mate < lone, "unverified: a haplotype-0 read is labelled in the unverified sites' block");
    check(chunk.island_haps.size() == chunk.reads.size() && chunk.island_haps[lone] != 0 &&
          chunk.island_haps[lone] == chunk.haps[mate] && chunk.island_phase_sets[lone] == chunk.phase_sets[mate],
          "unverified: the read takes its block's label from its unverified calls");
    check(apply_island_labels(chunk) == 1 && chunk.gap_haps[lone] == chunk.island_haps[lone],
          "unverified: applied through the gap channel");
}

static void test_length_window_sites() {
    // Five 10-A runs, each with an indel candidate in the alignment solve. At
    // the run at 300 half the reads carry +1 A and half +3 A; every other run
    // is read exactly, which teaches the stutter model. One window, over the
    // het run, splitting the reads by length.
    std::string ref;
    uint32_t state = 99;
    for (int i = 0; i < 700; ++i) {
        state = state * 1103515245u + 12345u;
        char b = "CGT"[(state >> 16) % 3];  // no A outside the runs
        ref += b;
    }
    for (int run : {100, 200, 300, 400, 500}) ref.replace(static_cast<size_t>(run), 10, std::string(10, 'A'));
    PhasingChunk bam;
    bam.ref_beg = 1;
    bam.ref_end = 700;
    bam.ref_seq = ref;
    for (int run : {100, 200, 300, 400, 500}) {
        CandidateVariant c;
        c.key.tid = 0;
        c.key.pos = run + 1;  // insertion after the run's first base (1-based anchor + 1)
        c.key.type = VariantType::Insertion;
        c.key.ref_len = 0;
        c.key.alt = "A";
        c.counts.category = VariantCategory::RepeatHetIndel;
        bam.candidates.push_back(c);
    }
    GraphChunkBuildResult gc;
    PhasingChunk& graph = gc.chunk;
    graph.ref_beg = 1;
    graph.ref_end = 700;
    for (int r = 0; r < 24; ++r) {
        const std::string name = "len" + std::to_string(r);
        const int extra = r % 2 == 0 ? 1 : 3;
        const std::string seq = ref.substr(0, 305) + std::string(static_cast<size_t>(extra), 'A') + ref.substr(305);
        ReadRecord rec = read_named(name, 60);
        rec.alignment = make_alignment(name, seq,
            {bam_cigar_gen(305, BAM_CMATCH), bam_cigar_gen(static_cast<uint32_t>(extra), BAM_CINS),
             bam_cigar_gen(395, BAM_CMATCH)});
        bam.reads.push_back(std::move(rec));
        graph.reads.push_back(read_named(name, 60));
        graph.read_var_profile.push_back(ReadVariantProfile{});
    }
    finish_chunk(graph);
    const std::vector<LocusWindowSite> windows = build_length_window_sites(bam, gc, {});
    check(windows.size() == 1, "length windows: one window, at the het run, got " + std::to_string(windows.size()));
    if (windows.size() != 1) return;
    const LocusWindowSite& w = windows[0];
    check(w.pos <= 301 && w.end >= 310, "length windows: the window covers the run");
    std::set<int> side_of[2];
    for (const auto& [ri, side, e] : w.observations) side_of[ri % 2].insert(side);
    check(w.observations.size() == 24, "length windows: every read called, got " + std::to_string(w.observations.size()));
    check(side_of[0].size() == 1 && side_of[1].size() == 1 && *side_of[0].begin() != *side_of[1].begin(),
          "length windows: reads split by their length");
}

static void test_em_length_window_islands() {
    // An island-only window is no EM site, but phases the reads the solve
    // leaves unlabelled, in a block of its own at the window.
    PhasingChunk chunk;
    for (int k = 0; k < 6; ++k) chunk.candidates.push_back(clean_snp(100 + 100 * k));
    std::vector<int> truth;
    add_tiled_reads(chunk, truth, 0, 5, 3);
    const size_t first_lone = chunk.reads.size();
    for (int r = 0; r < 4; ++r) {
        chunk.reads.push_back(read_named("lone" + std::to_string(r), 60));
        chunk.read_var_profile.push_back(ReadVariantProfile{});
        truth.push_back(r % 2);
    }
    finish_chunk(chunk);
    LocusWindowSite w;
    w.pos = 900;
    w.end = 940;
    w.island_only = true;
    for (int r = 0; r < 4; ++r) w.observations.emplace_back(first_lone + static_cast<size_t>(r), r % 2, 0.05f);
    const std::vector<LocusWindowSite> loci{w};
    phase_chunk_by_global_em(chunk, &loci, nullptr, nullptr);
    check(chunk.island_haps.size() == chunk.reads.size(), "length islands: island labels recorded");
    if (chunk.island_haps.size() != chunk.reads.size()) return;
    bool ok = true;
    for (int r = 0; r < 4; ++r) {
        const size_t ri = first_lone + static_cast<size_t>(r);
        ok = ok && chunk.haps[ri] == 0 && chunk.island_haps[ri] != 0 && chunk.island_phase_sets[ri] == 900 &&
             chunk.island_haps[ri] == chunk.island_haps[first_lone + static_cast<size_t>(r % 2)];
    }
    check(ok && chunk.island_haps[first_lone] != chunk.island_haps[first_lone + 1],
          "length islands: lone reads take their window side, in a block at the window");
}

int main() {
    test_em_phases_one_block();
    test_em_cuts_unlinked_blocks();
    test_em_leaves_noise_unphased();
    test_em_recovers_from_switched_start();
    test_em_labels_low_mapq_reads();
    test_em_window_links_blocks();
    test_inject_selects_private_sites();
    window_case(true);
    window_case(false);
    window_case(true, false);
    test_realign_recalls_indel();
    test_realign_recalls_merged_deletion();
    test_bridge_weak_sites();
    test_em_join_needs_strong_reads();
    test_msa_absence_of_alt();
    test_confirm_catalog_snp_calls();
    test_inject_island_only_sites();
    test_inject_low_mapq_clean_het();
    test_em_island_blocks();
    test_em_deferred_labels();
    test_em_echo_sites();
    test_em_unverified_fallback(false);
    test_em_unverified_fallback(true);
    test_length_window_sites();
    test_em_length_window_islands();
    if (failures != 0) {
        std::cerr << failures << " union phasing check(s) failed\n";
        return 1;
    }
    std::cout << "union phasing tests passed\n";
    return 0;
}
