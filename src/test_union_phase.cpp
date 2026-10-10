// Unit tests for union gap phasing: site injection, global EM, local windows.

#include "union_phase.hpp"

#include "collect_phase.hpp"

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
    add_tiled_reads(chunk, truth, 0, 7, 8);
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
    bam.candidates.push_back(clean_snp(300));                                    // 2: ALT only at low MAPQ
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
    const size_t changed = realign_indel_observations(chunk, bam);
    check(changed == 4, "realign: the four missing calls filled, got " + std::to_string(changed));
    bool ok = true;
    for (int r = 0; r < 8; ++r) ok = ok && chunk.read_var_profile[static_cast<size_t>(r)].alleles[0] == (r % 2 == 0 ? 1 : 0);
    check(ok, "realign: each read takes the allele its sequence carries");
}


static void test_deferred_labels_follow_anchor() {
    // A last-resort label is stored against an anchor read and applied after
    // stitching: it follows the anchor's final phase set and any flip.
    GraphChunkBuildResult gc;
    PhasingChunk& chunk = gc.chunk;
    for (int r = 0; r < 3; ++r) chunk.reads.push_back(read_named("d" + std::to_string(r), 60));
    finish_chunk(chunk);
    chunk.haps[0] = 1;
    chunk.phase_sets[0] = 100;
    gc.deferred_read_labels.emplace_back(1, 0, true);   // same haplotype as read 0
    gc.deferred_read_labels.emplace_back(2, 0, false);  // the other haplotype
    // Stitching flips the block and renames it.
    chunk.haps[0] = 2;
    chunk.phase_sets[0] = 50;
    check(apply_deferred_read_labels(gc) == 2, "deferred: both labels applied");
    check(chunk.haps[1] == 2 && chunk.phase_sets[1] == 50, "deferred: same-haplotype read follows the flip");
    check(chunk.haps[2] == 1 && chunk.phase_sets[2] == 50, "deferred: other-haplotype read follows the flip");
    check(gc.deferred_read_labels.empty(), "deferred: list consumed");
}


static void test_consensus_marker_labels() {
    // Four labelled reads per haplotype; haplotype 1 carries a SNP the callers
    // never reported. An unlabelled read carrying it is labelled haplotype 1,
    // one without it haplotype 2.
    std::string ref;
    uint32_t state = 4242;
    for (int i = 0; i < 600; ++i) {
        state = state * 1103515245u + 12345u;
        ref += "ACGT"[(state >> 16) % 4];
    }
    std::string alt_seq = ref;
    alt_seq[300] = ref[300] == 'A' ? 'C' : 'A';
    PhasingChunk bam;
    bam.ref_beg = 1;
    bam.ref_end = 600;
    bam.ref_seq = ref;
    GraphChunkBuildResult gc;
    PhasingChunk& chunk = gc.chunk;
    for (int r = 0; r < 10; ++r) {
        const std::string name = "m" + std::to_string(r);
        chunk.reads.push_back(read_named(name, 60));
        chunk.read_var_profile.push_back(profile(r, 0, {}));
        ReadRecord rec = read_named(name, 60);
        const bool carrier = r < 4 || r == 8;
        rec.alignment = make_alignment(name, carrier ? alt_seq : ref, {bam_cigar_gen(600, BAM_CMATCH)});
        bam.reads.push_back(std::move(rec));
    }
    finish_chunk(chunk);
    for (int r = 0; r < 8; ++r) {
        chunk.haps[static_cast<size_t>(r)] = r < 4 ? 1 : 2;
        chunk.phase_sets[static_cast<size_t>(r)] = 7;
    }
    const size_t decided = label_reads_from_haplotype_consensus(gc, bam);
    check(decided == 2, "markers: both unlabelled reads decided, got " + std::to_string(decided));
    check(apply_deferred_read_labels(gc) == 2, "markers: labels applied");
    check(chunk.haps[8] == 1 && chunk.haps[9] == 2 && chunk.phase_sets[8] == 7,
          "markers: carrier joins haplotype 1, non-carrier haplotype 2");
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
    test_deferred_labels_follow_anchor();
    test_consensus_marker_labels();
    test_bridge_weak_sites();
    if (failures != 0) {
        std::cerr << failures << " union phasing check(s) failed\n";
        return 1;
    }
    std::cout << "union phasing tests passed\n";
    return 0;
}
