#include "collect_phase_noisy.hpp"
#include "collect_types.hpp"
#include "graph_bam_adapter.hpp"
#include <iostream>
using namespace pgphase_collect;
int main(int argc, char** argv) {
    if (argc != 7) return 2;
    Options opts;
    opts.bam_files = {argv[1]};
    opts.ref_fasta = argv[2];
    WorkerContext context(opts);
    const int tid = sam_hdr_name2tid(context.primary_header(), argv[3]);
    const hts_pos_t left_pos = std::stoll(argv[4]);
    const hts_pos_t right_pos = std::stoll(argv[5]);
    CandidateVariant deletion;
    deletion.key.type = VariantType::Deletion;
    deletion.key.pos = std::stoll(argv[6]);
    deletion.key.ref_len = 1;
    std::unique_ptr<hts_itr_t, decltype(&hts_itr_destroy)> iterator(
        sam_itr_queryi(context.indexes.front().get(), tid, left_pos - 1, right_pos),
        &hts_itr_destroy);
    std::unique_ptr<bam1_t, AlignmentDeleter> read(bam_init1());
    std::cout << "qname\tmapq\tbeg\tend\tleft_snp\tleft_quality\tdeletion\tright_snp\tright_quality\n";
    while (sam_itr_next(context.bams.front()->get(), iterator.get(), read.get()) >= 0) {
        if (read->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FDUP | BAM_FQCFAIL)) continue;
        if (read->core.qual < 30 || read->core.qual == 255) continue;
        int left_quality = 0, right_quality = 0;
        const int left = physical_snp_call(read.get(), left_pos, 'G', 'C', &left_quality);
        const int right = physical_snp_call(read.get(), right_pos, 'G', 'A', &right_quality);
        const int allele = bam_equivalent_deletion_allele(
            read.get(), deletion, context.ref, tid, context.primary_header(), 30);
        if ((left != 0 && left != 2) || left_quality < 10 || left_quality == 255) continue;
        std::cout << bam_get_qname(read.get()) << '\t' << static_cast<int>(read->core.qual)
            << '\t' << read->core.pos << '\t' << bam_endpos(read.get()) << '\t'
            << left << '\t' << left_quality << '\t' << allele << '\t' << right
            << '\t' << right_quality << '\n';
    }
}
