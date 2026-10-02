#include "collect_phase_noisy.hpp"
#include "collect_types.hpp"
#include <array>
#include <iostream>
using namespace pgphase_collect;
int main() {
    Options opts;
    opts.bam_files = {"test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam"};
    opts.ref_fasta = "test_data/chm13v2.0.chr20.renamed.fa";
    WorkerContext context(opts);
    const int tid = sam_hdr_name2tid(context.primary_header(), "CHM13#0#chr20");
    std::unique_ptr<hts_itr_t, decltype(&hts_itr_destroy)> iterator(
        sam_itr_queryi(context.indexes.front().get(), tid, 50548244, 50562066), &hts_itr_destroy);
    std::unique_ptr<bam1_t, AlignmentDeleter> read(bam_init1());
    std::array<std::array<int, 3>, 2> old_counts{}, new_counts{};
    int spanning = 0;
    while (sam_itr_next(context.bams.front()->get(), iterator.get(), read.get()) >= 0) {
        if (read->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FDUP | BAM_FQCFAIL)) continue;
        if (read->core.qual < 30 || read->core.qual == 255 || read->core.pos > 50548244 || bam_endpos(read.get()) < 50562066) continue;
        ++spanning;
        for (int length = 1; length <= 2; ++length) {
            CandidateVariant deletion;
            deletion.key.type = VariantType::Deletion;
            deletion.key.pos = 50548246;
            deletion.key.ref_len = length;
            int qi = -1;
            const int old_call = bam_exact_indel_allele(read.get(), deletion, 30, &qi);
            const int new_call = bam_equivalent_deletion_allele(read.get(), deletion, context.ref, tid, context.primary_header(), 30);
            ++old_counts[length - 1][old_call + 1];
            ++new_counts[length - 1][new_call + 1];
        }
    }
    std::cout << "length\tspanning\texact_unknown\texact_ref\texact_alt\tequivalent_unknown\tequivalent_ref\tequivalent_alt\n";
    for (int length = 1; length <= 2; ++length) {
        std::cout << length << '\t' << spanning;
        for (int n : old_counts[length - 1]) std::cout << '\t' << n;
        for (int n : new_counts[length - 1]) std::cout << '\t' << n;
        std::cout << '\n';
    }
}
