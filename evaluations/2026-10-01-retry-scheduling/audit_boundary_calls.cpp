#include "collect_phase_noisy.hpp"
#include "collect_types.hpp"
#include <iostream>
using namespace pgphase_collect;
int main() {
 Options opts; opts.bam_files={"test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam"}; opts.ref_fasta="test_data/chm13v2.0.chr20.renamed.fa";
 WorkerContext context(opts); const int tid=sam_hdr_name2tid(context.primary_header(),"CHM13#0#chr20");
 std::unique_ptr<hts_itr_t,decltype(&hts_itr_destroy)> it(sam_itr_queryi(context.indexes.front().get(),tid,57854342,57866714),&hts_itr_destroy);
 std::unique_ptr<bam1_t,AlignmentDeleter> read(bam_init1());
 std::cout<<"read\tmapq\texact1\tequivalent1\texact2\tequivalent2\tnear_cigar\n";
 while(sam_itr_next(context.bams.front()->get(),it.get(),read.get())>=0) {
  if(read->core.flag&(BAM_FUNMAP|BAM_FSECONDARY|BAM_FSUPPLEMENTARY|BAM_FDUP|BAM_FQCFAIL))continue;
  if(read->core.qual<30||read->core.qual==255||read->core.pos>57854341||bam_endpos(read.get())<57866714)continue;
  std::cout<<bam_get_qname(read.get())<<'\t'<<int(read->core.qual);
  for(int len=1;len<=2;++len) {CandidateVariant d;d.key.type=VariantType::Deletion;d.key.pos=57854342;d.key.ref_len=len;int qi=-1;
   std::cout<<'\t'<<bam_exact_indel_allele(read.get(),d,30,&qi)<<'\t'<<bam_equivalent_deletion_allele(read.get(),d,context.ref,tid,context.primary_header(),30);}
  hts_pos_t p=read->core.pos;int q=0;const auto* cigar=bam_get_cigar(read.get());std::cout<<'\t';
  for(uint32_t i=0;i<read->core.n_cigar;++i){int op=bam_cigar_op(cigar[i]),len=bam_cigar_oplen(cigar[i]);
   if(p+len>=57854302&&p<=57854382)std::cout<<p<<':'<<len<<bam_cigar_opchr(cigar[i])<<"@"<<q<<',';
   if(bam_cigar_type(op)&1)q+=len;if(bam_cigar_type(op)&2)p+=len;}
  std::cout<<'\n';
 }
}
