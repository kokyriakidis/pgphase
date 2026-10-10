// Evaluation-only abPOA check; production discovery does not use these results.
#include "phasing_types.hpp"
#include "abpoa.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <vector>

static std::vector<std::string> fields(const std::string& line) {
    std::vector<std::string> result;
    std::istringstream input(line);
    for (std::string field; std::getline(input, field, '\t');) result.push_back(field);
    return result;
}

static void verify(const std::string& mode, const std::vector<std::pair<std::string,std::string>>& reads,
                   std::ostream& consensuses, std::ostream& membership) {
    const pgphase_collect::Options opts;
    std::unique_ptr<abpoa_t,decltype(&abpoa_free)> ab(abpoa_init(),abpoa_free);
    std::unique_ptr<abpoa_para_t,decltype(&abpoa_free_para)> para(abpoa_init_para(),abpoa_free_para);
    para->wb=-1; para->inc_path_score=1; para->out_msa=1; para->out_cons=1;
    para->cons_algrm=ABPOA_MF; para->max_n_cons=2; para->min_freq=opts.min_af;
    para->match=opts.match; para->mismatch=opts.mismatch;
    para->gap_open1=opts.gap_open1; para->gap_ext1=opts.gap_ext1;
    para->gap_open2=opts.gap_open2; para->gap_ext2=opts.gap_ext2;
    abpoa_post_set_para(para.get());
    std::vector<std::vector<uint8_t>> encoded;
    std::vector<int> lengths;
    for (const auto& read:reads) {
        std::vector<uint8_t> sequence;
        for(char base:read.second) {
            const size_t code=std::string("ACGT").find(base);
            if(code==std::string::npos) throw std::runtime_error("non-DNA input sequence");
            sequence.push_back(code);
        }
        encoded.push_back(std::move(sequence));lengths.push_back(read.second.size());
    }
    std::vector<uint8_t*> pointers;
    for(auto& sequence:encoded)pointers.push_back(sequence.data());
    abpoa_msa(ab.get(),para.get(),reads.size(),nullptr,lengths.data(),pointers.data(),nullptr,nullptr);
    const auto* result=ab->abc;
    if(result->n_cons<1) throw std::runtime_error("abPOA returned no consensus");
    for(size_t i=0;i<reads.size();++i) {
        std::string restored;
        for(int k=0;k<result->msa_len;++k)if(result->msa_base[i][k]<4) restored+="ACGT"[result->msa_base[i][k]];
        if(restored!=reads[i].second) throw std::runtime_error("MSA row does not preserve original read");
    }
    for(int ci=0;ci<result->n_cons;++ci) {
        std::string sequence;
        for(int k=0;k<result->cons_len[ci];++k)sequence+="ACGT"[result->cons_base[ci][k]];
        consensuses<<mode<<'\t'<<ci<<'\t'<<result->clu_n_seq[ci]<<'\t'<<sequence<<'\n';
        for(int k=0;k<result->clu_n_seq[ci];++k) {
            const int index=result->clu_read_ids[ci][k];
            if(index<0||index>=static_cast<int>(reads.size()))throw std::runtime_error("invalid cluster read index");
            membership<<mode<<'\t'<<ci<<'\t'<<reads[index].first<<'\n';
        }
    }
}

int main(int argc,char** argv) {
    if(argc!=4){std::cerr<<"Usage: msa_verify COHORT.tsv OUTPUT_PREFIX PHYSICAL_ID\n";return 1;}
    std::ifstream input(argv[1]);if(!input)throw std::runtime_error("cannot open cohort");
    std::vector<std::pair<std::string,std::string>> reads;
    std::string line;std::getline(input,line);
    while(std::getline(input,line)) {
        const auto row=fields(line);
        if(row[0]==argv[3])reads.emplace_back(row[1],row[5]);
    }
    if(reads.empty())throw std::runtime_error("selected physical context has no reads");
    std::ofstream consensuses(std::string(argv[2])+".consensuses.tsv"),membership(std::string(argv[2])+".membership.tsv");
    consensuses<<"mode\tconsensus\tmolecules\tsequence\n";membership<<"mode\tconsensus\tmolecule\n";
    verify("all",reads,consensuses,membership);
    std::reverse(reads.begin(),reads.end());verify("all-reversed",reads,consensuses,membership);
    std::map<std::string,std::vector<std::pair<std::string,std::string>>> exact;
    for(const auto& read:reads)exact[read.second].push_back(read);
    size_t group=0;
    for(const auto& [sequence,support]:exact)if(support.size()>=2) {
        verify("exact-"+std::to_string(group++),support,consensuses,membership);
    }
    auto longest=exact.end();
    for(auto it=exact.begin();it!=exact.end();++it)
        if(it->second.size()>=2 && (longest==exact.end() || it->first.size()>longest->first.size())) longest=it;
    if(longest!=exact.end())for(const auto& excluded:longest->second) {
        std::vector<std::pair<std::string,std::string>> heldout;
        for(const auto& read:reads)if(read.first!=excluded.first)heldout.push_back(read);
        verify("heldout:"+excluded.first,heldout,consensuses,membership);
    }
}
