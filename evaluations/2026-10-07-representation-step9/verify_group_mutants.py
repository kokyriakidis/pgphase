#!/usr/bin/env python3
"""Require the grouping fixtures to detect mutations of production identity."""
import json
from pathlib import Path
import subprocess

ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent
WORK=ROOT/'test_data/tmp_representation_step9/group-mutants';WORK.mkdir(parents=True,exist_ok=True)
source=(ROOT/'src/allele_context.cpp').read_text()
mutations={
 'ignore_left_bound':('std::tie(beg, end, alleles) < std::tie(other.beg, other.end, other.alleles)', 'std::tie(end, alleles) < std::tie(other.end, other.alleles)'),
 'ignore_right_bound':('std::tie(beg, end, alleles) < std::tie(other.beg, other.end, other.alleles)', 'std::tie(beg, alleles) < std::tie(other.beg, other.alleles)'),
 'ignore_other_hypotheses':('std::tie(beg, end, alleles) < std::tie(other.beg, other.end, other.alleles)', 'std::tie(beg, end) < std::tie(other.beg, other.end)'),
 'use_selected_pair_as_identity':('full_allele_sequences(context)}].push_back(i)', 'std::vector<std::string>(context.alleles.begin(), context.alleles.end())}].push_back(i)'),
 'retain_only_first_description':('groups[{context.beg, context.end, full_allele_sequences(context)}].push_back(i);',
  'auto& members = groups[{context.beg, context.end, full_allele_sequences(context)}]; if (members.empty()) members.push_back(i);'),
 'make_each_description_independent':('groups[{context.beg, context.end, full_allele_sequences(context)}]', 'groups[{context.beg + static_cast<hts_pos_t>(i), context.end, full_allele_sequences(context)}]'),
}
obj=WORK/'test.o'
subprocess.run(['g++','-O0','-std=c++17','-Isrc','-c','src/test_allele_context.cpp','-o',str(obj)],cwd=ROOT,check=True)
results={}
for name,(old,new) in mutations.items():
 assert source.count(old)==1,name
 path=WORK/(name+'.cpp');path.write_text(source.replace(old,new));mutant=WORK/(name+'.o');binary=WORK/name
 subprocess.run(['g++','-O0','-std=c++17','-Isrc','-Ithird_party/edlib/edlib/include','-c',str(path),'-o',str(mutant)],cwd=ROOT,check=True)
 subprocess.run(['g++','-o',str(binary),str(obj),str(mutant),'src/allele_identity.o','src/edlib.o','-lhts'],cwd=ROOT,check=True)
 r=subprocess.run([str(binary)],cwd=ROOT,text=True,capture_output=True)
 (WORK/(name+'.log')).write_text(r.stdout+r.stderr)
 assert r.returncode!=0 and 'FAIL:' in r.stderr,name
 results[name]=dict(exit=r.returncode,failed_checks=r.stderr.count('FAIL:'))
(OUT/'group-mutation-checks.json').write_text(json.dumps(results,indent=2)+'\n');print(json.dumps(results,indent=2))
