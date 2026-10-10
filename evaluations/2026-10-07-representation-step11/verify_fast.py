#!/usr/bin/env python3
"""Detect production-composer mutations and time saved-state compositions."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time

ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--runs',type=int,default=100)
args=parser.parse_args()
work=ROOT/'test_data/tmp_representation_step11/mutants';work.mkdir(parents=True,exist_ok=True)
source=(ROOT/'src/allele_context.cpp').read_text()
mutations={
 'omit_no_edit_choices':('        visit(parent, next + 1);\n        const auto& group', '        const auto& group'),
 'retain_only_first_alt':('alt < site.alts.size()', 'alt < std::min(size_t(1), site.alts.size())'),
 'prune_unsupported_length_prefix':('composed.status == AlleleCompositionStatus::Overlap', 'composed.status != AlleleCompositionStatus::Valid'),
 'expose_truncated_catalog':('result.paths.clear();', '(void)result.paths;'),
 'ignore_prefix_limit':('result.visited_prefixes == max_prefixes', 'false'),
 'discard_overlap_pruning':('composed.status == AlleleCompositionStatus::Overlap', 'false'),
 'collapse_parent_provenance':('result.paths.push_back({parent, std::move(provenance), composed.sequence});', 'result.paths.push_back({0, std::move(provenance), composed.sequence});'),
 'forget_selected_backtracking':('selected.resize(previous_selected);', '(void)selected;'),
 'ignore_recursion_depth_limit':('groups.size() > kMaxAllelePathSites', 'false'),
 'keep_duplicate_physical_choices':('std::make_tuple(site.pos, site.ref, alts)', 'std::make_tuple(site.pos + static_cast<hts_pos_t>(site.candidate), site.ref, alts)'),
 'ignore_full_alt_table':('std::make_tuple(site.pos, site.ref, alts)', 'std::make_tuple(site.pos, site.ref, std::vector<std::string>{})'),
 'use_input_neighbor_order':('std::sort(sites.begin(), sites.end(), [](const auto& left, const auto& right) {\n        return left.candidate < right.candidate;\n    });', '(void)sites;'),
}
obj=work/'test.o'
subprocess.run(['g++','-O0','-std=c++17','-Isrc','-c','src/test_allele_context.cpp','-o',str(obj)],cwd=ROOT,check=True)
results={}
for name,(old,new) in mutations.items():
 assert source.count(old)==1,name
 path=work/(name+'.cpp');path.write_text(source.replace(old,new));mutant=work/(name+'.o');binary=work/name
 subprocess.run(['g++','-O0','-std=c++17','-Isrc','-Ithird_party/edlib/edlib/include','-c',str(path),'-o',str(mutant)],cwd=ROOT,check=True)
 subprocess.run(['g++','-o',str(binary),str(obj),str(mutant),'src/allele_identity.o','src/edlib.o','-lhts'],cwd=ROOT,check=True)
 run=subprocess.run([str(binary)],cwd=ROOT,text=True,capture_output=True)
 (work/(name+'.log')).write_text(run.stdout+run.stderr)
 assert run.returncode!=0 and 'FAIL:' in run.stderr,name
 results[name]=dict(exit=run.returncode,failed_checks=run.stderr.count('FAIL:'))
 print(name,results[name],flush=True)
(OUT/'mutation-checks.json').write_text(json.dumps(results,indent=2)+'\n')
state=ROOT/'test_data/tmp_representation_step11/timing';state.mkdir(exist_ok=True)
commands={'fixtures':['./test_allele_context']}
for owner in ['matrix4-final','matrix5-final','matrix65-final','whole-final']:
 commands[owner]=['./test_allele_context','--paths',str(ROOT/'test_data/tmp_representation_step11/accepted-final'/owner),str(state/owner)]
timings={}
for name,command in commands.items():
 result=subprocess.run(command,cwd=ROOT,capture_output=True,text=True,check=True)
 start=time.monotonic()
 for _ in range(args.runs):subprocess.run(command,cwd=ROOT,stdout=subprocess.DEVNULL,check=True)
 timings[name]=dict(executions=args.runs,average_ms=1000*(time.monotonic()-start)/args.runs,result=result.stdout.strip())
 print(name,timings[name],flush=True)
(OUT/'fast-checks.json').write_text(json.dumps(dict(test_binary_sha256=hashlib.sha256((ROOT/'test_allele_context').read_bytes()).hexdigest(),timings=timings),indent=2)+'\n')
