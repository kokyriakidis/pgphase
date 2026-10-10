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
work=ROOT/'test_data/tmp_representation_step12/mutants';work.mkdir(parents=True,exist_ok=True)
source=(ROOT/'src/allele_context.cpp').read_text()
mutations={
 'publish_nonmandatory_matches':('consuming_edges == 1 && matched_alt < cols','consuming_edges >= 1 && matched_alt < cols'),
 'ignore_optimal_deletions':('== result.distance) ++consuming_edges;', '== result.distance) (void)consuming_edges;'),
 'allow_mismatch_edges_as_matches':("if (ref[i] == alt[j]) matched_alt = j;", "matched_alt = j;"),
 'shift_alt_offsets':('else result.matches.push_back({i, matched_alt, 1});', 'else result.matches.push_back({i, matched_alt + 1, 1});'),
 'break_subpath_coalescing':('++result.matches.back().length;', '(void)result.matches.back();'),
 'ignore_cell_budget':('rows * cols > max_cells', 'false && rows * cols > max_cells'),
 'lose_case_normalization':('base = context_base(base);', '(void)base;'),
 'reverse_insert_delete_costs':('back(i, j) = alt.size() - j;', 'back(i, j) = ref.size();'),
}

obj=work/'test.o'
subprocess.run(['g++','-O0','-std=c++17','-Isrc','-c','src/test_allele_context.cpp','-o',str(obj)],cwd=ROOT,check=True)
results={}
for name,(old,new) in mutations.items():
 assert source.count(old)>=1,name
 path=work/(name+'.cpp');path.write_text(source.replace(old,new,1));mutant=work/(name+'.o');binary=work/name
 subprocess.run(['g++','-O0','-std=c++17','-Isrc','-Ithird_party/edlib/edlib/include','-c',str(path),'-o',str(mutant)],cwd=ROOT,check=True)
 subprocess.run(['g++','-o',str(binary),str(obj),str(mutant),'src/allele_identity.o','src/edlib.o','-lhts'],cwd=ROOT,check=True)
 run=subprocess.run([str(binary)],cwd=ROOT,text=True,capture_output=True)
 (work/(name+'.log')).write_text(run.stdout+run.stderr)
 assert run.returncode!=0 and 'FAIL:' in run.stderr,name
 results[name]=dict(exit=run.returncode,failed_checks=run.stderr.count('FAIL:'))
 print(name,results[name],flush=True)
(OUT/'mutation-checks.json').write_text(json.dumps(results,indent=2)+'\n')
state=ROOT/'test_data/tmp_representation_step12/timing';state.mkdir(exist_ok=True)
commands={'fixtures':['./test_allele_context']}
for owner in ['matrix4-final','matrix5-final','matrix65-final','whole-final']:
 commands[owner]=['./test_allele_context','--maps',str(ROOT/'test_data/tmp_representation_step12/accepted-final'/owner),str(state/owner)]
timings={}
for name,command in commands.items():
 result=subprocess.run(command,cwd=ROOT,capture_output=True,text=True,check=True)
 start=time.monotonic()
 for _ in range(args.runs):subprocess.run(command,cwd=ROOT,stdout=subprocess.DEVNULL,check=True)
 timings[name]=dict(executions=args.runs,average_ms=1000*(time.monotonic()-start)/args.runs,result=result.stdout.strip())
 print(name,timings[name],flush=True)
(OUT/'fast-checks.json').write_text(json.dumps(dict(test_binary_sha256=hashlib.sha256((ROOT/'test_allele_context').read_bytes()).hexdigest(),timings=timings),indent=2)+'\n')
