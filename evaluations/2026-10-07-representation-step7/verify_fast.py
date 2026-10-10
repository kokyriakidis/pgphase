#!/usr/bin/env python3
"""Verify production genotype mutations and measure cached full-allele fits."""
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
work=ROOT/'test_data/tmp_representation_step7/mutants';work.mkdir(parents=True,exist_ok=True)
source=(ROOT/'src/allele_genotype.cpp').read_text()
mutations={
 'retain_first_tied_pair':('fit.alleles.reset();','(void)pair;'),
 'train_on_heldout_read':('total - own_cost','total'),
 'count_duplicate_molecules':('by_molecule.emplace(observation.molecule, observation.costs)','by_molecule.emplace(observation.molecule + ":" + std::to_string(by_molecule.size()), observation.costs)'),
 'overwrite_conflicting_molecule':('entry->second.reset();','entry->second = observation.costs;'),
 'force_heterozygous_pairs':('size_t second = first;','size_t second = first + 1;'),
 'overflow_cohort_cost':('int64_t total = 0;','int total = 0;'),
 'ignore_other_read_hypotheses':('const auto best = std::min_element(row.begin(), row.end());',
  'const auto best = row.begin() + (*heldout.fit.alleles)[row[(*heldout.fit.alleles)[0]] < row[(*heldout.fit.alleles)[1]] ? 0 : 1];'),
}
test_obj=work/'test.o'
subprocess.run(['g++','-O0','-std=c++17','-Isrc','-c','src/test_allele_genotype.cpp','-o',str(test_obj)],cwd=ROOT,check=True)
results={}
for name,(old,new) in mutations.items():
 assert source.count(old)==1,name
 path=work/(name+'.cpp');path.write_text(source.replace(old,new))
 obj,binary=work/(name+'.o'),work/name
 subprocess.run(['g++','-O0','-std=c++17','-Isrc','-c',str(path),'-o',str(obj)],cwd=ROOT,check=True)
 subprocess.run(['g++','-o',str(binary),str(test_obj),str(obj),'src/allele_context.o','src/allele_identity.o','src/edlib.o','-lhts'],cwd=ROOT,check=True)
 run=subprocess.run([str(binary)],cwd=ROOT,text=True,capture_output=True)
 (work/(name+'.log')).write_text(run.stdout+run.stderr)
 assert run.returncode!=0 and 'FAIL:' in run.stderr,name
 results[name]=dict(exit=run.returncode,failed_checks=run.stderr.count('FAIL:'))
(OUT/'mutation-checks.json').write_text(json.dumps(results,indent=2)+'\n')
commands={'fixtures':['./test_allele_genotype']}
state=ROOT/'test_data/tmp_representation_step7/timing';state.mkdir(exist_ok=True)
for name in ['matrix4-final','matrix5-final','matrix65-final','whole-final']:
 commands[name]=['./test_allele_genotype','--state',str(ROOT/'test_data/tmp_representation_step6/accepted'/name),str(state/name)]
timings={}
for name,command in commands.items():
 result=subprocess.run(command,cwd=ROOT,capture_output=True,text=True,check=True)
 start=time.monotonic()
 for _ in range(args.runs):subprocess.run(command,cwd=ROOT,stdout=subprocess.DEVNULL,check=True)
 timings[name]=dict(executions=args.runs,average_ms=1000*(time.monotonic()-start)/args.runs,result=result.stdout.strip())
(OUT/'fast-checks.json').write_text(json.dumps(dict(test_binary_sha256=hashlib.sha256((ROOT/'test_allele_genotype').read_bytes()).hexdigest(),timings=timings),indent=2)+'\n')
print(json.dumps(dict(mutations=results,timings=timings),indent=2))
