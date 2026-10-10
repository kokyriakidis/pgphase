#!/usr/bin/env python3
"""Kill production cohort admission mutants with the focused context fixtures."""
import json
from pathlib import Path
import subprocess

ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent
WORK=ROOT/'test_data/tmp_representation_step8/cohort-mutants'
WORK.mkdir(parents=True,exist_ok=True)
source=(ROOT/'src/allele_context.cpp').read_text()
mutations={
    'admit_secondary':('BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY','BAM_FUNMAP | BAM_FSUPPLEMENTARY'),
    'ignore_mapq_floor':('alignment->core.qual < min_mapq','false'),
    'ignore_filtered_policy':('!include_filtered && (alignment->core.flag & (BAM_FQCFAIL | BAM_FDUP)) != 0','false'),
    'revive_duplicate_primary':('if (!inserted) entry->second = nullptr;','if (!inserted) entry->second = alignment;'),
    'coverage_before_uniqueness':('const auto [entry, inserted] = unique.emplace',
        'if (alignment->core.pos + 1 > context.beg || bam_endpos(alignment) < context.end) continue;\n        const auto [entry, inserted] = unique.emplace'),
    'admit_reference_skip':('if (op == BAM_CREF_SKIP && ref_pos <= context.end && ref_pos + length > context.beg)', 'if (false)'),
}
obj=WORK/'test.o'
subprocess.run(['g++','-O0','-std=c++17','-Isrc','-c','src/test_allele_context.cpp','-o',str(obj)],cwd=ROOT,check=True)
results={}
for name,(old,new) in mutations.items():
    assert source.count(old)==1,name
    path=WORK/(name+'.cpp');path.write_text(source.replace(old,new))
    mutant=WORK/(name+'.o');binary=WORK/name
    subprocess.run(['g++','-O0','-std=c++17','-Isrc','-Ithird_party/edlib/edlib/include','-c',str(path),'-o',str(mutant)],cwd=ROOT,check=True)
    subprocess.run(['g++','-o',str(binary),str(obj),str(mutant),'src/allele_identity.o','src/edlib.o','-lhts'],cwd=ROOT,check=True)
    result=subprocess.run([str(binary)],cwd=ROOT,text=True,capture_output=True)
    (WORK/(name+'.log')).write_text(result.stdout+result.stderr)
    assert result.returncode!=0 and 'FAIL:' in result.stderr,name
    results[name]=dict(exit=result.returncode,failed_checks=result.stderr.count('FAIL:'))
(OUT/'cohort-mutation-checks.json').write_text(json.dumps(results,indent=2)+'\n')
print(json.dumps(results,indent=2))
