#!/usr/bin/env python3
"""Compare the native owner and verify one/four-thread output identity."""
from collections import Counter
import hashlib
import importlib.util
import json
from pathlib import Path
import pysam

ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('audit',ROOT/'evaluations/2026-10-05-largest-hiphase-block/audit_reads.py')
audit=importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)
truth={f[0]:f[1]=='PATERNAL' for line in (ROOT/'test_data/derived/chr20_truth_hap.tsv').open()
       if len(f:=line.rstrip().split('\t'))==2}
base=ROOT/'test_data/tmp_gap_fix74/baseline/52'
fixed=ROOT/'test_data/tmp_gap_fix74/frozen_final/52'
single=ROOT/'test_data/tmp_gap_fix74/frozen_single/52'
result={'preservation':audit.audit(base,fixed,truth),
        'binary_sha256':hashlib.sha256((ROOT/'pgphase').read_bytes()).hexdigest()}
for name in ('candidates.tsv','phased.vcf'):
    assert (fixed/name).read_bytes()==(single/name).read_bytes(),name
assert audit.assignments(fixed/'phased.bam',truth)==audit.assignments(single/'phased.bam',truth)
result['one_four_thread_identity']=True
with pysam.AlignmentFile(str(ROOT/'test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam')) as bam:
    names={r.query_name for r in bam.fetch('CHM13#0#chr20',52696939,52711825)
           if not r.is_secondary and not r.is_supplementary and r.query_name in truth}
result['seam']={}
for label,directory in [('before',base),('pgphase',fixed)]:
    tags,status=audit.assignments(directory/'phased.bam',truth)
    counts=Counter(status.get(q,'unphased') for q in names)
    core=Counter(tags[q][1] for q in names if status.get(q)=='correct' and tags[q][1]<1000000000)
    result['seam'][label]={'scorable':len(names),'counts':dict(counts),'core_correct':max(core.values(),default=0)}
result['hiphase_seam']={'scorable':136,'correct':117,'core_correct':117,
    'geometry_verified_source':'evaluations/2026-10-05-fourth-largest-block-target/next-block.json'}
assert result['seam']['pgphase']['scorable']==136
assert result['seam']['pgphase']['counts']['correct']>=117
assert result['seam']['pgphase']['core_correct']>=117
(OUT/'owner-results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
