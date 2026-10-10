#!/usr/bin/env python3
"""Validate independent whole-cohort MSA agreement and discovery provenance."""
from collections import defaultdict
import csv
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent
WORK=ROOT/'test_data/tmp_representation_step14'
def rows(path):
    with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))
folder=WORK/'accepted-final/matrix65-final'
cohort={r['read']:r['sequence'] for r in rows(folder/'matrix.chunk0.physical-cohort.tsv') if r['physical']=='3'}
assert len(cohort)==39
exact=defaultdict(set)
for name,sequence in cohort.items():exact[sequence].add(name)
supported={s:names for s,names in exact.items() if len(names)>=2}
long=max(supported,key=len);short=max(supported,key=lambda s:len(supported[s]))
assert len(long)==312 and len(supported[long])==6 and len(short)==51 and len(supported[short])==12
consensuses=defaultdict(dict);membership=defaultdict(set)
for r in rows(WORK/'cached/msa-context3.membership.tsv'):membership[r['mode'],r['consensus']].add(r['molecule'])
for r in rows(WORK/'cached/msa-context3.consensuses.tsv'):
    names=membership[r['mode'],r['consensus']]
    assert len(names)==int(r['molecules'])
    consensuses[r['mode']][r['sequence']]=names
assert set(consensuses['all'])=={short,long}
assert consensuses['all']==consensuses['all-reversed']
assert len(consensuses['all'][short])==23 and len(consensuses['all'][long])==16
assert set.union(*consensuses['all'].values())==set(cohort)
assert not (consensuses['all'][short]&consensuses['all'][long])
assert supported[long]<=consensuses['all'][long] and supported[short]<=consensuses['all'][short]
exact_modes=[m for m in consensuses if m.startswith('exact-')]
assert len(exact_modes)==len(supported)==4
assert {s for m in exact_modes for s in consensuses[m]}==set(supported)
for mode in exact_modes:
    assert len(consensuses[mode])==1
    sequence,names=next(iter(consensuses[mode].items()))
    assert names==supported[sequence]
heldout_modes={m for m in consensuses if m.startswith('heldout:')}
assert heldout_modes=={'heldout:'+name for name in supported[long]}
for mode in heldout_modes:
    excluded=mode.removeprefix('heldout:')
    assert consensuses[mode]=={s:names-{excluded} for s,names in consensuses['all'].items()}
original=[r['sequence'] for r in rows(folder/'matrix.chunk0.nested-alleles.tsv') if r['physical']=='3']
assert max(map(len,original))==55 and long not in original and short not in original
actual={r['sequence'] for r in rows(folder/'matrix.chunk0.read-hypotheses.tsv') if r['physical']=='3'}
assert actual==set(supported)
truth={f[0]:f[1] for line in (ROOT/'test_data/derived/chr20_truth_hap.tsv').read_text().splitlines() if len(f:=line.split('\t'))==2}
def parental(names):
    return {label:sum(truth.get(name)==label for name in names) for label in ['MATERNAL','PATERNAL']}
exact_parental={str(len(s)):parental(names) for s,names in supported.items()}
msa_parental={str(len(s)):parental(names) for s,names in consensuses['all'].items()}
assert exact_parental['312']=={'MATERNAL':6,'PATERNAL':0}
assert exact_parental['51']=={'MATERNAL':0,'PATERNAL':12}
assert msa_parental['312']=={'MATERNAL':11,'PATERNAL':5}
assert msa_parental['51']=={'MATERNAL':4,'PATERNAL':19}
report=dict(owner='matrix65-final',physical=3,original_molecules=39,seed_max_length=55,long_consensus_length=312,short_consensus_length=51,exact_long_supporters=6,exact_short_supporters=12,msa_long_cluster_molecules=16,msa_short_cluster_molecules=23,reversed_order_consensuses_and_membership_identical=True,long_supporter_exclusions_verified=6,all_msa_input_rows_preserved=True,exact_group_consensuses_verified=4,old_catalog_missing_both_msa_consensuses=True,read_hypotheses_match_independent_grouping=True,production_uses_msa_results=False)
report.update(posthoc_exact_support_parental_counts=exact_parental,posthoc_msa_cluster_parental_counts=msa_parental,truth_used_for_discovery=False)
(OUT/'msa-checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
