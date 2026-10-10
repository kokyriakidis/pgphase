#!/usr/bin/env python3
"""Verify complete allele costs and held-out diploid fits against saved original evidence."""
import argparse
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--work', type=Path, default=ROOT/'test_data/tmp_representation_step7/cached')
parser.add_argument('--production', action='store_true', help='Read production matrix files and compare them to cached replay outputs')
parser.add_argument('--scalar-costs', action='store_true', help='Independently verify every individual allele distance with scalar DP')
args = parser.parse_args()

def rows(path):
    with path.open() as h:
        return list(csv.DictReader(h, delimiter='\t'))

def scalar_distance(a, b):
    row = list(range(len(b)+1))
    for i, x in enumerate(a):
        diagonal, row[0] = row[0], i+1
        for j, y in enumerate(b):
            old = row[j+1]
            row[j+1] = min(row[j]+1, row[j+1]+1, diagonal+(x!=y))
            diagonal = old
    return row[-1]

def fit(costs):
    if len(costs)==0:
        return None, None, None, 0
    n = costs.shape[1]
    pairs = [(a,b) for a in range(n) for b in range(a,n)]
    scores = np.array([np.minimum(costs[:,a],costs[:,b]).sum(dtype=np.int64) for a,b in pairs], dtype=np.int64)
    order = np.argsort(scores, kind='stable')
    best = int(scores[order[0]])
    ties = int(sum(scores==best))
    second = int(scores[order[1]]) if len(order)>1 else None
    return pairs[order[0]] if ties==1 else None, best, second, ties

def saved_fit(row):
    pair = None if row['allele0']=='.' else (int(row['allele0']),int(row['allele1']))
    return pair, None if row['cost']=='.' else int(row['cost']), None if row['runner_up']=='.' else int(row['runner_up']), int(row['tied_pairs'])

truth = {f[0]:f[1]=='PATERNAL' for line in (ROOT/'test_data/derived/chr20_truth_hap.tsv').open()
         if len(f:=line.rstrip().split('\t'))==2 and f[1] in ('PATERNAL','MATERNAL')}
_, hi_status = json.loads((ROOT/'test_data/tmp_representation_step7/hiphase-molecule-tags.json').read_text())
reports=[]
calibration_rows=[]
started=time.monotonic()
for name in ['matrix4-final','matrix5-final','matrix65-final','whole-final']:
    prefix = args.work/name/'matrix.chunk0.joint' if args.production else args.work/name
    def state(kind):
        return Path(str(prefix)+('-' if args.production else '.')+kind+'.tsv')
    hypotheses=defaultdict(list)
    for r in rows(state('alleles')):
        assert int(r['allele'])==len(hypotheses[r['locus']])
        hypotheses[r['locus']].append(r['sequence'])
    contexts={r['locus']:r for r in rows(ROOT/'test_data/tmp_representation_step6/accepted'/name/'matrix.chunk0.joint-contexts.tsv')}
    original={ (r['locus'],r['read']):r for r in rows(ROOT/'test_data/tmp_representation_step6/accepted'/name/'matrix.chunk0.joint-sequences.tsv') if r['status']=='scored'}
    cost_rows=rows(state('costs'))
    assert len(cost_rows)==len(original) and {(r['locus'],r['read']) for r in cost_rows}==set(original)
    grouped=defaultdict(list)
    scalar_checks=0
    for r in cost_rows:
        seqs=hypotheses[r['locus']]
        c=contexts[r['locus']]
        assert seqs==sorted({c['allele0'],c['allele1']}|set(c['other_alleles'].split(','))-{''})
        costs=list(map(int,r['costs'].split(',')))
        assert len(costs)==len(seqs) and all(v>=0 for v in costs)
        if args.scalar_costs:
            query=original[r['locus'],r['read']]['sequence']
            assert costs==[scalar_distance(query,s) for s in seqs]
            scalar_checks+=len(seqs)
        grouped[r['locus']].append((r['read'],costs))
    full={r['locus']:r for r in rows(state('genotypes'))}
    heldout={(r['locus'],r['read']):r for r in rows(state('heldout'))}
    assert set(full)==set(hypotheses) and set(heldout)==set(original)
    changed=ambiguous=stable_contexts=alt_alt=0
    partitions=Counter()
    conflicts=Counter()
    locus_reports=[]
    for locus, observations in grouped.items():
        names=[r[0] for r in observations]
        assert len(names)==len(set(names))
        matrix=np.asarray([r[1] for r in observations],dtype=np.int64)
        expected=fit(matrix)
        assert saved_fit(full[locus])==expected and int(full[locus]['molecules'])==len(names)
        assert full[locus]['conflicting_molecules']=='0'
        pair=expected[0]
        seqs=hypotheses[locus]
        old_pair=tuple(sorted([seqs.index(contexts[locus]['allele0']),seqs.index(contexts[locus]['allele1'])]))
        changed+=pair is not None and pair!=old_pair
        ambiguous+=pair is None
        alt_alt+=pair is not None and all(seqs[a]!=contexts[locus]['reference'] for a in pair)
        stable=0
        scored_predictions=[]
        for i,name_ in enumerate(names):
            training=np.delete(matrix,i,axis=0)
            actual=heldout[locus,name_]
            expected_heldout=fit(training)
            assert saved_fit(actual)==expected_heldout
            training_pair=expected_heldout[0]
            same=pair is not None and pair==training_pair
            winner=int(np.argmin(matrix[i]))
            call=winner if training_pair is not None and int(sum(matrix[i]==matrix[i,winner]))==1 and winner in training_pair else -1
            assert int(actual['allele'])==call and int(actual['same_pair'])==same
            stable+=same
            reason='ambiguous_pair' if training_pair is None else 'unstable_pair' if not same else 'outside_or_tied_sequence' if call<0 else 'stable_member'
            partitions[reason]+=1
            if original[locus,name_]['source_allele']=='-2':
                conflicts[reason]+=1
            # Parental truth is evaluation-only. Orient using other molecules;
            # the held-out read's truth never sets its own allele gauge.
            prediction=None
            if reason=='stable_member' and name_ in truth:
                votes={a:Counter() for a in training_pair}
                for j,other in enumerate(names):
                    if i==j or other not in truth:
                        continue
                    allele=int(np.argmin(matrix[j]))
                    if int(sum(matrix[j]==matrix[j,allele]))==1 and allele in votes:
                        votes[allele][truth[other]]+=1
                orientation={a:counts.most_common(1)[0][0] for a,counts in votes.items()
                             if counts and counts[True]!=counts[False]}
                if len(orientation)==2 and len(set(orientation.values()))==2:
                    prediction=orientation[call]
            scored_predictions.append((name_,prediction))
            source=original[locus,name_]
            quality_bytes=bytes.fromhex(source['qualities_hex'])
            ordered=sorted(matrix[i].tolist())
            outcome='unassigned' if prediction is None else 'correct' if prediction==truth[name_] else 'discordant'
            calibration_rows.append(dict(owner=name,locus=int(locus),read=name_,partition=reason,
                source_conflict=source['source_allele']=='-2',outcome=outcome,
                mapq=int(source['mapq']),minimum_base_quality=min(quality_bytes),
                unknown_base_qualities=sum(q==255 for q in quality_bytes),
                best_sequence_distance=ordered[0],nearest_margin=ordered[1]-ordered[0],
                genotype_cost_per_molecule=expected[1]/len(names)))
        stable_contexts+=pair is not None and stable==len(names)
        scorable=[(n,p) for n,p in scored_predictions if n in truth]
        correct=sum(p==truth[n] for n,p in scorable if p is not None)
        wrong=sum(p!=truth[n] for n,p in scorable if p is not None)
        hi_correct=sum(hi_status.get(n)=='correct' for n,p in scorable)
        locus_reports.append(dict(locus=int(locus),alleles=len(seqs),molecules=len(names),cost=expected[1],tied_pairs=expected[3],
            inferred_pair=pair,original_pair=old_pair,stable_heldout=stable,
            evaluation_only_parental_correct=correct,evaluation_only_parental_discordant=wrong,
            evaluation_only_parental_unassigned=len(scorable)-correct-wrong,
            identical_molecule_hiphase_correct=hi_correct))
    if args.production:
        for kind in ['alleles','costs','genotypes','heldout']:
            cached=ROOT/'test_data/tmp_representation_step7/cached'/f'{name}.{kind}.tsv'
            assert state(kind).read_bytes()==cached.read_bytes(), (name,kind)
    report=dict(owner=name,contexts=len(hypotheses),molecules=len(cost_rows),different_unique_pairs=changed,
                ambiguous_pairs=ambiguous,stable_for_every_heldout_molecule=stable_contexts,inferred_alt_alt_pairs=alt_alt,
                partitions=dict(partitions),conflict_partitions=dict(conflicts),scalar_distance_checks=scalar_checks,
                loci=locus_reports)
    reports.append(report)
    print(json.dumps({k:v for k,v in report.items() if k!='loci'}),flush=True)
file='production-genotype-checks.json' if args.production else 'cached-genotype-checks.json'
(OUT/file).write_text(json.dumps(dict(seconds=time.monotonic()-started,owners=reports),indent=2)+'\n')

if not args.production:
    bins=defaultdict(Counter)
    for r in calibration_rows:
        key=(r['owner'],r['minimum_base_quality'],r['nearest_margin'],r['partition'])
        bins[key][r['outcome']]+=1
    (OUT/'calibration-checks.json').write_text(json.dumps(dict(
        interpretation='Evaluation-only held-out parental separability; no production confidence threshold inferred',
        bins=[dict(owner=k[0],minimum_base_quality=k[1],nearest_margin=k[2],partition=k[3],counts=dict(v)) for k,v in sorted(bins.items())]
    ),indent=2)+'\n')
    with (OUT/'heldout-parental-audit.tsv').open('w') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(calibration_rows[0]),delimiter='\t')
        writer.writeheader();writer.writerows(calibration_rows)
