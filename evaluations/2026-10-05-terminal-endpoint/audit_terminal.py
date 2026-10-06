#!/usr/bin/env python3
"""Inspect original terminal alignments, physical alleles, and both tools' tags."""
from collections import Counter, defaultdict
import json
from pathlib import Path
import pysam

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
TRUTH = {f[0]: f[1] for line in (ROOT/'test_data/derived/chr20_truth_hap.tsv').open()
         if len(f := line.rstrip().split('\t')) == 2}
BAM = ROOT/'test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam'
PG = ROOT/'test_data/tmp_gap_fix69/frozen_final/0/phased.bam'
HI = ROOT/'test_data/tmp_gap_fix48/competitor/hiphase_dv/phased.bam'
CONTIG = 'CHM13#0#chr20'
LEFT, END = 66194204, 66206480
with pysam.AlignmentFile(str(BAM)) as bam:
    reads = [r for r in bam.fetch(CONTIG, LEFT-1, END)
             if not r.is_secondary and not r.is_supplementary]
names = {r.query_name for r in reads}
tags = {'pgphase': {}, 'hiphase': {}}
votes = defaultdict(Counter)
with pysam.AlignmentFile(str(PG), check_sq=False) as bam:
    for r in bam:
        if not r.is_secondary and not r.is_supplementary and r.query_name in TRUTH and r.has_tag('HP') and r.has_tag('PS'):
            hp, ps = r.get_tag('HP'), r.get_tag('PS')
            if hp in (1, 2) and ps > 0:
                votes[ps][(hp == 1) == (TRUTH[r.query_name] == 'MATERNAL')] += 1
        if r.query_name in names and not r.is_secondary and not r.is_supplementary:
            tags['pgphase'][r.query_name] = {t: r.get_tag(t) if r.has_tag(t) else 0 for t in ('HP','PS')}
original = {r.query_name: r for r in reads}
with pysam.AlignmentFile(str(HI)) as bam:
    contig = next(c for c in bam.references if c.endswith('chr20'))
    for r in bam.fetch(contig, LEFT-1, END):
        if r.query_name in names and not r.is_secondary and not r.is_supplementary:
            source = original[r.query_name]
            assert (r.reference_start, r.reference_end, r.cigarstring, r.query_sequence) == (source.reference_start, source.reference_end, source.cigarstring, source.query_sequence)
            tags['hiphase'][r.query_name] = {t: r.get_tag(t) if r.has_tag(t) else 0 for t in ('HP','PS')}
# Orient every pgphase PS on its full output cohort. HiPhase's only PS here
# is 63182011, independently oriented on its whole cohort in the preceding
# largest-block audit, where HP1 is maternal. Never orient on the target reads.
orientation = {ps: c.most_common(1)[0][0] for ps, c in votes.items()}
assert all(t['PS'] in (0, 63182011) for t in tags['hiphase'].values())
rows=[]
for r in reads:
    bases={p+1: (r.query_sequence[q], r.query_qualities[q])
           for q,p in r.get_aligned_pairs() if q is not None and p is not None and p+1 in (LEFT, END)}
    row={'qname':r.query_name,'truth':TRUTH.get(r.query_name),'mapq':r.mapping_quality,
         'start':r.reference_start+1,'end':r.reference_end,'bases':bases}
    for tool in tags:
        t=tags[tool].get(r.query_name, {'HP':0,'PS':0})
        row[tool]=t
        row[tool]['status']='unphased' if t['HP'] not in (1,2) or not t['PS'] else (
            'correct' if ((t['HP']==1)==(row['truth']=='MATERNAL')) == (orientation[t['PS']] if tool == 'pgphase' else True) else 'discordant')
    rows.append(row)
summary={}
for name,beg,end in [('tail',LEFT,END),('terminal_site',END,END)]:
    eligible=[r for r in rows if r['truth'] and r['start']<=end and r['end']>=beg]
    cores = {tool: Counter(r[tool]['PS'] for r in eligible if r[tool]['status']=='correct' and r[tool]['PS']<1_000_000_000) for tool in tags}
    summary[name]={'scorable':len(eligible), 'parental_counts': dict(Counter(r['truth'] for r in eligible)),
                   'dominant_correct_core': {tool: max(c.values(), default=0) for tool, c in cores.items()}, 'pgphase':dict(Counter(r['pgphase']['status'] for r in eligible)),
                   'hiphase':dict(Counter(r['hiphase']['status'] for r in eligible)),
                   'physical_terminal':dict(Counter(str((r['truth'],r['bases'].get(END, ('no_base',0))[0])) for r in eligible if name=='terminal_site'))}
result={'pgphase_full_cohort_orientation':{str(ps): dict(votes[ps]) for ps in {r['pgphase']['PS'] for r in rows}}, 'hiphase_global_hp1_parent':'MATERNAL','interval':[LEFT,END], 'summary':summary,'reads':rows}
(OUT/'terminal-reads.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(summary,indent=2))
# Preserve the actual catalog alleles and HiPhase records in the tail.
for label,path,contig in [('catalog','test_data/chr20.sites.striped.vcf.gz',CONTIG),
                          ('hiphase','test_data/tmp_gap_fix48/competitor/hiphase_dv/phased.vcf.gz','chr20')]:
    with pysam.VariantFile(str(ROOT/path)) as vcf:
        records=[str(v) for v in vcf.fetch(contig,LEFT-1,END)]
    selected = records if label == 'hiphase' else [s for s in records if int(s.split('\t')[1]) == END]
    (OUT/(label+'-terminal.vcf.txt' if label == 'catalog' else label+'-tail.vcf.txt')).write_text(''.join(selected))
    print(label, len(records), 'tail records; exact terminal:', ''.join(s for s in records if int(s.split('\t')[1])==END).strip())
