#!/usr/bin/env python3
"""Run six preselected HiPhase-bridged gaps with fixed recovery settings."""
import csv, hashlib, json, shlex, subprocess
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
OUT = Path('/tmp/pgphase-singleton-link-validation')
DATA = Path.home() / 'Downloads/pgphase-eval-data'
REPORT = Path(__file__).resolve().parent
OUT.mkdir(exist_ok=True)
(REPORT / 'fingerprint.txt').write_text(hashlib.sha256((ROOT/'pgphase').read_bytes()).hexdigest()+'\n')
manifest=[]
for chrom in ('chr12','chr18','chr20'):
    rows=list(csv.DictReader((ROOT/f'evaluations/2026-09-12-chr12-18-20-comparison/{chrom}.correct_bridges.tsv').open(), delimiter='\t'))
    rows=[r for r in rows if r['tool'].lower()=='hiphase']
    if chrom=='chr20': rows=rows[2:]
    for row in rows[:2]:
        manifest.append(dict(row, name=f"{chrom}_{row['gap_start']}", region=f"{'CHM13#0#' if chrom=='chr20' else ''}{chrom}:{int(row['gap_start'])-50000}-{int(row['gap_end'])+50000}"))
(REPORT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
def run(cmd, log):
    with open(log,'w') as f: subprocess.run([str(x) for x in cmd],stdout=f,stderr=subprocess.STDOUT,check=True)
for row in manifest:
    chrom=row['chromosome']; case=OUT/row['name']; case.mkdir(exist_ok=True)
    print(row['name'], flush=True)
    for arm in ('clean','auto'):
        sites=f'{chrom}.sites.striped.vcf.gz' if chrom=='chr20' else f'{chrom}.sites.vcf.gz'
        cmd=[ROOT/'pgphase','collect-hybrid-variation','--ref',ROOT/f'test_data/chm13v2.0.{chrom}.renamed.fa','--bam',ROOT/f'test_data/HG002_{chrom}_hifi_mapped_to_CHM13_{chrom}_annotated.bam','--graph-sites',ROOT/'test_data'/sites,'--gaf',ROOT/f'test_data/HG002.{chrom}.annotated.coord.gaf.gz','-r',row['region'],'--chunk-size','500000','--threads','2','--link-by-alleles','--block-link-window','8','--min-read-margin','2','-o',case/f'{arm}.tsv','--phased-vcf-out',case/f'{arm}.vcf','-b',case/f'{arm}.bam']
        if arm=='auto': cmd+=['--min-block-link-reads','1']
        if arm=='auto': cmd+=['--phase-matrix-dump',str(case/'matrix')]
        if arm=='auto': cmd+=['--recover-gaps','--gap-recovery-report',case/'auto.tiers.tsv']
        (case/f'{arm}.command.sh').write_text(shlex.join(map(str,cmd))+'\n')
        run(cmd,case/f'{arm}.log')
    import pysam
    names=set()
    for arm in ('clean','auto'):
        with pysam.AlignmentFile(case/f'{arm}.bam') as bam: names.update(r.query_name for r in bam)
    (case/'names.txt').write_text('\n'.join(sorted(names))+'\n')
    with (case/'truth.bam').open('wb') as f:
        subprocess.run(['samtools','view','-b','-N',str(case/'names.txt'),str(DATA/f'truth/{chrom}/diplinator_merged.bam')],stdout=f,check=True)
    for arm in ('clean','auto'):
        dest=case/f'{arm}.eval'; dest.mkdir(exist_ok=True)
        run(['python3',ROOT/'scripts/evaluate_phase_accuracy.py',case/f'{arm}.bam',case/'truth.bam','0','0','5','',dest,'samtools','','','',''],case/f'{arm}.eval.log')
    print('evaluated '+row['name'],flush=True)
