#!/usr/bin/env python3
"""Retest every truth-concordant, multiply-spanned clean-panel break locally.

Local runs diagnose a failure class; they do not replace full-chromosome results.
"""
import argparse,csv,hashlib,json,shlex,shutil,subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import pysam
ROOT=Path(__file__).resolve().parents[2]
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--chromosomes',nargs='+',default=['chr20'])
parser.add_argument('--label',default='chr20_link_panel')
parser.add_argument('--workers',type=int,default=2)
parser.add_argument('--all-concordant',action='store_true',help='Include singleton and indirect links too')
parser.add_argument('--cases',nargs='+',help='Restrict to named chromosome_left_right cases')
a=parser.parse_args()
REPORT=Path(__file__).resolve().parent/a.label
OUT=Path('/tmp')/('pgphase-'+a.label)
DATA=Path.home()/'Downloads/pgphase-eval-data'
REPORT.mkdir(exist_ok=True);OUT.mkdir(exist_ok=True)
binary=OUT/'pgphase'
if not binary.exists():shutil.copy2(ROOT/'pgphase',binary)
(REPORT/'binary.sha256').write_text(hashlib.sha256(binary.read_bytes()).hexdigest()+'\n')
manifest={}
for chrom in a.chromosomes:
    path=REPORT.parent/f'{chrom}.clean.missed.tsv'
    if not path.exists():raise FileNotFoundError(path)
    for r in csv.DictReader(path.open(),delimiter='\t'):
        if r['kind']!='split_block' or r['competitor_pair_truth']!='concordant' or (not a.all_concordant and r['reason']!='multiple_mapq30_spanning_reads_need_link_audit'):continue
        left,right=int(r['left']),int(r['right'])
        manifest[chrom,left,right]={'chromosome':chrom,'left':left,'right':right,'name':f'{chrom}_{left}_{right}','region':f"{'CHM13#0#' if chrom=='chr20' else ''}{chrom}:{max(1,left-50000)}-{right+50000}"}
manifest=list(manifest.values())
if a.cases:
    manifest=[row for row in manifest if row['name'] in a.cases]
    missing=set(a.cases)-{row['name'] for row in manifest}
    if missing:raise ValueError('Unknown cases: '+', '.join(sorted(missing)))
(REPORT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
def run(cmd,path):
    path.with_suffix('.command.txt').write_text(shlex.join(map(str,cmd))+'\n')
    with path.open('w') as log:subprocess.run(list(map(str,cmd)),stdout=log,stderr=log,check=True)
def phase(row):
    chrom=row['chromosome'];case=OUT/row['name'];case.mkdir(exist_ok=True)
    for arm in ('clean','recovery2'):
        d=case/arm;d.mkdir(exist_ok=True)
        if (d/'done').exists():continue
        sites=f'{chrom}.sites.striped.vcf.gz' if chrom=='chr20' else f'{chrom}.sites.vcf.gz'
        cmd=[binary,'collect-hybrid-variation','--ref',ROOT/f'test_data/chm13v2.0.{chrom}.renamed.fa','--bam',ROOT/f'test_data/HG002_{chrom}_hifi_mapped_to_CHM13_{chrom}_annotated.bam','--graph-sites',ROOT/'test_data'/sites,'--gaf',ROOT/f'test_data/HG002.{chrom}.annotated.coord.gaf.gz','-r',row['region'],'--chunk-size','500000','--threads','2','--link-by-alleles','--block-link-window','8','--min-read-margin','2','-o',d/'candidates.tsv','--phased-vcf-out',d/'native.vcf','-b',d/'phased.bam']
        if arm=='recovery2':cmd+=['--recover-gaps','--gap-recovery-report',d/'tiers.tsv']
        run(cmd,d/'phase.log');run(['samtools','index',d/'phased.bam'],d/'index.log')
        # Restrict to the same caller records in the solve window.
        region=row['region'].replace('CHM13#0#','')
        run(['python3',ROOT/'scripts/phase_vcf_from_hp.py',d/'phased.bam',DATA/f'shared_calls/{chrom}/deepvariant.vcf.gz',d/'shared.vcf','--region',region,'--support-cache',d/'support.tsv','--rebuild-support-cache'],d/'projection.log')
        (d/'done').touch()
    print('phased '+row['name'],flush=True)
with ThreadPoolExecutor(max_workers=a.workers) as pool:list(pool.map(phase,manifest))
for chrom in sorted({r['chromosome'] for r in manifest}):
    names=set()
    for row in manifest:
        if row['chromosome']!=chrom:continue
        for arm in ('clean','recovery2'):
            with pysam.AlignmentFile(str(OUT/row['name']/arm/'phased.bam')) as bam:names.update(r.query_name for r in bam)
    namefile=OUT/f'{chrom}.names.txt';namefile.write_text('\n'.join(sorted(names))+'\n')
    truth=OUT/f'{chrom}.truth.bam'
    with truth.open('wb') as out:subprocess.run(['samtools','view','-b','-N',str(namefile),str(DATA/f'truth/{chrom}/diplinator_merged.bam')],stdout=out,check=True)
    def evaluate(job):
        row,arm=job;d=OUT/row['name']/arm;e=d/'read_eval';e.mkdir(exist_ok=True)
        if (e/'summary.json').exists():return
        run(['python3',ROOT/'scripts/evaluate_phase_accuracy.py',d/'phased.bam',truth,'0','0','5','',e,'samtools','','','',''],d/'eval.log')
    jobs=[(row,arm) for row in manifest if row['chromosome']==chrom for arm in ('clean','recovery2')]
    with ThreadPoolExecutor(max_workers=a.workers) as pool:list(pool.map(evaluate,jobs))
print('finished '+str(len(manifest))+' cases',flush=True)
