#!/usr/bin/env python3
"""Evaluate completed chromosome runs against the same frozen variant truth."""
import argparse,csv,json,shlex,shutil,subprocess
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--chromosome',default='chr20');p.add_argument('--arms',nargs='+',default=['clean','recovery2','graph_support2'])
a=p.parse_args();chrom=a.chromosome
ROOT=Path(__file__).resolve().parents[2];REPORT=Path(__file__).resolve().parent
OUT=Path('/tmp/pgphase-panel-gap-audit')/chrom
FROZEN=Path.home()/f'Downloads/pgphase-eval-data/results/chr12-18-20-comparison/{chrom}'
WHATSHAP=Path.home()/'micromamba/envs/bench-phasers/bin/whatshap'
length=int((ROOT/f'test_data/chm13v2.0.{chrom}.renamed.fa.fai').read_text().splitlines()[0].split()[1])
def run(cmd,path):
    path.with_suffix('.command.txt').write_text(shlex.join(map(str,cmd))+'\n')
    with path.open('w') as log:subprocess.run(list(map(str,cmd)),stdout=log,stderr=log,check=True)
rows=[]
for arm in a.arms:
    d=OUT/arm
    if not (d/'done').exists():raise RuntimeError(f'Run not complete: {d}')
    run([WHATSHAP,'compare','--tsv-pairwise',d/'compare.tsv','--switch-error-bed',d/'switch_errors.bed','--names','truth,'+arm,FROZEN/'truth.vcf.gz',d/'shared.vcf.gz'],d/'compare.log')
    run([WHATSHAP,'stats','--tsv',d/'stats.tsv','--block-list',d/'blocks.tsv',d/'shared.vcf.gz'],d/'stats.log')
    run(['python3',ROOT/'scripts/compute_ngc50.py',d/'blocks.tsv',d/'switch_errors.bed',d/'ngc50.json','--genome-size',length],d/'ngc50.log')
    for f in ('compare.tsv','stats.tsv','ngc50.json','read_eval/summary.json','tiers.tsv'):
        if (d/f).exists():shutil.copy2(d/f,REPORT/(chrom+'.'+arm+'.'+Path(f).name))
for name in a.arms+['hiphase','whatshap','whatshap_opt','longphase']:
    if name in a.arms:
        d=OUT/name;compare=d/'compare.tsv';ngc=d/'ngc50.json';stats=d/'stats.tsv';read=d/'read_eval/summary.json'
    else:
        compare=FROZEN/(name+'.compare.tsv');ngc=FROZEN/(name+'.ngc50.json');stats=FROZEN/(name+'.stats.tsv');read=FROZEN/f'eval/{name}_reads/summary.json'
    c=next(csv.DictReader(compare.open(),delimiter='\t'));n=json.loads(ngc.read_text());r=json.loads(read.read_text())
    rows.append({'method':name,'truth_assessed_variants':c['covered_variants'],'variant_switches':c['all_switches'],'variant_hamming':c['blockwise_hamming'],'variant_hamming_rate':c['blockwise_hamming_rate'],'ngc50_bp':n['ngc50_bp'],'read_truth_evaluated':r['total_reads_evaluated'],'discordant_reads':r['discordant_reads']})
with (REPORT/(chrom+'.metrics.tsv')).open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
print((REPORT/(chrom+'.metrics.tsv')).read_text())
