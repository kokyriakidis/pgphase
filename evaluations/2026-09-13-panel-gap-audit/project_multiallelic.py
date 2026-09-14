#!/usr/bin/env python3
"""Extend the frozen biallelic projection cache with genotype-aware multi-allelic calls."""
import argparse,csv,json,os,shlex,subprocess
from pathlib import Path
import pysam
ROOT=Path(__file__).resolve().parents[2]
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--chromosome',required=True)
parser.add_argument('--case',type=Path,required=True)
parser.add_argument('--projector',type=Path,default=ROOT/'scripts/phase_vcf_from_hp.py')
a=parser.parse_args();d=a.case;chrom=a.chromosome
DATA=Path.home()/'Downloads/pgphase-eval-data'
source=DATA/f'shared_calls/{chrom}/deepvariant.vcf.gz'
with pysam.VariantFile(str(source)) as vin:
    with pysam.VariantFile(str(d/'multi.input.vcf'),'w',header=vin.header) as vout:
        for r in vin:
            if len(r.alts or ())>1:vout.write(r)
pysam.tabix_compress(str(d/'multi.input.vcf'),str(d/'multi.input.vcf.gz'),force=True)
pysam.tabix_index(str(d/'multi.input.vcf.gz'),preset='vcf',force=True)
env=dict(os.environ,PYTHONPATH=str(ROOT/'scripts'))
def run(cmd,name):
    (d/(name+'.command.txt')).write_text(shlex.join(map(str,cmd))+'\n')
    with (d/(name+'.log')).open('w') as f:subprocess.run(list(map(str,cmd)),stdout=f,stderr=f,env=env,check=True)
run(['python3',a.projector,d/'phased.bam',d/'multi.input.vcf.gz',d/'multi.phased.vcf','--region',chrom,'--support-cache',d/'multi.support.tsv','--rebuild-support-cache'],'multi_projection')
with (d/'combined.support.tsv').open('w') as out:
    out.write((d/'shared.support.tsv').read_text())
    with (d/'multi.support.tsv').open() as f:
        next(f);out.writelines(f)
run(['python3',a.projector,d/'phased.bam',source,d/'shared.multiallelic.vcf','--region',chrom,'--support-cache',d/'combined.support.tsv'],'combined_projection')
pysam.tabix_compress(str(d/'shared.multiallelic.vcf'),str(d/'shared.multiallelic.vcf.gz'),force=True)
pysam.tabix_index(str(d/'shared.multiallelic.vcf.gz'),preset='vcf',force=True)
whatshap=Path.home()/'micromamba/envs/bench-phasers/bin/whatshap'
for label,vcf in [('original',d/'shared.vcf.gz'),('multiallelic',d/'shared.multiallelic.vcf.gz')]:
    run([whatshap,'compare','--tsv-pairwise',d/(label+'.compare.tsv'),'--switch-error-bed',d/(label+'.switch_errors.bed'),'--names','truth,'+label,DATA/f'results/chr12-18-20-comparison/{chrom}/truth.vcf.gz',vcf],label+'.compare')
    run([whatshap,'stats','--tsv',d/(label+'.stats.tsv'),'--block-list',d/(label+'.blocks.tsv'),vcf],label+'.stats')
    length=int((ROOT/f'test_data/chm13v2.0.{chrom}.renamed.fa.fai').read_text().splitlines()[0].split()[1])
    run(['python3',ROOT/'scripts/compute_ngc50.py',d/(label+'.blocks.tsv'),d/(label+'.switch_errors.bed'),d/(label+'.ngc50.json'),'--genome-size',length],label+'.ngc50')
print((d/'combined_projection.log').read_text())
print((d/'original.compare.tsv').read_text())
print((d/'multiallelic.compare.tsv').read_text())
