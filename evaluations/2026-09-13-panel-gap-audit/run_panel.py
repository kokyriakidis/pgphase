#!/usr/bin/env python3
"""Run the current code over the full available three-chromosome panel."""
import argparse,hashlib,json,shlex,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
REPORT=Path(__file__).resolve().parent
OUT=Path('/tmp/pgphase-panel-gap-audit')
DATA=Path.home()/'Downloads/pgphase-eval-data'
OUT.mkdir(exist_ok=True)
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--chromosomes',nargs='+',default=['chr20'])
parser.add_argument('--arms',nargs='+',choices=['clean','recovery1','recovery2','baseline_recovery2','graph2','graph_support2','nested_seed','owner_clean','observed_owner'],default=['clean','recovery2'])
a=parser.parse_args()
# Freeze the executable so later fixes cannot change a running comparison.
import shutil
binary=OUT/'pgphase.baseline'
if not binary.exists():shutil.copy2(ROOT/'pgphase',binary)
(REPORT/'baseline.sha256').write_text(hashlib.sha256(binary.read_bytes()).hexdigest()+'\n')
current=OUT/'pgphase.context'
if not current.exists():shutil.copy2(ROOT/'pgphase',current)
(REPORT/'context.sha256').write_text(hashlib.sha256(current.read_bytes()).hexdigest()+'\n')
projector=OUT/'phase_vcf_from_hp.baseline.py'
if not projector.exists():shutil.copy2(ROOT/'scripts/phase_vcf_from_hp.py',projector)
def run(cmd,path):
    path.with_suffix('.command.txt').write_text(shlex.join(map(str,cmd))+'\n')
    with path.open('w') as log:subprocess.run(list(map(str,cmd)),stdout=log,stderr=log,check=True)
def run_chromosome(chrom):
    contig='CHM13#0#chr20' if chrom=='chr20' else chrom
    sites=f'{chrom}.sites.striped.vcf.gz' if chrom=='chr20' else f'{chrom}.sites.vcf.gz'
    for arm in a.arms:
        case=OUT/chrom/arm;case.mkdir(parents=True,exist_ok=True)
        if (case/'done').exists():continue
        print(chrom,arm,flush=True)
        solver=OUT/'pgphase.fast' if arm=='recovery1' and (OUT/'pgphase.fast').exists() else binary
        if arm=='recovery2':solver=current
        if arm=='graph2':
            solver=OUT/'pgphase.graph'
            if not solver.exists():shutil.copy2(ROOT/'pgphase',solver)
            (REPORT/'graph.sha256').write_text(hashlib.sha256(solver.read_bytes()).hexdigest()+'\n')
        if arm=='graph_support2':
            solver=OUT/'pgphase.graph_support2'
            if not solver.exists():shutil.copy2(ROOT/'pgphase',solver)
            (REPORT/'graph_support2.sha256').write_text(hashlib.sha256(solver.read_bytes()).hexdigest()+'\n')
        if arm=='nested_seed':
            solver=OUT/'pgphase.nested_seed'
            if not solver.exists():shutil.copy2(ROOT/'pgphase',solver)
            (REPORT/'nested_seed.sha256').write_text(hashlib.sha256(solver.read_bytes()).hexdigest()+'\n')
        if arm in ('owner_clean','observed_owner'):
            solver=OUT/'pgphase.observed_owner'
            if not solver.exists():shutil.copy2(ROOT/'pgphase',solver)
            (REPORT/'observed_owner.sha256').write_text(hashlib.sha256(solver.read_bytes()).hexdigest()+'\n')
        cmd=[solver,'collect-hybrid-variation','--ref',ROOT/f'test_data/chm13v2.0.{chrom}.renamed.fa','--bam',ROOT/f'test_data/HG002_{chrom}_hifi_mapped_to_CHM13_{chrom}_annotated.bam','--graph-sites',ROOT/'test_data'/sites,'--gaf',ROOT/f'test_data/HG002.{chrom}.annotated.coord.gaf.gz','-r',contig,'--chunk-size','500000','--threads','8','--link-by-alleles','--block-link-window','8','--min-read-margin','2','-o',case/'candidates.tsv','--phased-vcf-out',case/'native.vcf','-b',case/'phased.bam']
        if arm not in ('clean','owner_clean'):cmd+=['--recover-gaps','--min-block-link-reads','1' if arm=='recovery1' else '2','--gap-recovery-report',case/'tiers.tsv']
        if not (case/'phased.done').exists():
            run(['/usr/bin/time','-v']+cmd,case/'phase.log')
            (case/'phased.done').touch()
        run(['samtools','index','-@','4',case/'phased.bam'],case/'index.log')
        e=case/'read_eval';e.mkdir(exist_ok=True)
        run(['python3',ROOT/'scripts/evaluate_phase_accuracy.py',case/'phased.bam',DATA/f'truth/{chrom}/diplinator_merged.bam','0','0','5','',e,'samtools','','','',''],case/'read_eval.log')
        run(['python3',projector,case/'phased.bam',DATA/f'shared_calls/{chrom}/deepvariant.vcf.gz',case/'shared.vcf','--region',chrom,'--min-reads','2','--min-ratio','0.70','--support-cache',case/'shared.support.tsv','--rebuild-support-cache'],case/'projection.log')
        run(['bgzip','-f',case/'shared.vcf'],case/'compress.log')
        run(['tabix','-p','vcf',case/'shared.vcf.gz'],case/'tabix.log')
        (case/'done').touch()
        print('finished',chrom,arm,flush=True)

from concurrent.futures import ThreadPoolExecutor
with ThreadPoolExecutor(max_workers=2) as pool:
    list(pool.map(run_chromosome, a.chromosomes))
