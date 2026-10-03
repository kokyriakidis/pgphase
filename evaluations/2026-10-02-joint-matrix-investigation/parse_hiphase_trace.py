"""Recover genuine local HiPhase calls from serial trace output.

Match each complete allele vector to BAM fetch order, with the exact HiPhase
flag/MAPQ filter. Assert record and variant counts before accepting the mapping.
"""
import re,json
from pathlib import Path
import pysam
root=Path('test_data/tmp_gap_next29')
for c,w in [(7,'7264321-7280346'),(24,'24121713-24131707')]:
 blocks=[];variants=[]
 for line in (root/f'hiphase{c}.trace.log').read_text().splitlines():
  if 'Solving problem: PhaseBlock' in line:
   m=re.search(r'coordinates: "(.+):(\d+)-(\d+)".*num_variants: (\d+)',line)
   block={'chrom':m[1],'beg':int(m[2]),'end':int(m[3]),'nvars':int(m[4]),'rows':[]};blocks.append(block)
  elif 'TRACE hiphase::read_parsing] Variant {' in line:
   m=re.search(r'variant_type: (\w+), position: (\d+), ref_len: (\d+), prefix_len: (\d+), postfix_len: (\d+), allele0: \[(.*?)\], allele1: \[(.*?)\].*index_allele0: (\d+), index_allele1: (\d+)',line)
   assert m,line
   p,q=int(m[4]),int(m[5]);a=[]
   for j in [6,7]:
    s=bytes(map(int,m[j].split(', '))).decode();a.append(s[p:len(s)-q if q else len(s)])
   variants.append({'type':m[1],'pos':int(m[2]),'alleles':a,'vcf_indices':[int(m[8]),int(m[9])]})
  elif 'All alleles [' in line:
   alleles=line.split('All alleles [')[1].split(']')[0].split(', ')
   assert len(variants)==block['nvars']==len(alleles)
   if 'variants' not in block:block['variants']=variants
   else:assert variants==block['variants']
   block['rows'].append(alleles);variants=[]
 with pysam.AlignmentFile(f'test_data/tmp_gap_next28/{w}/input.bam') as bam:
  for b in blocks:
   reads=[r for r in bam.fetch(b['chrom'],b['beg'],b['end']+1) if not(r.flag & (4|256|512|1024)) and r.mapping_quality>=5]
   assert len(reads)==len(b['rows']),(c,b['beg'],len(reads),len(b['rows']))
   b['reads']=[{'name':r.query_name,'mapq':r.mapping_quality,'flag':r.flag,'beg':r.reference_start,'end':r.reference_end,'alleles':a} for r,a in zip(reads,b.pop('rows'))]
 (root/f'hiphase{c}.matrix.json').write_text(json.dumps(blocks,indent=2)+'\n')
 print(c,[(b['beg'],len(b['variants']),len(b['reads'])) for b in blocks])
