from pathlib import Path
from collections import Counter,defaultdict
import importlib.util,pysam,json
spec=importlib.util.spec_from_file_location('a','evaluations/2026-10-06-next-source-component-block/read_audit.py');a=importlib.util.module_from_spec(spec);spec.loader.exec_module(a)
truth={f[0]:f[1]=='PATERNAL' for l in open('test_data/derived/chr20_truth_hap.tsv') if len(f:=l.rstrip().split('\t'))==2}
tags,status=a.assignments(Path('test_data/tmp_gap_fix90/baseline/4/phased.bam'),truth)
fa=pysam.FastaFile('test_data/chm13v2.0.chr20.renamed.fa');chrom='CHM13#0#chr20'
def dist(a,b):
 row=list(range(len(b)+1))
 for i,x in enumerate(a,1):
  nxt=[i]
  for j,y in enumerate(b,1):nxt.append(min(nxt[-1]+1,row[j]+1,row[j-1]+(x!=y)))
  row=nxt
 return row[-1]
def snp(read,pos):
 for qi,rp in read.get_aligned_pairs(matches_only=True):
  if rp==pos-1:return read.query_sequence[qi],read.query_qualities[qi]
 return None,None
reports=[]
with pysam.AlignmentFile('test_data/HG002_chr20_hifi_mapped_to_CHM13_chr20_annotated.bam') as bam:
 for pos,lengths in [(4625183,[4,8]),(4637299,[2,3])]:
  longest=max(lengths);motif=fa.fetch(chrom,pos-1,pos-1+longest)
  period=next(i for i in range(1,len(motif)+1) if len(motif)%i==0 and all(b==motif[k%i] for k,b in enumerate(motif)))
  left=pos;right=pos+longest
  while fa.fetch(chrom,left-2,left-1)==fa.fetch(chrom,left-2+period,left-1+period):left-=1
  while fa.fetch(chrom,right-1,right)==fa.fetch(chrom,right-1-period,right-period):right+=1
  beg,end=left-16,right+16
  ref=fa.fetch(chrom,beg-1,end).upper();expected=[ref[:pos-beg]+ref[pos-beg+n:] for n in lengths]
  counts=Counter();rows=[]
  for r in bam.fetch(chrom,pos-1,pos+longest):
   if r.is_secondary or r.is_supplementary or r.is_duplicate or r.is_qcfail or r.is_unmapped or r.mapping_quality<20 or r.mapping_quality==255:continue
   pairs={rp+1:qi for qi,rp in r.get_aligned_pairs(matches_only=True)}
   qb,qe=pairs.get(beg),pairs.get(end)
   if qb is None or qe is None:continue
   qual=min(r.query_qualities[qb],r.query_qualities[qe]);query=r.query_sequence[qb:qe+1]
   ds=[dist(query,e) for e in expected];observed=len(ref)-len(query)
   diffs=[abs(observed-l) for l in lengths]
   allele=0 if diffs[0]<diffs[1] else 1 if diffs[1]<diffs[0] else -1
   physical=allele>=0 and diffs[allele]<=2 and ds[allele]<ds[1-allele] and ds[allele]<=diffs[allele]+1
   if not physical or qual<10:continue
   row={'qname':r.query_name,'parent':'PAT' if truth.get(r.query_name) else 'MAT','allele':allele,'qual':qual,'flank_qualities':[r.query_qualities[qb],r.query_qualities[qe]],'mapq':r.mapping_quality,'observed':observed,'distance':ds,'before_status':status.get(r.query_name),'before_tag':tags.get(r.query_name),'left_snp':snp(r,4618685),'right_snp':snp(r,4642430),'span':[r.reference_start+1,r.reference_end]}
   rows.append(row);counts[(row['parent'],allele,row['before_status'])]+=1
  print(pos,'period',period,'repeat',left,right,'counts',dict(counts))
  print('unphased',[(r['qname'],r['parent'],r['allele'],r['left_snp'],r['right_snp']) for r in rows if r['before_status']=='unphased'])
  reports.append({'pos':pos,'lengths':lengths,'period':period,'repeat':[left,right],'rows':rows})
(Path(__file__).resolve().parent/'physical-calls.json').write_text(json.dumps(reports,indent=2)+'\n')
