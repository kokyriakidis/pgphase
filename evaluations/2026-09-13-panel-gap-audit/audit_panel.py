#!/usr/bin/env python3
"""Inventory exact competitor-phased sites and block links missing from current pgphase."""
import argparse,bisect,csv,json,sys
from collections import Counter,defaultdict
from pathlib import Path
import pysam
csv.field_size_limit(sys.maxsize)
ROOT=Path(__file__).resolve().parents[2]
REPORT=Path(__file__).resolve().parent
DATA=Path.home()/'Downloads/pgphase-eval-data'
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--chromosome',required=True)
parser.add_argument('--arm',default='recovery1')
parser.add_argument('--out-root',type=Path,default=Path('/tmp/pgphase-panel-gap-audit'))
parser.add_argument('--multiallelic',action='store_true',help='Audit the genotype-aware projection and combined cache')
a=parser.parse_args();chrom=a.chromosome;case=a.out_root/chrom/a.arm
label=a.arm+('.multiallelic' if a.multiallelic else '')
sys.path.insert(0,str(ROOT/'scripts'))
from phase_vcf_from_hp import alignment_allele

def key(r):return (r.pos,r.ref,tuple(r.alts or ()))
def records(path):
    with pysam.VariantFile(str(path)) as f:
        for r in f:
            s=next(iter(r.samples.values()));gt=s.get('GT')
            if gt and len(gt)==2 and None not in gt and gt[0]!=gt[1]:
                yield key(r),s.get('PS') if s.phased else None
pg=dict(records(case/('shared.multiallelic.vcf.gz' if a.multiallelic else 'shared.vcf.gz')))
native=dict(records(case/'native.vcf'))
candidates=defaultdict(list)
with (case/'candidates.tsv').open() as f:
    for r in csv.DictReader(f,delimiter='\t'):candidates[int(r['POS'])].append(r['CATEGORY'])
positions=sorted(candidates)
tiers=[]
if (case/'tiers.tsv').exists():
    with (case/'tiers.tsv').open() as f:tiers=list(csv.DictReader(f,delimiter='\t'))
recovery_gaps=defaultdict(list)
for r in tiers:recovery_gaps[int(r['GAP_LEFT']),int(r['GAP_RIGHT'])].append(r)
support=defaultdict(lambda: defaultdict(lambda: defaultdict(Counter)))
with (case/('combined.support.tsv' if a.multiallelic else 'shared.support.tsv')).open() as f:
    for r in csv.DictReader(f,delimiter='\t'):
        support[int(r['POS0'])+1,r['REF'],tuple(r['ALT'].split(','))][int(r['PS'])][int(r['HP'])].update({0:int(r['REF_COUNT']),1:int(r['ALT_COUNT'])})
raw_bam=pysam.AlignmentFile(str(ROOT/f'test_data/HG002_{chrom}_hifi_mapped_to_CHM13_{chrom}_annotated.bam'))
contig=next(c for c in raw_bam.references if c.split('#')[-1]==chrom)
def phased_details(path):
    details={}
    with pysam.VariantFile(str(path)) as f:
        for r in f:
            sample=next(iter(r.samples.values()));gt=sample.get('GT')
            if sample.phased and gt and len(gt)==2 and None not in gt and gt[0]!=gt[1]:
                # This assembly-based truth has chromosome-wide paternal/maternal GT
                # and no PS field; explicit PS records retain their block boundaries.
                ps=sample.get('PS') if 'PS' in f.header.formats else r.contig
                details[key(r)]=(tuple(gt),ps)
    return details
truth=phased_details(DATA/f'results/chr12-18-20-comparison/{chrom}/truth.vcf.gz')
def pair_truth(left,right,competitor):
    if left not in truth or right not in truth:return 'endpoints_not_both_in_truth'
    tl,tr=truth[left],truth[right];cl,cr=competitor[left],competitor[right]
    if tl[1] is None or tl[1]!=tr[1]:return 'different_truth_blocks'
    if sorted(tl[0])!=sorted(cl[0]) or sorted(tr[0])!=sorted(cr[0]):return 'genotype_mismatch'
    return 'concordant' if (tl[0][0]!=cl[0][0])==(tr[0][0]!=cr[0][0]) else 'discordant'
def support_reason(k):
    if len(k[2])!=1 and not a.multiallelic:return 'projection_multiallelic_unsupported'
    if not support[k]:return 'no_tagged_allele_support'
    counts=max(support[k].values(),key=lambda hp:sum(sum(c.values()) for c in hp.values()))
    n1,n2=sum(counts[1].values()),sum(counts[2].values())
    if min(n1,n2)<2:return 'fewer_than_two_tagged_reads_on_one_haplotype'
    a1,c1=counts[1].most_common(1)[0];a2,c2=counts[2].most_common(1)[0]
    if c1/n1<0.7 or c2/n2<0.7:return 'mixed_alleles_within_haplotype'
    if a1==a2:return 'both_haplotypes_favor_same_allele'
    return 'other_projection_failure'
cache={};rows=[]
for tool in ('hiphase','whatshap','whatshap_opt','longphase'):
    frozen=DATA/f'results/chr12-18-20-comparison/{chrom}'
    accuracy={}
    with (frozen/f'eval/{tool}_reads/per_phase_set.tsv').open() as f:
        for r in csv.DictReader(f,delimiter='\t'):accuracy[int(r['phase_set'])]=float(r['accuracy'])
    competitor=phased_details(frozen/tool/'phased.vcf.gz')
    blocks=defaultdict(list)
    for k,ps in records(frozen/tool/'phased.vcf.gz'):
        if ps is not None:blocks[ps].append(k)
    for ps,sites in blocks.items():
        sites.sort()
        for i,k in enumerate(sites):
            if pg.get(k) is None:
                reason=support_reason(k)
                rows.append(dict(chromosome=chrom,tool=tool,competitor_ps=ps,competitor_accuracy=accuracy.get(ps),kind='unphased_site',competitor_pair_truth='not_applicable',left=k[0],right=k[0],left_pg_ps='',right_pg_ps='',reason=reason,primary_spanning='',spanning_mapq30='',native_left_ps='',native_right_ps='',snp_pair_callable='',snp_pair_dominant='',recovery_status='',candidate_categories=';'.join(candidates[k[0]]),ref=k[1],alt=','.join(k[2])))
        present=[k for k in sites if pg.get(k) is not None]
        for left,right in zip(present,present[1:]):
            if pg[left]==pg[right]:continue
            pair=left,right
            if pair not in cache:
                names=set();qualified=set();pairs={}
                for r in raw_bam.fetch(contig,left[0]-1,right[0]):
                    if r.is_secondary or r.is_supplementary or r.is_unmapped or r.is_duplicate or r.is_qcfail:continue
                    if r.reference_start>left[0]-1 or r.reference_end<right[0]:continue
                    names.add(r.query_name)
                    if r.mapping_quality < 30: continue
                    qualified.add(r.query_name)
                    if len(left[1])==len(right[1])==1 and len(left[2])==len(right[2])==1 and len(left[2][0])==len(right[2][0])==1:
                        la=alignment_allele(r,left[0]-1,left[1],left[2][0]);ra=alignment_allele(r,right[0]-1,right[1],right[2][0])
                        if la is not None and ra is not None:pairs[r.query_name]=la^ra
                counts=Counter(pairs.values());cache[pair]=(len(names),len(qualified),len(pairs),max(counts.values(),default=0))
            span,qualified,callable_reads,dominant=cache[pair]
            overlapping=[v[-1]['STATUS'] for (x,y),v in recovery_gaps.items() if x<=right[0] and y>=left[0]]
            labels=Counter(label for pos in positions[bisect.bisect_left(positions,left[0]):bisect.bisect_right(positions,right[0])] for label in candidates[pos])
            reason='no_spanning_read_between_projected_endpoints' if span==0 else 'all_spanning_reads_below_mapq30' if qualified==0 else 'singleton_mapq30_between_projected_endpoints' if qualified==1 else 'multiple_mapq30_spanning_reads_need_link_audit'
            if native.get(left) is not None and native.get(left)==native.get(right):reason='read_tag_projection_break_inside_native_block'
            rows.append(dict(chromosome=chrom,tool=tool,competitor_ps=ps,competitor_accuracy=accuracy.get(ps),kind='split_block',competitor_pair_truth=pair_truth(left,right,competitor),left=left[0],right=right[0],left_pg_ps=pg[left],right_pg_ps=pg[right],reason=reason,primary_spanning=span,spanning_mapq30=qualified,native_left_ps=native.get(left),native_right_ps=native.get(right),snp_pair_callable=callable_reads,snp_pair_dominant=dominant,recovery_status=';'.join(overlapping) or 'not_attempted',candidate_categories=json.dumps(labels,sort_keys=True),ref='',alt=''))
raw_bam.close()
output=REPORT/f'{chrom}.{label}.missed.tsv'
with output.open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]) if rows else ['chromosome','kind'],delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
summary=Counter((r['kind'],r['reason']) for r in rows)
(REPORT/f'{chrom}.{label}.summary.json').write_text(json.dumps({'rows':len(rows),'unique_sites':len({r['left'] for r in rows if r['kind']=='unphased_site'}),'unique_links':len({(r['left'],r['right']) for r in rows if r['kind']=='split_block'}),'unique_block_pairs':len({tuple(sorted((r['left_pg_ps'],r['right_pg_ps']))) for r in rows if r['kind']=='split_block'}),'counts':{'/'.join(k):v for k,v in summary.items()}},indent=2)+'\n')
print((REPORT/f'{chrom}.{label}.summary.json').read_text())
