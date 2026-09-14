#!/usr/bin/env python3
"""Compare exact chr12 target/flank boundaries, not regional block counts."""
import csv,json
from pathlib import Path
import pysam
ROOT=Path(__file__).resolve().parents[2]
REPORT=Path(__file__).resolve().parent
DATA=Path.home()/'Downloads/pgphase-eval-data/results/chr12-18-20-comparison/chr12'
positions=(46725702,46838918,46842224,46874196)
rows=[]
for tool in ('hiphase','whatshap','whatshap_opt','longphase'):
    with pysam.VariantFile(str(DATA/tool/'phased.vcf.gz')) as f:
        observed={}
        for r in f.fetch('chr12',positions[0]-1,positions[-1]):
            if r.pos not in positions:continue
            s=next(iter(r.samples.values()))
            if s.phased and len(set(s['GT']))>1:observed[r.pos]=s.get('PS')
        rows.extend([tool,p,observed.get(p)] for p in positions)
with (Path('/tmp/pgphase-singleton-link-validation/chr12_46725702/auto.tsv')).open() as f:
    observed={int(r['POS']):int(r['PHASE_SET']) for r in csv.DictReader(f,delimiter='\t')
              if r['TYPE']=='SNP' and int(r['POS']) in positions}
rows.extend(['pgphase_link1',p,observed.get(p)] for p in positions)
with (REPORT/'endpoint_phase_sets.tsv').open('w') as f:
    w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['tool','position','phase_set']);w.writerows(rows)
with pysam.AlignmentFile(str(ROOT/'test_data/HG002_chr12_hifi_mapped_to_CHM13_chr12_annotated.bam')) as f:
    names={r.query_name for r in f.fetch('chr12',46842223,46874196)
           if not r.is_secondary and not r.is_supplementary
           and r.reference_start<=46842223 and r.reference_end>=46874196}
(REPORT/'spanning_reads.json').write_text(json.dumps({'left_snp_1based':46842224,'right_snp_1based':46874196,'primary_spanning_reads':sorted(names)},indent=2)+'\n')
print((REPORT/'endpoint_phase_sets.tsv').read_text())
print('Primary spanning reads:',len(names))
