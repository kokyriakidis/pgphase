from pathlib import Path
import json
rows=json.loads(Path('evaluations/2026-10-05-source-quality-path/hiphase-variants.json').read_text());source={}
for l in Path('test_data/tmp_gap_fix60/source-chain/13/matrix.recovery.chunk0.window0.chunk-1.recovery-source.tsv').read_text().splitlines():
 f=l.split('\t')
 if f[0]=='VAR':source[(int(f[2]),f[3],int(f[6]),f[7])]=f
for row in rows:
 pos=row['pos'];ref,alt=row['alleles']
 while ref and alt and ref[-1]==alt[-1]:ref,alt=ref[:-1],alt[:-1]
 while ref and alt and ref[0]==alt[0]:pos+=1;ref,alt=ref[1:],alt[1:]
 kind='X' if len(ref)==len(alt)==1 else 'D' if not alt else 'I' if not ref else 'M'
 key=(pos,kind,len(ref),alt);f=source[key]
 row['normalized_source_key']=key;row['source_category_bitmask']=int(f[4]);row['source_ps']=int(f[8]);row['source_gt']=[int(f[9]),int(f[10])]
 assert row['source_gt']==row['gt'] and row['source_ps']==row['ps']
Path('evaluations/2026-10-05-source-quality-path/hiphase-variants.json').write_text(json.dumps(rows,indent=2)+'\n')
