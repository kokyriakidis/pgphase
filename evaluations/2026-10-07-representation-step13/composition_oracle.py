from collections import Counter

def placements(pos,ref,alt):
 if ref==alt:return [None]
 left=0
 while left<min(len(ref),len(alt)) and ref[left]==alt[left]:left+=1
 right=0
 while right<min(len(ref),len(alt)) and ref[-1-right]==alt[-1-right]:right+=1
 trimmed=min(len(ref),len(alt),left+right)
 return [(pos+p,len(ref)-trimmed,alt[p:len(alt)-(trimmed-p) if trimmed-p else len(alt)])
         for p in range(max(0,trimmed-right),min(left,trimmed)+1)]

def render(reference,beg,events):
 cursor=beg;parts=[];insertions=set()
 for pos,length,alt in sorted(set(e for e in events if e is not None)):
  if pos<cursor or not length and pos in insertions:return None
  parts.extend([reference[cursor-beg:pos-beg],alt]);cursor=pos+length
  if not length:insertions.add(pos)
 return ''.join(parts)+reference[cursor-beg:]

class Oracle:
 def __init__(self,beg,reference):self.beg=beg;self.reference=reference;self.cache={}
 def prepare(self,raw):
  if raw in self.cache:return self.cache[raw]
  pos,ref,alt=raw;ref=ref.upper();alt=alt.upper()
  if pos<self.beg or pos+len(ref)>self.beg+len(self.reference) or len(ref)>4096 or len(alt)>4096 or set(ref+alt)-set('ACGT') or self.reference[pos-self.beg:pos-self.beg+len(ref)]!=ref:
   result=('unsupported',None,None,None,None)
  else:
   choices=placements(pos,ref,alt);key=choices[0]
   if key is None:result=('noop',None,None,None,choices)
   else:
    effect=render(self.reference,self.beg,[key]);region=None
    if choices[0][0]!=choices[-1][0]:region=(choices[0][0],choices[-1][0]+key[1])
    result=('edit',key,region,effect,choices)
  self.cache[raw]=result;return result
 def compose(self,raw_edits):
  canonical=[(pos,ref.upper(),alt.upper()) for pos,ref,alt in raw_edits]
  data=[self.prepare(raw) for raw in dict.fromkeys(canonical)]
  if any(r[0]=='unsupported' for r in data):return 'unsupported',None
  data=[r for r in data if r[0]=='edit'];counts=Counter(r[1] for r in data)
  for _,key,region,_,_ in data:
   if region:
    if counts[key]>1:return 'overlap',None
    if any(k!=key and k[0]<=region[1] and k[0]+k[1]>=region[0] for k in counts):return 'overlap',None
  effects={}
  for _,key,_,effect,_ in data:
   if effect in effects and effects[effect]!=key:return 'overlap',None
   effects[effect]=key
  sequence=render(self.reference,self.beg,list(counts))
  if sequence is None:return 'overlap',None
  if len(sequence)>4096:return 'unsupported',None
  return 'valid',sequence
