from pathlib import Path
import concurrent.futures
import json,os,subprocess,time
import sys
root=Path('test_data/tmp_test_cache') / ('native-' + sys.argv[1]);root.mkdir(parents=True,exist_ok=True)
names=Path('test_data/tmp_test_cache/test-names.txt').read_text().splitlines()
windows=[n for n in names if n in ('chr20 gap windows','chr20 gap windows: panel totals')]
assert len(windows)==2
other=[n for n in names if n not in windows]
shards=[windows]+[other[i::3] for i in range(3)]
assert len(set(n for shard in shards for n in shard))==len(names)==86
(root/'cases.json').write_text(json.dumps(shards,indent=2)+'\n')
def run(i,shard):
 selector=','.join('"'+name.replace('\\','\\\\').replace('"','\\"')+'"' for name in shard)
 env=dict(os.environ,LD_LIBRARY_PATH='/usr/lib/x86_64-linux-gnu',PGPHASE_BIN='./pgphase',PGPHASE_TEST_CACHE='test_data/tmp_test_cache/real-cache',PGPHASE_TEST_WORKDIR=str(root/f'shard{i}'))
 t=time.monotonic()
 with (root/f'shard{i}.log').open('w') as log:
  r=subprocess.run(['./test_gap_windows',selector],env=env,stdout=log,stderr=log)
 report={'shard':i,'cases':len(shard),'returncode':r.returncode,'seconds':time.monotonic()-t}
 (root/f'shard{i}.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps(report),flush=True)
 return r.returncode
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
 codes=list(pool.map(lambda pair:run(*pair),enumerate(shards)))
raise SystemExit(any(codes))
