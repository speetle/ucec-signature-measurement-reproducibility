from pathlib import Path
import json,hashlib,urllib.request,re,concurrent.futures
R=Path(__file__).resolve().parents[1]
def fetch(item):
 dest=R/item['file'];dest.parent.mkdir(parents=True,exist_ok=True)
 if dest.exists() and hashlib.sha256(dest.read_bytes()).hexdigest()==item['sha256']:return
 url=item['url'];req=urllib.request.Request(url,headers={'User-Agent':'UCEC-reproduction/1.0'})
 with urllib.request.urlopen(req,timeout=60) as response:data=response.read()
 if len(data)!=item['bytes'] or hashlib.sha256(data).hexdigest()!=item['sha256']:
  # Some source servers return cached partial responses. Reassemble exact requested ranges.
  ranges=[(x,min(item['bytes']-1,x+524287)) for x in range(0,item['bytes'],524288)]
  def part(ab):
   lo,hi=ab
   with urllib.request.urlopen(urllib.request.Request(url,headers={'Range':f'bytes={lo}-{hi}','User-Agent':'UCEC-reproduction/1.0'}),timeout=90) as response:q=response.read()
   if len(q)==item['bytes']:return q[lo:hi+1]
   if len(q)!=hi-lo+1:raise ValueError('Range length mismatch; source must be re-audited')
   return q
  with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:data=b''.join(pool.map(part,ranges))
 if hashlib.sha256(data).hexdigest()!=item['sha256']:raise ValueError('Source checksum mismatch: '+item['file'])
 tmp=dest.with_suffix(dest.suffix+'.part');tmp.write_bytes(data);tmp.replace(dest);print('Fetched and verified',item['file'])
if __name__=='__main__':
 for x in json.loads((R/'optional_sources.json').read_text()):fetch(x)
