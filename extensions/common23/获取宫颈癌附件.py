from pathlib import Path
import argparse,hashlib,urllib.request,os
url='https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-024-53830-0/MediaObjects/41467_2024_53830_MOESM6_ESM.xlsx'
expected='d57e4ecc948710bebe7293f43d5ab40a434e7ea35d94ebc70806ad0d200d22f9';size=172804197
p=argparse.ArgumentParser();p.add_argument('destination',type=Path);a=p.parse_args();dest=a.destination.expanduser()
if dest.exists():raise SystemExit('目标已存在，不覆盖；请使用新文件名')
dest.parent.mkdir(parents=True,exist_ok=True);temp=dest.with_name(dest.name+'.partial')
if temp.exists():raise SystemExit('临时文件已存在，请检查后处理')
h=hashlib.sha256();n=0
try:
 with urllib.request.urlopen(url,timeout=60) as source,temp.open('xb') as f:
  while True:
   chunk=source.read(1024*1024)
   if not chunk:break
   n+=len(chunk);h.update(chunk);f.write(chunk)
 if n!=size or h.hexdigest()!=expected:raise ValueError(f'校验失败: bytes={n} SHA256={h.hexdigest()}')
 os.replace(temp,dest);print('校验通过',dest,h.hexdigest())
except Exception:
 temp.unlink(missing_ok=True);raise
