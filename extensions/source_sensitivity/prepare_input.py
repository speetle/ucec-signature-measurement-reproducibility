from pathlib import Path
import hashlib,urllib.request
p=Path(__file__).resolve().parent/'input'/'GSE235046_Count_table.txt.gz';p.parent.mkdir(exist_ok=True)
url='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE235nnn/GSE235046/suppl/GSE235046_Count_table.txt.gz';sha='683ec07b3d81ea8cb6b834b08e071ba9ec2b89f49c61ac303029864b1a110906'
if not p.exists():
 with urllib.request.urlopen(url,timeout=120) as r: data=r.read()
 if hashlib.sha256(data).hexdigest()!=sha:raise ValueError('Downloaded source matrix hash mismatch')
 p.write_bytes(data)
if hashlib.sha256(p.read_bytes()).hexdigest()!=sha:raise ValueError('Existing source matrix hash mismatch')
print('verified',p.name,sha)
