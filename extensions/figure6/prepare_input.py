from pathlib import Path
import hashlib,urllib.request
p=Path(__file__).resolve().parent/'input'/'GSE316734_Raw_count_file_rsemgenes.txt.gz';p.parent.mkdir(exist_ok=True)
url='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE316nnn/GSE316734/suppl/GSE316734_Raw_count_file_rsemgenes.txt.gz';sha='392af9c884b59df16e9036654125a63c73878ee171d81f3e7c63dc799a392aae'
if not p.is_file():
 with urllib.request.urlopen(url,timeout=120) as r: data=r.read()
 if hashlib.sha256(data).hexdigest()!=sha:raise ValueError('Downloaded GSE316734 file SHA-256 mismatch; stop.')
 p.write_bytes(data)
actual=hashlib.sha256(p.read_bytes()).hexdigest()
if actual!=sha:raise ValueError('Existing GSE316734 file SHA-256 mismatch; stop.')
print('verified',p.name,actual)
