import zipfile,re
import pandas as pd
def read_xlsx(path,sheet_name=None,header=0):
 import xml.etree.ElementTree as ET
 ns={'x':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
 with zipfile.ZipFile(path) as zz:
  shared=[]
  if 'xl/sharedStrings.xml' in zz.namelist():
   shared=[''.join(t.text or '' for t in si.iter('{'+ns['x']+'}t')) for si in ET.fromstring(zz.read('xl/sharedStrings.xml'))]
  wb=ET.fromstring(zz.read('xl/workbook.xml'));sheets=wb.find('x:sheets',ns);sheet=next((s for s in sheets if s.attrib['name']==sheet_name),sheets[0]) if sheet_name else sheets[0];rid=sheet.attrib['{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id']
  rels=ET.fromstring(zz.read('xl/_rels/workbook.xml.rels'));target=next(x.attrib['Target'] for x in rels if x.attrib['Id']==rid);target=target.lstrip('/') if target.startswith('/') else 'xl/'+target
  rows=[]
  for row in ET.fromstring(zz.read(target)).findall('.//x:sheetData/x:row',ns):
   values={}
   for c in row:
    label=re.match(r'[A-Z]+',c.attrib['r']).group();j=0
    for cc in label:j=j*26+ord(cc)-64
    v=c.find('x:v',ns);typ=c.attrib.get('t');val=v.text if v is not None else ''
    if typ=='s':val=shared[int(val)]
    if typ=='inlineStr':val=''.join(t.text or '' for t in c.findall('.//x:t',ns))
    if typ not in ['s','inlineStr','str'] and val!='':
     try:val=float(val);val=int(val) if val.is_integer() else val
     except (ValueError,AttributeError):pass
    values[j-1]=None if val=='' else val
   rows.append(values)
  head=rows[header];width=max(head)+1;names=[head.get(k,'Unnamed_'+str(k)) for k in range(width)]
  return pd.DataFrame([[row.get(k,None) for k in range(width)] for row in rows[header+1:]],columns=names)
