from pathlib import Path
import requests,json,hashlib
from bs4 import BeautifulSoup
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'fsqca_preparation_20261016';EV=OUT/'audit_evidence';RAW=OUT/'institutional_sources';TOOLS=OUT/'tools'
for p in [OUT,EV,RAW,TOOLS]:p.mkdir(parents=True,exist_ok=True)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def fetch(name,url):
 p=RAW/name
 if not p.exists():
  r=requests.get(url,timeout=150);p.write_bytes(r.content);(RAW/(name+'.receipt.json')).write_text(json.dumps(dict(url=url,status=r.status_code,bytes=len(r.content),sha256=sha(p)),indent=2));r.raise_for_status()
 return p
if __name__=='__main__':
 # Snapshot all existing research/result artifacts touched by this scoped review.
 protected=[]
 for folder in [ROOT/'econometric_results',ROOT/'reports/final_model_selection_20261016']:
  for p in folder.iterdir():
   if p.is_file():protected.append(dict(path=str(p),sha256=sha(p)))
 (EV/'protected_existing_files.json').write_text(json.dumps(protected,indent=2))
 page=fetch('R_windows_base.html','https://cran.r-project.org/bin/windows/base/');soup=BeautifulSoup(page.read_bytes(),'html.parser');link=next(a['href'] for a in soup.find_all('a',href=True) if a['href'].endswith('-win.exe'));url=requests.compat.urljoin('https://cran.r-project.org/bin/windows/base/',link);p=fetch('R_official_installer.exe',url);fetch('R_MD5','https://cran.r-project.org/bin/windows/base/md5sum.R-4.6.1.txt');print('Installer',p,'bytes',p.stat().st_size,flush=True)
 fetch('R_installer_script.iss','https://raw.githubusercontent.com/wch/r-source/trunk/src/gnuwin32/installer/R.iss')
