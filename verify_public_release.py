"""Read-only, standard-library checks of public aggregate outputs; no model rerun."""
from pathlib import Path
import csv, hashlib, json, math
root=Path(__file__).resolve().parent
with (root/'MANIFEST.csv').open() as f:
    files=list(csv.DictReader(f))
for row in files:
    p=root/row['path']
    assert p.is_file(),row['path']
    assert p.stat().st_size==int(row['bytes']),row['path']
    assert hashlib.sha256(p.read_bytes()).hexdigest()==row['sha256'],row['path']
def rows(path):
    with (root/path).open(encoding='utf-8-sig') as f:return list(csv.DictReader(f))
co=rows('evidence/econometric_results/audit_evidence/coefficients.csv')
base={r['term']:float(r['estimate']) for r in co if r['model']=='Primary Baseline'}
assert math.isclose(base['FD'],-0.010955379185357828,abs_tol=1e-10)
assert math.isclose(base['Governance'],2.3028463099040692,abs_tol=1e-10)
meta=rows('evidence/econometric_results/audit_evidence/model_metadata.csv')
b=next(r for r in meta if r['model']=='Primary Baseline')
assert int(b['n'])==1104 and int(b['country_clusters'])==137
sol=rows('evidence/fsqca_results_20261016/audit_evidence/Python_solutions.csv')
assert all(not r['expression'] for r in sol if r['scenario']=='primary_Ypc')
alt=[r for r in sol if r['scenario']=='primary_Yabs']
assert alt and all(r['expression']=='FD + RQ' for r in alt)
for r in alt:
    assert math.isclose(float(r['consistency']),0.8433721379673498,abs_tol=1e-12)
    assert math.isclose(float(r['coverage']),0.558149900216072,abs_tol=1e-12)
print(json.dumps({'status':'PASS','manifest_files_checked':len(files),'checks':'File integrity and selected aggregate output identities only','models_reestimated':False},indent=2))
