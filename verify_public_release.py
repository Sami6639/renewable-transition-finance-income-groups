"""Read-only, standard-library checks of public aggregate outputs; no model rerun."""
from pathlib import Path
import csv, hashlib, json, math
root=Path(__file__).resolve().parent
with (root/'MANIFEST.csv').open(encoding='utf-8', newline='') as f:
    files=list(csv.DictReader(f))
# Git checkouts can translate LF and CRLF. Preserve the release manifest's
# original byte counts and SHA-256 hashes, and permit only that translation
# for known UTF-8 text formats. Binary and unknown formats remain byte-exact.
TEXT_SUFFIXES = {'.csv', '.py', '.json', '.R', '.md', '.txt'}
TEXT_NAMES = {'.gitignore'}

def matches_manifest(data, row):
    return (len(data) == int(row['bytes']) and
            hashlib.sha256(data).hexdigest() == row['sha256'])

def verify_file(path, row):
    data = path.read_bytes()
    if matches_manifest(data, row):
        return False
    if path.suffix in TEXT_SUFFIXES or path.name in TEXT_NAMES:
        try:
            data.decode('utf-8')
        except UnicodeDecodeError:
            pass
        else:
            lf = data.replace(b'\r\n', b'\n')
            # Do not strip whitespace, alter a BOM, or add a final newline.
            if (matches_manifest(lf, row) or
                    matches_manifest(lf.replace(b'\n', b'\r\n'), row)):
                return True
    raise AssertionError('File integrity mismatch: ' + row['path'])

line_ending_matches = 0
for row in files:
    p = root / row['path']
    if not p.is_file():
        raise AssertionError('Missing manifest file: ' + row['path'])
    line_ending_matches += verify_file(p, row)

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
print(json.dumps({'status':'PASS','manifest_files_checked':len(files),'line_ending_normalized_matches':line_ending_matches,'checks':'File integrity and selected aggregate output identities only','models_reestimated':False},indent=2))
