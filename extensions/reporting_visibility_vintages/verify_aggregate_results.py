#!/usr/bin/env python3
"""Public aggregate consistency checks; does NOT re-estimate source-dependent models."""
from pathlib import Path
import csv, json, math, sys
ROOT=Path(__file__).resolve().parent
def tab(name):
    with (ROOT/'results'/name).open(encoding='utf-8-sig',newline='') as f:
        return list(csv.DictReader(f))
def check(cond,msg):
    if not cond: raise AssertionError(msg)
def close(a,b): return math.isclose(float(a),b,abs_tol=1e-9,rel_tol=1e-10)
r=json.loads((ROOT/'results'/'analysis_receipt.json').read_text())
c=tab('visibility_coverage.csv')
check(len(c)==5 and sum(int(x['grid']) for x in c)==2170,'Visibility grid must be 2,170')
check(sum(int(x['observed']) for x in c)==1399,'Recorded published amounts must be 1,399')
check(sum(int(x['missing']) for x in c)==771,'Missing source records must be 771')
v=tab('vintage_coverage.csv')
check(sum(int(x['n']) for x in v)==1953,'Matched source-vintage grid must be 1,953')
for key,target in [('both',1121),('2025 only',142),('2026 only',13),('neither',677)]:
    check(sum(int(x[key]) for x in v)==target,'Vintage state mismatch '+key)
hi=next(x for x in v if x['group']=='High income')
check(int(hi['2025 only'])==136,'High-income 2025-only mismatch')
l=tab('visibility_logit.csv')
rq=next(x for x in l if x['term']=='rq01')
check(close(rq['p_cluster'],0.016856185194254674),'RQ logit inference mismatch')
b=next(x for x in l if x['term']=='bank10')
check(close(b['p_cluster'],0.7775805294778066),'Bank logit inference mismatch')
w=tab('visibility_within_lpm.csv')
rqw=next(x for x in w if x['term']=='rq01')
check(close(rqw['p'],0.2074809755414807),'LPM RQ inference mismatch')
matched=tab('matched_vintage_FE.csv')
check(len(matched)==4 and all(int(x['n'])==863 and int(x['clusters'])==116 for x in matched),'Matched FE sample mismatch')
check(close(next(x for x in matched if x['vintage']=='y25' and x['term']=='rq')['estimate'],2.9204223133846368),'2025-vintage matched FE RQ mismatch')
check(close(next(x for x in matched if x['vintage']=='y26' and x['term']=='rq')['estimate'],2.5323016751103142),'2026-vintage matched FE RQ mismatch')
check(r['visibility']['N']==1522 and r['vintage']['positive_both']==1099,'Receipt mismatch')
# Verify immutable byte identities of all published result tables and receipt.
import hashlib
hash_rows=list(csv.DictReader((ROOT/'RESULTS_SHA256.csv').open(encoding='utf-8',newline='')))
check(len(hash_rows)==9,'Expected nine published frozen result artifacts')
for item in hash_rows:
    data=(ROOT/item['path']).read_bytes()
    check(len(data)==int(item['bytes']),'Published result byte-length mismatch: '+item['path'])
    check(hashlib.sha256(data).hexdigest()==item['sha256'],'Published result SHA256 mismatch: '+item['path'])

print(json.dumps({'status':'PASS','scope':'aggregate identities only, not a raw-data re-estimation','published_result_files_checked':8,'n_visibility':1522,'n_matched_vintage':863},indent=2))
