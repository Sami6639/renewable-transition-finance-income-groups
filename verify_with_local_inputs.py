from pathlib import Path
import hashlib,json,csv,re
import pandas as pd,numpy as np
from threadpoolctl import threadpool_limits
import argparse
parser=argparse.ArgumentParser(description='Read-only partial checks of the selected replication handoff; not a full rerun.')
parser.add_argument('--package',type=Path,required=True,help='Directory containing MANIFEST.csv and evidence/')
parser.add_argument('--output',type=Path,default=Path('verification_results.json'))
args=parser.parse_args();p=args.package.resolve();r={}
manifest=list(csv.DictReader((p/'MANIFEST.csv').open()))
r['manifest']={'entries':len(manifest),'missing':[],'mismatch':[]}
for row in manifest:
 f=p/row['path']
 if not f.exists():r['manifest']['missing'].append(row['path'])
 elif hashlib.sha256(f.read_bytes()).hexdigest()!=row['sha256'] or f.stat().st_size!=int(row['bytes']):r['manifest']['mismatch'].append(row['path'])
d=pd.read_csv(p/'evidence/reports/final_model_selection_20261016/four_income_group_final_dataset.csv')
r['panel']={'rows':len(d),'countries':d.country_code.nunique(),'years':[int(d.year.min()),int(d.year.max())],'duplicate_country_year':int(d.duplicated(['country_code','year']).sum()),'available_finance':int(d.irena2025_pdf.notna().sum()),'missing_finance':int(d.irena2025_pdf.isna().sum())}
cols=['FD.AST.PRVT.GD.ZS','regulatory_quality_estimate','gdp_per_capita_constant_2015_usd','population','trade_pct_gdp','inflation_cpi_pct','irena2025_pdf']
z=d.loc[np.isfinite(d[cols]).all(axis=1)&d.group.isin(['Low income','Lower middle income','Upper middle income','High income'])&(d.population>0)&(d.gdp_per_capita_constant_2015_usd>0)&(d.irena2025_pdf>=0)].copy()
while True:
 n=len(z);z=z.loc[(z.groupby('country_code').year.transform('size')>=2)&(z.groupby('year').country_code.transform('size')>=2)]
 if len(z)==n:break
r['primary']={'rows':len(z),'countries':z.country_code.nunique(),'zeros':int((z.irena2025_pdf==0).sum()),'group_rows':z.groupby('group').size().to_dict(),'sample_flag_agrees':set(zip(z.country_code,z.year))==set(zip(d.loc[d.sample_B_Bank_credit,'country_code'],d.loc[d.sample_B_Bank_credit,'year']))}
X=np.column_stack([pd.get_dummies(z.country_code,dtype=float),pd.get_dummies(z.year,drop_first=True,dtype=float),z[cols[:2]].to_numpy()-z[cols[:2]].mean().to_numpy(),np.log(z.gdp_per_capita_constant_2015_usd),np.log(z.population),z.trade_pct_gdp,z.inflation_cpi_pct])
with threadpool_limits(limits=1):beta=np.linalg.lstsq(X,np.log1p(z.irena2025_pdf),rcond=None)[0][-6:]
c=pd.read_csv(p/'evidence/econometric_results/audit_evidence/coefficients.csv');c=c.loc[c.model=='Primary Baseline'].set_index('term')
terms=['FD','Governance','log_real_GDP_pc','log_population','trade_GDP','inflation'];r['independent_baseline_OLS']={'coefficients':dict(zip(terms,beta.tolist())),'max_difference_vs_stored':float(np.max(np.abs(beta-c.loc[terms,'estimate'].to_numpy()))),'scope':'Direct full country/year dummy least-squares. No CR2, bootstrap or PPML rerun.'}
q=pd.read_csv(p/'evidence/fsqca_results_20261016/fsqca_calibrated_dataset.csv');r['qca']={'cases':len(q),'unique_countries':q.country_code.nunique()}
plan=json.loads((p/'evidence/fsqca_results_20261016/audit_evidence/approved_execution_plan_frozen.json').read_text())
for key,(col,(lo,ce,hi)) in plan['sets'].items():
 raw=q[col].to_numpy();v=1/(1+np.exp(-np.log(19)*(raw-ce)/np.where(raw<ce,ce-lo,hi-ce)))
 r['qca'][key+'_calibration_max_difference']=float(np.max(np.abs(v-q[key])))
x=np.maximum(q.FD,q.RQ);y=q.Yabs;xy=np.minimum(x,y);contra=np.minimum(np.minimum(x,y),1-y)
r['qca']['alternative_FD_or_RQ']={'consistency':float(xy.sum()/x.sum()),'coverage':float(xy.sum()/y.sum()),'PRI':float((xy.sum()-contra.sum())/(x.sum()-contra.sum()))}
r['absolute_path_files']=[str(f.relative_to(p)) for f in p.rglob('*') if f.is_file() and f.suffix in ['.py','.md','.R','.json','.txt','.csv','.mjs'] and re.search(r'C:[/\\]+Users',f.read_text(errors='replace'))]
args.output.write_text(json.dumps(r,indent=2));print(json.dumps(r,indent=2))
