from planb_setup import *
import numpy as np,pandas as pd
from scipy.linalg import qr
from threadpoolctl import threadpool_limits
_blas_limits=threadpool_limits(limits=1)
d=pd.read_csv(SOURCE);assert len(d)==2170 and not d.duplicated(['country_code','year']).any()
G=['Low income','Lower middle income','Upper middle income','High income'];X=['FD','RQ','lnGDP','lnPop','Trade','Inflation']
d['FD']=d['FD.AST.PRVT.GD.ZS'];d['RQ']=d.regulatory_quality_estimate;d['lnGDP']=np.log(d.gdp_per_capita_constant_2015_usd.where(d.gdp_per_capita_constant_2015_usd>0));d['lnPop']=np.log(d.population.where(d.population>0));d['Trade']=d.trade_pct_gdp;d['Inflation']=d.inflation_cpi_pct;d['Yabs']=d.irena2025_pdf;d['Ypc']=d.Yabs*1e6/d.population;d['t']=d.year
def eligible(z):return z[np.isfinite(z[X+['Yabs','Ypc']]).all(axis=1)&z.Yabs.ge(0)&z.group.isin(G)].copy()
def core(z):
 while True:
  n=len(z);z=z[z.groupby('country_code').t.transform('size').ge(2)&z.groupby('t').country_code.transform('size').ge(2)].copy()
  if len(z)==n:return z.sort_values(['country_code','t']).reset_index(drop=True)
a=core(eligible(d));assert len(a)==1104 and a.country_code.nunique()==137
assert set(zip(a.country_code,a.year))==set(zip(d.loc[d.sample_B_Bank_credit,'country_code'],d.loc[d.sample_B_Bank_credit,'year']))
models=[];support=[];audit=[];window_rows=[];bases={}
def add(name,z,method,outcome,scope):
 z=core(z);z['y']=np.log1p(z['Yabs' if outcome=='absolute' else 'Ypc']);z['model']=name;z.to_csv(EV/(name+'_data.csv'),index=False)
 models.append(dict(model=name,method=method,outcome=outcome,scope=scope,n=len(z),clusters=z.country_code.nunique(),periods=z.t.nunique(),first_period=int(z.t.min()),last_period=int(z.t.max())))
 for g in ['Pooled']+G:
  zz=z if g=='Pooled' else z[z.group==g];support.append(dict(model=name,income_group=g,observations=len(zz),countries=zz.country_code.nunique(),periods=zz.t.nunique()))
 bases[name]=z
for outcome in ['absolute','per_capita']:
 for method in ['TWFE','CRE']+(['MundlakOLS'] if outcome=='absolute' else []):add('annual_'+method+'_'+outcome,a,method,outcome,'Annual native benchmark')
bal=a[a.groupby('country_code').t.transform('size').eq(10)].copy()
for method in ['TWFE','CRE']:add('balanced_'+method,bal,method,'absolute','Complete ten-year histories')
lag=d.copy();prev=d[['country_code','year']+X].copy();prev['year']+=1;prev=prev.rename(columns={x:'lag_'+x for x in X});lag=lag.merge(prev,on=['country_code','year'],how='left',validate='one_to_one');lag[X]=lag[['lag_'+x for x in X]].values;lag=core(eligible(lag))
for method in ['TWFE','CRE']:add('lagged_'+method,lag,method,'absolute','One-year lag all X; no pre2014 values')
for horizon,windows in [(2,[(2014,2015),(2016,2017),(2018,2019),(2020,2021),(2022,2023)]),(3,[(2014,2016),(2017,2019),(2020,2022)])]:
 records=[]
 for code,z in d.groupby('country_code'):
  for start,end in windows:
   zz=z[z.year.between(start,end)].sort_values('year');fin=bool(len(zz)==horizon and zz.Yabs.notna().all() and zz.Yabs.ge(0).all());cov=bool(len(zz)==horizon and np.isfinite(zz[X+['Ypc']]).all(axis=1).all());groupsok=bool(zz.group.isin(G).all());ok=fin and cov and groupsok
   window_rows.append(dict(horizon=horizon,country_code=code,start_year=start,end_year=end,required_years=horizon,observed_financing_years=int(zz.Yabs.notna().sum()),complete_financing=fin,complete_covariates=cov,eligible_before_FE=ok,income_group=zz.iloc[0].group,income_changed=zz.group.nunique()>1))
   if ok:
    row=dict(country_code=code,country_name=zz.iloc[0].country_name,group=zz.iloc[0].group,t=start,year=start,end_year=end,horizon=horizon,Yabs=zz.Yabs.mean(),Ypc=zz.Ypc.mean(),sovereign_scope=zz.iloc[0].sovereign_scope,income_changed=zz.group.nunique()>1)
    row.update({x:zz[x].mean() for x in X});records.append(row)
 w=core(pd.DataFrame(records));w.to_csv(EV/f'complete_{horizon}year_windows.csv',index=False)
 keys={(r.country_code,y) for r in w.itertuples() for y in range(r.t,r.end_year+1)};match=a[[ (r.country_code,r.year) in keys for r in a.itertuples() ]].copy();assert len(match)==horizon*len(w)
 for outcome in ['absolute','per_capita']:
  add(f'window{horizon}_TWFE_'+outcome,w,'TWFE',outcome,f'Complete nonoverlap {horizon}-year windows')
  add(f'matched{horizon}_annual_TWFE_'+outcome,match,'TWFE',outcome,f'Annual exact country-year support of {horizon}-year windows')
  if horizon==2:add(f'window2_CRE_'+outcome,w,'CRE',outcome,'Complete two-year windows')
pd.DataFrame(models).to_csv(EV/'model_plan.csv',index=False);pd.DataFrame(support).to_csv(EV/'income_support.csv',index=False);pd.DataFrame(window_rows).to_csv(EV/'window_eligibility.csv',index=False)
# Freeze full model designs before any coefficient estimation. Drop only exact aliases.
designs={};designmeta=[]
for m in models:
 z=bases[m['model']];T=pd.get_dummies(z.t,prefix='year',drop_first=True,dtype=float)
 if m['method']=='TWFE':D=pd.concat([pd.get_dummies(z.country_code,prefix='country',dtype=float),T,z[X]],axis=1)
 else:
  B=z[X].groupby(z.country_code).transform('mean');W=z[X]-B;bt=T.groupby(z.country_code).transform('mean')
  D=pd.concat([pd.DataFrame({'Intercept':np.ones(len(z))}),W.add_prefix('W_'),B.add_prefix('B_'),T,bt.add_prefix('profile_')],axis=1)
 # Sequential QR-independent-column selection preserves scientific covariate terms first.
 keep=[];r=0
 for c in D:
  mat=D[keep+[c]].to_numpy();norm=np.linalg.norm(mat,axis=0);mat=mat/np.where(norm>0,norm,1)
  nr=np.linalg.matrix_rank(mat,tol=1e-10)
  if nr>r:keep.append(c);r=nr
 omitted=[c for c in D if c not in keep];D=D[keep];D.to_csv(EV/(m['model']+'_design.csv'),index=False);designs[m['model']]=dict(columns=keep,omitted_exact_aliases=omitted)
 designmeta.append(dict(model=m['model'],n=len(z),rank=r,residual_df=len(z)-r,condition_index=float(np.linalg.cond(D/np.linalg.norm(D,axis=0))),omitted_exact_aliases=';'.join(omitted)))
(EV/'frozen_design_columns.json').write_text(json.dumps(designs,indent=2),encoding='utf-8');pd.DataFrame(designmeta).to_csv(EV/'design_metadata.csv',index=False)
(EV/'design_freeze_receipt.json').write_text(json.dumps(dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),source_sha256=sha(SOURCE),preanalysis_sha256=sha(EV/'preanalysis_plan_frozen.json'),models=len(models),design_sha256={m['model']:sha(EV/(m['model']+'_design.csv')) for m in models}),indent=2),encoding='utf-8')
print(pd.DataFrame(models).to_string(index=False),flush=True)
