"""Authorized exploratory second stage; archived original engine is untouched."""
from stage2_setup import *
import sys,importlib.util,warnings,platform
import pandas as pd,numpy as np
from scipy import stats
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from threadpoolctl import threadpool_limits
spec=importlib.util.spec_from_file_location('stage2_engine',SCRIPTS/'engine_frozen_copy.py');eng=importlib.util.module_from_spec(spec);spec.loader.exec_module(eng)
SOURCE=ROOT/'second_stage_investigation_20261016/processed_data/second_stage_candidate_grid.csv'
d=pd.read_csv(SOURCE,low_memory=False);d['group']=d.income_group
bench=eng.core(d[d.sample_B_Bank_credit.eq(True)].copy());assert len(bench)==1104 and bench.country_code.nunique()==137
raw=ROOT/'second_stage_investigation_20261016/raw_data'
rise=pd.read_csv(raw/'WB_RISE.csv',low_memory=False);pillars=['WB_RISE_RE_GOV','WB_RISE_RE_ELEC','WB_RISE_RE_HE_CO','WB_RISE_RE_TRNS','WB_RISE_RE_LVL_PLYNG_FLD'];r=rise[rise.INDICATOR.isin(pillars+['WB_RISE_RE_ALL'])].pivot(index=['REF_AREA','TIME_PERIOD'],columns='INDICATOR',values='OBS_VALUE');diff=r[pillars].mean(axis=1)-r.WB_RISE_RE_ALL
gate={'source':'https://rise.esmap.org/methodology','official_methodology_verified':'When indicators change scores are back-calculated under new methodology; series starts 2010; renewable pillar revised from seven to five in 2024; legal-policy desk research and partner databases','downloaded_data_edition':'Data360 WB_RISE, dataset 2025-04-01 / series 2025-04-07 metadata, annual2010–2023','rows_reconciled':len(r),'all_five_components_complete':bool(r[pillars].notna().all().all()),'max_difference_total_vs_equal_mean_five':float(diff.abs().max()),'within_1e_8':bool(diff.abs().le(1e-8).all()),'metadata_discrepancy':'Legacy metadata uses three pillars/25–75 traffic lights, current four pillars/33–67; no traffic-light calibration used','caution':'Retrospective back-calculation is publisher construction, not independent annual surveys or real-time observability'}
(EV/'RISE_definition_validation.json').write_text(json.dumps(gate,indent=2),encoding='utf-8');assert gate['within_1e_8'],gate
basecols=['FD','Governance','log_real_GDP_pc','log_population','trade_GDP','inflation']
target={'A':'RISE','B':'FDI','C':'Deficit'}
FITS={};COEF=[];META=[];MEM=[];SUP=[];BOOT=[];JOINT=[];LOCO=[];PPML=[];FAIL=[];COMPARE=[];COV=[];CHECK=[]
save=lambda rows,n:pd.DataFrame(rows).to_csv(EV/(n+'.csv'),index=False)
def support(family,timing='main'):
 cols={'A':['L1_RISE_RE'],'B':['L1_FDI_GDP','gdp_growth_pct'],'C':['L1_access_deficit']}[family]
 if timing=='all_lag':cols += ['L1_FD.AST.PRVT.GD.ZS','L1_regulatory_quality_estimate','L1_gdp_per_capita_constant_2015_usd','L1_population','L1_trade_pct_gdp','L1_inflation_cpi_pct'];cols=[('L1_gdp_growth_pct' if x=='gdp_growth_pct' else x) for x in cols]
 z=bench[np.isfinite(bench[cols].astype(float)).all(axis=1)].copy()
 if timing=='all_lag':z=z[z.L1_population.gt(0)&z.L1_gdp_per_capita_constant_2015_usd.gt(0)]
 return eng.core(z)
def design(z,family,baseline=False,interaction=False,asinh_infl=False,asinh_fdi=False,all_lag=False,gov=None,centers=None):
 pre='L1_' if all_lag else '';fd=z[pre+'FD.AST.PRVT.GD.ZS'];g=z[gov or pre+'regulatory_quality_estimate'];centers=centers or {'FD':float(fd.mean()),'RQ':float(g.mean()),'RISE':float(z.L1_RISE_RE.mean())}
 X=pd.DataFrame(dict(FD=fd-centers['FD'],Governance=g-centers['RQ'],log_real_GDP_pc=np.log(z[pre+'gdp_per_capita_constant_2015_usd']),log_population=np.log(z[pre+'population']),trade_GDP=z[pre+'trade_pct_gdp'],inflation=z[pre+'inflation_cpi_pct'])).reset_index(drop=True)
 if asinh_infl:X['inflation']=np.arcsinh(X.inflation)
 if not baseline:
  if family=='A':X['RISE']=z.L1_RISE_RE.to_numpy()-centers['RISE']
  elif family=='B':X['FDI']=np.arcsinh(z.L1_FDI_GDP.to_numpy()) if asinh_fdi else z.L1_FDI_GDP.to_numpy();X['growth']=z['L1_gdp_growth_pct' if all_lag else 'gdp_growth_pct'].to_numpy()
  elif family=='C':X['Deficit']=z.L1_access_deficit.to_numpy()
 if interaction:X['RQ_x_RISE']=X.Governance*X.RISE
 return X,centers
def fit(z,family,**kwargs):
 z=eng.core(z);X,centers=design(z,family,**{k:v for k,v in kwargs.items() if k!='inference'});eng.variables=lambda *args,**kw:(X,(centers['FD'],centers['RQ']));f=eng.Fit(z,inference=kwargs.get('inference',True));f.design_centers=centers;f.design_options={k:v for k,v in kwargs.items() if k not in ['inference','centers']};f.family=family;return f
def record(label,z,family,role,**kwargs):
 f=fit(z,family,**kwargs);FITS[label]=f
 COEF.extend(dict(model=label,family=family,role=role,**x) for x in f.coefficients());META.append(dict(model=label,family=family,role=role,n=f.N,clusters=f.G,rank=f.rank,residual_df=f.N-f.rank,first_year=int(f.z.year.min()),last_year=int(f.z.year.max()),within_R2=f.within_R2,overall_R2=f.full_R2,SSE=f.sse,max_within_VIF=float(f.VIF.max()),condition_index=float(f.condition),centers=json.dumps(f.design_centers),options=json.dumps(f.design_options)))
 MEM.extend(dict(model=label,**r) for r in f.z[['country_code','year','group']].to_dict('records'))
 for g in ['Pooled']+eng.G:
  zz=f.z if g=='Pooled' else f.z[f.z.group.eq(g)];SUP.append(dict(model=label,income_group=g,n=len(zz),clusters=zz.country_code.nunique(),years=zz.year.nunique()))
 for i,a in enumerate(f.names):
  for j,b in enumerate(f.names):COV.append(dict(model=label,term1=a,term2=b,covariance_CR2=f.cov[i,j]))
 X=np.column_stack([f.D,f.V/f.norm]);ref=sm.OLS(f.y,X).fit();beta=ref.params[-f.p:]/f.norm;assert np.allclose(beta,f.beta,rtol=1e-7,atol=1e-8);assert np.allclose(ref.resid,f.e,atol=1e-7)
 CHECK.append(dict(model=label,max_coefficient_discrepancy=float(np.max(np.abs(beta-f.beta))),max_residual_discrepancy=float(np.max(np.abs(ref.resid-f.e))),CR2_identity_error=f.working_model_covariance_error))
 print('OLS',label,f.N,f.G,flush=True);return f
def compare(ext,base):
 a,b=FITS[ext],FITS[base];assert set(zip(a.z.country_code,a.z.year))==set(zip(b.z.country_code,b.z.year))
 for n in basecols:
  c=a.contrast(a.unit(n),n);r=b.contrast(b.unit(n),n);COMPARE.append(dict(extension=ext,matched_baseline=base,term=n,extension_estimate=c['estimate'],baseline_estimate=r['estimate'],coefficient_change=c['estimate']-r['estimate'],extension_CI95_low=c['CI95_low'],extension_CI95_high=c['CI95_high'],baseline_CI95_low=r['CI95_low'],baseline_CI95_high=r['CI95_high'],extension_within_R2=a.within_R2,baseline_within_R2=b.within_R2,within_R2_change=a.within_R2-b.within_R2,partial_R2=(b.sse-a.sse)/b.sse,n=a.N,clusters=a.G,interpretation='Same rows and nuisance definition; change is not a causal mediation estimate'))
with threadpool_limits(limits=1):
 for family,expected in [('A',(914,109)),('B',(1094,135)),('C',(1094,136))]:
  z=support(family);assert (len(z),z.country_code.nunique())==expected
  record(family+' matched baseline',z,family,'Matched baseline',baseline=True);f=record(family+' primary',z,family,'Primary exploratory');compare(family+' primary',family+' matched baseline')
  JOINT.append(dict(model=family+' primary',**f.joint([f.unit(target[family])],target[family]+' zero')))
 record('A readiness interaction',support('A'),'A','Secondary interaction',interaction=True);f=FITS['A readiness interaction'];JOINT.append(dict(model='A readiness interaction',**f.joint([f.unit('RQ_x_RISE')],'RQ × lagged RISE zero')))
 # Secondary transforms and timing are frozen before these fits.
 for family in ['A','B','C']:
  z=support(family)
  record(family+' asinh inflation baseline',z,family,'Sensitivity',baseline=True,asinh_infl=True);record(family+' asinh inflation',z,family,'Sensitivity',asinh_infl=True);compare(family+' asinh inflation',family+' asinh inflation baseline')
  z=support(family,'all_lag');record(family+' all-X lag baseline',z,family,'Timing sensitivity',baseline=True,all_lag=True);record(family+' all-X lag',z,family,'Timing sensitivity',all_lag=True);compare(family+' all-X lag',family+' all-X lag baseline')
 record('B asinh signed FDI',support('B'),'B','FDI-scale sensitivity',asinh_fdi=True);compare('B asinh signed FDI','B matched baseline')
 for gov,name in [('VA_EST','Voice and accountability'),('PV_EST','Political stability')]:record(name+' benchmark',bench,'A','Institutional-definition sensitivity',baseline=True,gov=gov)
 # Exactly three primary tests; bootstrap same three plus the labelled secondary interaction.
 for label,term in [('A primary','RISE'),('B primary','FDI'),('C primary','Deficit'),('A readiness interaction','RQ_x_RISE')]:
  f=FITS[label];BOOT.append(dict(model=label,**f.bootstrap([f.unit(term)],term+' zero',draws=9999,seed=20261016)));print('BOOT complete',label,flush=True)
  save(BOOT,'wild_bootstrap')
 for label in ['A primary','B primary','C primary','A readiness interaction']:
  f=FITS[label]
  for code in f.codes:
   ff=fit(f.z[f.z.country_code.ne(code)],f.family,**f.design_options,centers=f.design_centers,inference=False)
   for n in f.names:
    LOCO.append(dict(model=label,omitted_country=code,term=n,estimate=ff.beta[ff.names.index(n)],full_sample_estimate=f.beta[f.names.index(n)],change=ff.beta[ff.names.index(n)]-f.beta[f.names.index(n)],n=ff.N,clusters=ff.G,centering='Fixed parent centers; no observations trimmed from parent'))
  print('LOCO complete',label,flush=True)
 # PPML: all three extensions plus their exact-support baselines; nonlinear CR1 only.
 for label in ['A matched baseline','A primary','B matched baseline','B primary','C matched baseline','C primary']:
  try:
   parent=FITS[label];z=eng.core(parent.z,ppml=True);f=fit(z,parent.family,**parent.design_options,centers=parent.design_centers,inference=False);X=np.column_stack([f.D,f.V/f.norm]);pos=z[eng.Y].gt(0).to_numpy();_,rpos,_=eng.basis(X[pos]);_,rank,_=eng.basis(X)
   if rpos!=rank:raise RuntimeError('Positive-outcome design rank deficient; general separation not cleared; PPML withheld')
   with warnings.catch_warnings(record=True) as ww:res=sm.GLM(z[eng.Y].to_numpy(),X,family=sm.families.Poisson()).fit(maxiter=200,tol=1e-10,cov_type='cluster',cov_kwds={'groups':z.country_code.to_numpy(),'use_correction':True})
   if not res.converged:raise RuntimeError('IRLS failed to converge')
   score=float(np.max(np.abs(X.T@(z[eng.Y].to_numpy()-res.mu)))/max(1,z[eng.Y].sum()));assert score<1e-6
   beta=res.params[-f.p:]/f.norm;V=res.cov_params()[-f.p:,-f.p:]/f.norm[:,None]/f.norm[None,:]
   for j,n in enumerate(f.names):
    se=float(np.sqrt(V[j,j]));p=2*stats.t.sf(abs(beta[j]/se),f.G-1);crit=stats.t.ppf(.975,f.G-1);PPML.append(dict(model=label,term=n,estimate=beta[j],SE_country_CR1=se,p_value=p,reference_df=f.G-1,CI95_low=beta[j]-crit*se,CI95_high=beta[j]+crit*se,n=f.N,clusters=f.G,rank=rank,positive_rank=rpos,converged=True,iterations=res.fit_history['iteration'],relative_score=score,deviance=res.deviance,inference='Country CR1/t(G−1); nonlinear CR2 and wild bootstrap not implemented',warnings='; '.join(str(w.message) for w in ww)))
   print('PPML complete',label,flush=True)
  except Exception as e:FAIL.append(dict(model=label,procedure='PPML',reason=repr(e)));print('PPML failed',label,repr(e),flush=True)
coef=pd.DataFrame(COEF);mask=np.array([r.model in ['A primary','B primary','C primary'] and r.term==target[r.family] for r in coef.itertuples()]);coef['p_Holm_primary3']=np.nan;coef.loc[mask,'p_Holm_primary3']=multipletests(coef.loc[mask,'p_value'],method='holm')[1]
bt=pd.DataFrame(BOOT);mask=bt.model.isin(['A primary','B primary','C primary']);bt['p_wild_Holm_primary3']=np.nan;bt.loc[mask,'p_wild_Holm_primary3']=multipletests(bt.loc[mask,'p_wild_restricted'],method='holm')[1]
for name,rows in [('coefficients',coef),('model_metadata',META),('sample_membership',MEM),('group_support',SUP),('wild_bootstrap',bt),('joint_tests',JOINT),('LOCO_coefficients',LOCO),('PPML_results',PPML),('matched_comparisons',COMPARE),('CR2_covariances',COV),('implementation_checks',CHECK)]:save(rows,name)
save(FAIL if FAIL else [dict(model='All authorized procedures',procedure='OLS/CR2/bootstrap/PPML',reason='No execution failures')],'execution_status')
old=pd.read_csv(ROOT/'econometric_results/audit_evidence/coefficients.csv');old[old.model.isin(['Government effectiveness Baseline','Rule of law Baseline','Control of corruption Baseline'])].to_csv(EV/'reused_original_governance_results.csv',index=False)
extreme=bench[['country_code','country_name','year','group','FDI_GDP','L1_FDI_GDP']].copy();extreme['abs_lagged_FDI']=extreme.L1_FDI_GDP.abs();extreme.nlargest(20,'abs_lagged_FDI').to_csv(EV/'largest_FDI_observations.csv',index=False)
(EV/'run_receipt.json').write_text(json.dumps(dict(status='Completed',OLS_models=len(META),wild_bootstrap_tests=len(BOOT),PPML_models=len(set(r['model'] for r in PPML)),failures=FAIL,source_hash=sha(SOURCE),protocol_hash=sha(EV/'estimation_protocol_frozen.json'),python=platform.python_version(),numpy=np.__version__,fsQCA_run=False),indent=2),encoding='utf-8')
assert all(sha(Path(p))==s for p,s in json.loads((EV/'protected_before.json').read_text()).items())
print(coef.loc[mask if len(mask)==len(coef) else coef.role.eq('Primary exploratory')].head().to_string(index=False),flush=True)
print('COMPLETED',len(META),'OLS',len(BOOT),'bootstraps',len(PPML),'PPML coefficient rows',len(FAIL),'failures',flush=True)
