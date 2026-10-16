"""Approved exploratory estimates, frozen IRENA 2025 only. Run from project results.
See execution_instructions.md; no fsQCA, fabricated data or replacement of missings.
"""
from pathlib import Path
import json,hashlib,itertools,platform,warnings,sys
import numpy as np,pandas as pd
from scipy import stats
from scipy.optimize import linprog
from threadpoolctl import threadpool_limits
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from econometric_engine import *
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'econometric_results';EV=OUT/'audit_evidence';FIG=OUT/'econometric_figures'
for p in [OUT,EV,FIG]:p.mkdir(parents=True,exist_ok=True)
SOURCE=ROOT/'reports/final_model_selection_20261016/four_income_group_final_dataset.csv';sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();source_hash=sha(SOURCE)
d=pd.read_csv(SOURCE);assert len(d)==2170 and not d.duplicated(['country_code','year']).any() and d.year.between(2014,2023).all();assert not any('2026' in c for c in d);assert d.loc[d.year.gt(2020),'imf_official'].isna().all()
save=lambda df,n:df.to_csv(EV/(n+'.csv'),index=False)
COEF=[];META=[];CONTR=[];JOINT=[];MEM=[];COV=[];BOOT=[];FAIL=[];CHECK=[];SUPPORT=[];FITS={};COUNTRY=[];OBS=[];PPML=[]
def reportfit(label,z,fd=FD,gov=RQ,pre='',block='Baseline',outcome='log1p',groupcoding='Historical annual',cache=False,unit=1):
 pre_rows=len(z);a=core(z)
 try:
  f=Fit(a,fd,gov,pre,block,outcome,unit=unit)
 except Exception as e:FAIL.append(dict(model=label,procedure='OLS/CR2',status='FAILED',reason=repr(e)));print('FAILED',label,repr(e),flush=True);return
 tag=dict(model=label,block=block,FD_measure=fd,governance=gov,timing='Lagged X' if pre else 'Current',outcome=outcome,income_coding=groupcoding)
 COEF.extend(dict(**tag,**r) for r in f.coefficients());META.append(dict(**tag,complete_rows=pre_rows,n=f.N,country_clusters=f.G,years=f.year_count,first_year=a.year.min(),last_year=a.year.max(),model_rank=f.rank,residual_df=f.N-f.rank,within_R2=f.within_R2,overall_R2=f.full_R2,FD_center=f.centers[0],governance_center=f.centers[1],max_VIF=f.VIF.max(),condition_index=f.condition,absorbed_terms=', '.join(f.absorbed),FinData_rows_any_imputed=int(a.imputed_components.gt(0).sum()) if fd=='findata_official' else None,working_cov_identity_error=f.working_model_covariance_error))
 for g in ['Pooled']+G:
  z=a if g=='Pooled' else a[a.group.eq(g)];SUPPORT.append(dict(model=label,income_group=g,rows=len(z),clusters=z.country_code.nunique(),years=z.year.nunique()))
 for r in a[['country_code','year','group']].to_dict('records'):MEM.append(dict(model=label,**r))
 for i,n1 in enumerate(f.names):
  for j,n2 in enumerate(f.names):COV.append(dict(model=label,term1=n1,term2=n2,covariance_CR2=f.cov[i,j]))
 if 'FD_x_Governance' in f.names:
  govvalues=a[pre+gov];qs=sorted(set(np.quantile(govvalues,[0,.1,.25,.5,.75,.9,1])))
  for r in qs:
   c=f.unit('FD')+(r-f.centers[1])*f.unit('FD_x_Governance');CONTR.append(dict(**tag,evaluation='FD slope at governance',governance_value=r,observed_governance_min=govvalues.min(),observed_governance_max=govvalues.max(),**f.contrast(c,f'FD slope at {r:.6g}')))
 if 'Income' in block or block=='Combined':
  slopes={g:f.unit('FD')+(f.unit('FD_x_'+g) if g!=G[0] else 0) for g in G}
  for g,c in slopes.items():CONTR.append(dict(**tag,evaluation='Group FD slope',income_group=g,**f.contrast(c,g)))
  pairs=[]
  for a1,a2 in itertools.combinations(G,2):pairs.append(dict(**tag,evaluation='Pairwise FD slope contrast',contrast=f'{a2} minus {a1}',**f.contrast(slopes[a2]-slopes[a1],f'{a2} minus {a1}')))
  adjusted=multipletests([r['p_value'] for r in pairs],method='holm')[1]
  for r,p in zip(pairs,adjusted):r['p_Holm_6_pair_family']=p;CONTR.append(r)
  R=np.array([f.unit('FD_x_'+g) for g in G[1:]]);JOINT.append(dict(model=label,**f.joint(R,'Three income-interaction coefficients jointly zero')))
 if 'FD_x_Governance' in f.names:JOINT.append(dict(model=label,**f.joint([f.unit('FD_x_Governance')],'FD × governance interaction zero')))
 # Independent full FE-dummy OLS coefficient and residual verification.
 X=np.column_stack([f.D,f.V/f.norm]);ols=sm.OLS(f.y,X).fit();beta=ols.params[-f.p:]/f.norm;assert np.allclose(beta,f.beta,rtol=1e-6,atol=1e-8);assert np.allclose(ols.resid,f.e,atol=1e-7)
 CHECK.append(dict(model=label,full_OLS_beta_max_difference=float(np.max(np.abs(beta-f.beta))),full_OLS_residual_max_difference=float(np.max(np.abs(ols.resid-f.e))),CR2_working_model_identity_error=f.working_model_covariance_error))
 if cache:FITS[label]=f
 print(label,f.N,f.G,'FD',f.beta[0],flush=True)
 return f
with threadpool_limits(limits=1):
 # Common support first; definitions and eligible periods fixed before estimation.
 plans=[('AB matched',2023,[('Bank depth',FD),('FinData','findata_official')]),('AC matched',2020,[('FinData','findata_official'),('IMF','imf_official')]),('BC matched',2020,[('Bank depth',FD),('IMF','imf_official')]),('ABC matched',2020,[('Bank depth',FD),('FinData','findata_official'),('IMF','imf_official')])]
 for support,end,measures in plans:
  z=d[d.year.le(end)]
  for _,fd in measures:z=eligible(z,fd)
  z=core(z)
  for name,fd in measures:
   for block in ['Baseline','RQ interaction','Income interaction']:reportfit(f'{support}: {name}: {block}',z,fd,block=block)
 # Primary native support and core frozen-membership check.
 z=eligible(d);a=core(z);assert len(a)==1104 and a.country_code.nunique()==137;assert set(zip(a.country_code,a.year))==set(zip(d.loc[d.sample_B_Bank_credit,'country_code'],d.loc[d.sample_B_Bank_credit,'year']))
 for block in ['Baseline','RQ interaction','Income interaction','Combined']:reportfit('Primary '+block,z,block=block,cache=True)
 for fd,label,end in [('findata_official','FinData native',2023),('imf_official','IMF native',2020)]:
  for block in ['Baseline','RQ interaction','Income interaction']:reportfit(label+' '+block,eligible(d[d.year.le(end)],fd),fd,block=block)
 for block in ['Baseline','RQ interaction','Income interaction']:reportfit('Lagged '+block,eligible(d,pre='L1_'),pre='L1_',block=block)
 joint=core(eligible(eligible(d),pre='L1_'))
 for pre,name in [('', 'Current matched lag'),('L1_','Lag matched current')]:reportfit(name,joint,pre=pre)
 for gov,label in [('government_effectiveness_estimate','Government effectiveness'),('rule_of_law_estimate','Rule of law'),('control_of_corruption_estimate','Control of corruption')]:
  for block in ['Baseline','RQ interaction']:reportfit(label+' '+block,eligible(d,gov=gov),gov=gov,block=block)
 for restriction,zz in [('Sovereign scope',d[d.sovereign_scope.eq(True)]),('Stable income',d[d.country_code.isin(d.groupby('country_code').filter(lambda t:t.group.isin(G).all() and t.group.nunique()==1).country_code)])]:
  for block in ['Baseline','RQ interaction','Income interaction']:reportfit(restriction+' '+block,eligible(zz),block=block)
 b=d.copy();lookup=d[d.year.eq(2014)].set_index('country_code').group;b['group']=b.country_code.map(lookup);reportfit('Fixed 2014 income',eligible(b),block='Income interaction',groupcoding='Fixed at 2014')
 for block in ['Baseline','RQ interaction','Income interaction']:reportfit('Positive only '+block,eligible(d[d[Y].gt(0)]),block=block,outcome='positive_log')
 # Prespecified scale-dependence investigation, not a preferred-model search.
 for unit,label in [(1e6,'log1p in USD'),(.001,'log1p in billion USD')]:reportfit(label,z,unit=unit)
 # Primary bootstrap hypotheses fixed before observing their estimates.
 f=FITS['Primary RQ interaction'];BOOT.append(dict(model='Primary RQ interaction',**f.bootstrap([f.unit('FD_x_Governance')],'FD × RQ zero')))
 f=FITS['Primary Income interaction'];R=np.array([f.unit('FD_x_'+g) for g in G[1:]]);BOOT.append(dict(model='Primary Income interaction',**f.bootstrap(R,'Income interactions jointly zero')))
 for g in G[1:]:BOOT.append(dict(model='Primary Income interaction',**f.bootstrap([f.unit('FD_x_'+g)],g+' minus low-income slope zero')))
 # All-country deletion coefficient checks, no data trimming in primary estimates.
 for label in ['Primary Baseline','Primary RQ interaction','Primary Income interaction']:
  f=FITS[label]
  for code in f.codes:
   cut=core(f.z[f.z.country_code.ne(code)]);ff=Fit(cut,block=f.block,centers=f.centers,inference=False)
   for term in ['FD']+[n for n in f.names if n.startswith('FD_x_')]:
    if term in ff.names:COUNTRY.append(dict(model=label,omitted_country=code,term=term,estimate=ff.beta[ff.names.index(term)],full_sample_estimate=f.beta[f.names.index(term)],change=ff.beta[ff.names.index(term)]-f.beta[f.names.index(term)],rows=ff.N,clusters=ff.G,centering='Fixed full-parent centers for coefficient comparability'))
  h=np.sum(f.Q*f.Q,axis=1);sigma2=f.sse/(f.N-f.rank);cook=f.e*f.e*h/(f.rank*sigma2*np.maximum(1-h,1e-12)**2);dfb=f.B*(f.e/np.maximum(1-h,1e-12))[:,None]/np.sqrt(np.diag(f.cov))
  for i,r in f.z.iterrows():OBS.append(dict(model=label,country_code=r.country_code,country_name=r.country_name,year=r.year,income_group=r.group,financing_musd2022=r[Y],log1p_financing=f.y[i],fitted=f.fitted[i],residual=f.e[i],leverage=h[i],Cook_distance_descriptive=cook[i],max_abs_CR2_scaled_DFBETA=float(np.max(np.abs(dfb[i]))),FD_DFBETA_CR2_scaled=dfb[i,0]))
 # Real-data validation: singleton-cluster limit reduces to independently computed HC2.
 f=FITS['Primary Baseline'];X=np.column_stack([f.D,f.V/f.norm]);ols=sm.OLS(f.y,X).fit();hc2=ols.get_robustcov_results(cov_type='HC2').cov_params()[-f.p:,-f.p:]/f.norm[:,None]/f.norm[None,:];h=np.sum(f.Q*f.Q,axis=1);ours=(f.B*(f.e*f.e/(1-h))[:,None]).T@f.B;assert np.allclose(hc2,ours,rtol=1e-6,atol=1e-10);CHECK.append(dict(model='Actual primary singleton-cluster HC2 identity',full_OLS_beta_max_difference=0,full_OLS_residual_max_difference=0,CR2_working_model_identity_error=np.max(np.abs(hc2-ours))))
 for block in ['Baseline','RQ interaction','Income interaction']:
  try:
   z=core(eligible(d),True);f=Fit(z,block=block,inference=False);X=np.column_stack([f.D,f.V/f.norm]);positive=z[Y].gt(0).to_numpy();_,rank_pos,_=basis(X[positive]);_,rank_all,_=basis(X);sep='No linear Poisson separation: positive-outcome design has full column rank' if rank_pos==rank_all else 'Positive design rank deficient; separation requires LP'
   if rank_pos<rank_all:
    lp=linprog(X[~positive].sum(axis=0),A_ub=X[~positive],b_ub=np.zeros((~positive).sum()),A_eq=X[positive],b_eq=np.zeros(positive.sum()),bounds=[(-1,1)]*X.shape[1],method='highs')
    if not lp.success or lp.fun<-1e-7:raise RuntimeError('PPML unresolved/separated: '+str(lp.message)+' objective '+str(lp.fun))
    sep='LP found no separating direction within normalized bounded search'
   with warnings.catch_warnings(record=True) as ww:
    res=sm.GLM(z[Y].to_numpy(),X,family=sm.families.Poisson()).fit(maxiter=200,tol=1e-10,cov_type='cluster',cov_kwds={'groups':z.country_code.to_numpy(),'use_correction':True})
   if not res.converged:raise RuntimeError('Poisson IRLS failed to converge')
   beta=res.params[-f.p:]/f.norm;vc=res.cov_params()[-f.p:,-f.p:]/f.norm[:,None]/f.norm[None,:];score=X.T@(z[Y].to_numpy()-res.mu);rel_score=np.max(np.abs(score))/max(1,z[Y].sum());assert rel_score<1e-6
   for i,n in enumerate(f.names):
    se=np.sqrt(vc[i,i]);crit=stats.t.ppf(.975,f.G-1);PPML.append(dict(model='PPML '+block,term=n,estimate=beta[i],SE_cluster_CR1=se,reference_df=f.G-1,p_value=2*stats.t.sf(abs(beta[i]/se),f.G-1),CI95_low=beta[i]-crit*se,CI95_high=beta[i]+crit*se,n=f.N,clusters=f.G,years=f.year_count,rank=rank_all,residual_df=f.N-rank_all,converged=res.converged,iterations=res.fit_history.get('iteration'),positive_design_rank=rank_pos,separation_check=sep,max_scaled_score=rel_score,minimum_mu=float(res.mu.min()),maximum_mu=float(res.mu.max()),deviance=res.deviance,inference='Country cluster CR1, t(G−1); nonlinear CR2/Satterthwaite unavailable',warnings='; '.join(str(w.message) for w in ww)))
  except Exception as e:FAIL.append(dict(model='PPML '+block,procedure='PPML',status='FAILED / not interpretable',reason=repr(e)));print('PPML failure',repr(e),flush=True)
for n,r in [('coefficients',COEF),('model_metadata',META),('marginal_effects_and_contrasts',CONTR),('joint_tests',JOINT),('sample_membership',MEM),('CR2_covariances',COV),('wild_bootstrap',BOOT),('LOCO_coefficients',COUNTRY),('observation_influence',OBS),('implementation_checks',CHECK),('group_support',SUPPORT),('PPML_results',PPML),('failures',FAIL)]:save(pd.DataFrame(r),n)
assert sha(SOURCE)==source_hash
(EV/'run_receipt.json').write_text(json.dumps(dict(source=str(SOURCE),source_sha256=source_hash,source_unchanged=True,python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,seed=20261016,bootstrap_draws=9999,OLS_models=len(META),bootstrap_hypotheses=len(BOOT),PPML_models=len(set(r['model'] for r in PPML)),failed_procedures=FAIL,fsQCA_run=False),indent=2))
print('COMPLETE',len(META),'OLS models',len(BOOT),'bootstraps',len(FAIL),'failures',flush=True)
