from pathlib import Path
import sys,json,itertools,hashlib,importlib.metadata,platform
import numpy as np,pandas as pd
from scipy import stats
from threadpoolctl import threadpool_limits
import statsmodels.api as sm
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0,str((Path(__file__).resolve().parents[1]/'econometric_results')))
from econometric_engine import *
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'econometric_results';EV=OUT/'audit_evidence';FIG=OUT/'econometric_figures';read=lambda n:pd.read_csv(EV/(n+'.csv'));save=lambda x,n:x.to_csv(EV/(n+'.csv'),index=False)
d=pd.read_csv(ROOT/'reports/final_model_selection_20261016/four_income_group_final_dataset.csv');a=core(eligible(d));coef=read('coefficients');meta=read('model_metadata');contr=read('marginal_effects_and_contrasts');joint=read('joint_tests');boot=read('wild_bootstrap');loco=read('LOCO_coefficients');obs=read('observation_influence');support=read('group_support');ppml=read('PPML_results')
# Every planned procedure retained; inference methods are named, never substituted silently.
boot['p_Holm_5_declared_bootstrap_tests']=multipletests(boot.p_wild_restricted,method='holm')[1];save(boot,'wild_bootstrap')
contr.loc[contr.model.eq('Primary Combined')&contr.evaluation.eq('FD slope at governance'),'income_group']='Low income reference';save(contr,'marginal_effects_and_contrasts')
desc=[];variables=[Y,FD,RQ,'log_real_GDP_pc','log_population','trade_pct_gdp','inflation_cpi_pct','findata_official','imf_official']
for scope,z in [('Primary FE sample',a),('All published financing rows',d[d[Y].notna()])]:
 for group in ['Pooled']+G:
  zz=z if group=='Pooled' else z[z.group.eq(group)]
  for c in variables+['log1p_financing']:
   s=np.log1p(zz[Y]) if c=='log1p_financing' else zz[c];v=s.dropna();desc.append(dict(scope=scope,income_group=group,variable=c,rows=len(zz),nonmissing=len(v),missing=int(s.isna().sum()),countries=zz.loc[s.notna(),'country_code'].nunique(),mean=v.mean(),SD=v.std(),minimum=v.min(),p10=v.quantile(.1),median=v.median(),p90=v.quantile(.9),maximum=v.max(),skewness=v.skew(),numeric_zeros=int(v.eq(0).sum())))
save(pd.DataFrame(desc),'descriptive_statistics')
largest=d[d[Y].notna()].nlargest(30,Y)[['country_code','country_name','year','group',Y,'sample_B_Bank_credit']];save(largest,'largest_published_commitments')
finflags=d.groupby(['year','group']).agg(economies=('country_code','size'),published_financing=(Y,'count'),primary_FE_rows=('sample_B_Bank_credit','sum'),FinData_with_any_imputation=('imputed_components',lambda s:int((s>0).sum()))).reset_index();save(finflags,'availability_and_provenance')
# Curves and all primary model VIFs.
curves=[];vifs=[];primary={}
with threadpool_limits(limits=1):
 for block in ['Baseline','RQ interaction','Income interaction','Combined']:
  f=Fit(a,block=block);primary[block]=f
  vifs.extend(dict(model='Primary '+block,term=n,VIF_after_FE=v) for n,v in zip(f.names,f.VIF))
 f=primary['RQ interaction']
 for rq in np.linspace(a[RQ].min(),a[RQ].max(),61):curves.append(dict(RQ=rq,**f.contrast(f.unit('FD')+(rq-f.centers[1])*f.unit('FD_x_Governance'),'FD marginal slope')))
 # PPML contrasts use the explicitly labelled nonlinear CR1 inference.
 ppmarg=[];ppjoint=[];ppcov=[];ppinfluence=[]
 for block in ['Baseline','RQ interaction','Income interaction']:
  z=core(eligible(d),True);ff=Fit(z,block=block,inference=False);X=np.column_stack([ff.D,ff.V/ff.norm]);res=sm.GLM(z[Y].to_numpy(),X,family=sm.families.Poisson()).fit(maxiter=200,tol=1e-10,cov_type='cluster',cov_kwds={'groups':z.country_code.to_numpy(),'use_correction':True});b=res.params[-ff.p:]/ff.norm;V=res.cov_params()[-ff.p:,-ff.p:]/ff.norm[:,None]/ff.norm[None,:]
  for i,n1 in enumerate(ff.names):
   for j,n2 in enumerate(ff.names):ppcov.append(dict(model='PPML '+block,term1=n1,term2=n2,covariance_cluster_CR1=V[i,j]))
  def contrast(c,label,kind):
   est=float(c@b);se=np.sqrt(c@V@c);df=ff.G-1;t=est/se;crit=stats.t.ppf(.975,df);return dict(model='PPML '+block,kind=kind,term=label,estimate=est,SE_cluster_CR1=se,reference_df=df,p_value=2*stats.t.sf(abs(t),df),CI95_low=est-crit*se,CI95_high=est+crit*se,inference='Nonlinear CR1 cluster sandwich, t(G−1), not CR2',percent_mean_change_per_10pp=100*np.expm1(10*est))
  if block=='RQ interaction':
   for rq in np.quantile(z[RQ],[0,.1,.25,.5,.75,.9,1]):ppmarg.append(dict(RQ=rq,**contrast(ff.unit('FD')+(rq-ff.centers[1])*ff.unit('FD_x_Governance'),f'RQ={rq:.5g}','FD marginal slope')))
  if block=='Income interaction':
   C={g:ff.unit('FD')+(ff.unit('FD_x_'+g) if g!=G[0] else 0) for g in G}
   for g,c in C.items():ppmarg.append(contrast(c,g,'Group FD slope'))
   rr=[]
   for g1,g2 in itertools.combinations(G,2):rr.append(contrast(C[g2]-C[g1],g2+' minus '+g1,'Pairwise slope contrast'))
   for r,p in zip(rr,multipletests([r['p_value'] for r in rr],method='holm')[1]):r['p_Holm_6_pair_family']=p;ppmarg.append(r)
   R=np.array([ff.unit('FD_x_'+g) for g in G[1:]]);th=R@b;F=float(th@np.linalg.solve(R@V@R.T,th))/3;ppjoint.append(dict(model='PPML Income interaction',hypothesis='Income interactions jointly zero',F_cluster_CR1=F,df1=3,df2=ff.G-1,p_cluster_CR1=stats.f.sf(F,3,ff.G-1),inference='Nonlinear CR1, approximate F(3,G−1); not HTZ/CR2'))
  if block=='Baseline':
   for country in z.groupby('country_code')[Y].max().nlargest(5).index:
    cut=core(z[z.country_code.ne(country)],True);cf=Fit(cut,inference=False);xx=np.column_stack([cf.D,cf.V/cf.norm]);cr=sm.GLM(cut[Y].to_numpy(),xx,family=sm.families.Poisson()).fit(maxiter=200,tol=1e-10);ppinfluence.append(dict(omitted_country=country,FD_coefficient=cr.params[-cf.p]/cf.norm[0],converged=cr.converged,rows=len(cut),note='Additional PPML influence diagnostic; primary keeps all countries'))
save(pd.DataFrame(curves),'RQ_marginal_curve');save(pd.DataFrame(vifs),'primary_VIF');save(pd.DataFrame(ppmarg),'PPML_marginal_effects');save(pd.DataFrame(ppjoint),'PPML_joint_tests');save(pd.DataFrame(ppcov),'PPML_covariances');save(pd.DataFrame(ppinfluence),'PPML_selected_country_deletions')
lsummary=[]
for (model,term),z in loco.groupby(['model','term']):
 full=z.full_sample_estimate.iloc[0];worst=z.loc[z.change.abs().idxmax()];lsummary.append(dict(model=model,term=term,full_estimate=full,minimum_leave_one_out=z.estimate.min(),maximum_leave_one_out=z.estimate.max(),sign_reversals=int((np.sign(z.estimate)!=np.sign(full)).sum()),largest_change_country=worst.omitted_country,largest_absolute_change=abs(worst.change),countries_deleted=z.omitted_country.nunique()))
save(pd.DataFrame(lsummary),'LOCO_summary')
influence=pd.concat([z.nlargest(20,'Cook_distance_descriptive') for _,z in obs.groupby('model')]);save(influence,'largest_influence_cases')
# Scientific figures, standalone PNG and vector PDF. All intervals pointwise 95%.
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.bbox':'tight'})
def figsave(fig,name):fig.savefig(FIG/(name+'.png'),dpi=300);fig.savefig(FIG/(name+'.pdf'));plt.close(fig)
fig,ax=plt.subplots(1,2,figsize=(10,3.8));ax[0].hist(a[Y],bins=40,color='#2B667D');ax[0].set(xlabel='Commitments, million constant 2022 USD',ylabel='Country-years',title='Published financing: primary sample');ax[1].hist(np.log1p(a[Y]),bins=30,color='#617E43');ax[1].set(xlabel='ln(1 + commitments in million USD)',ylabel='Country-years',title='Prespecified transformation');fig.tight_layout();figsave(fig,'01_financing_distribution')
curve=pd.DataFrame(curves);fig,ax=plt.subplots(figsize=(7,4.5));ax.fill_between(curve.RQ,curve.CI95_low*10,curve.CI95_high*10,color='#B5D7E3',label='Pointwise 95% CR2/Satterthwaite CI');ax.plot(curve.RQ,curve.estimate*10,color='#204B60',label='Conditional FD slope');ax.axhline(0,color='#777',lw=.8);ax.plot(a[RQ],np.full(len(a),ax.get_ylim()[0]),'|',alpha=.12,color='#555');ax.set(xlabel='Observed Regulatory Quality estimate',ylabel='Change in ln(1 + financing) per 10 pp bank credit/GDP',title='Financial depth × Regulatory Quality');ax.legend(fontsize=8);fig.tight_layout();figsave(fig,'02_RQ_marginal_effect')
fig,ax=plt.subplots(figsize=(8,4.5));ps=contr[contr.model.eq('Primary Income interaction')&contr.evaluation.eq('Group FD slope')].set_index('income_group').loc[G];ys=np.arange(4);ax.errorbar(ps.estimate*10,ys,xerr=np.array([ps.estimate-ps.CI95_low,ps.CI95_high-ps.estimate])*10,fmt='o',color='#244F68',capsize=4);ax.axvline(0,color='#777',lw=.8);ax.set(yticks=ys,yticklabels=G,xlabel='Change in ln(1 + financing) per 10 pp bank credit/GDP',title='Historical income-group slopes: pointwise 95% CR2 intervals');ax.invert_yaxis();fig.tight_layout();figsave(fig,'03_income_group_slopes')
fig,axs=plt.subplots(1,3,figsize=(12,3.8));targets=[('Primary Baseline','FD'),('Primary RQ interaction','FD_x_Governance'),('Primary Income interaction','FD_x_High income')]
for ax,(model,term) in zip(axs,targets):
 zz=loco[loco.model.eq(model)&loco.term.eq(term)];ax.hist(zz.estimate,bins=18,color='#547E91');ax.axvline(zz.full_sample_estimate.iloc[0],color='#A94F32',label='Full sample');ax.set(title=term,xlabel='Leave-one-country-out coefficient',ylabel='Country deletions');ax.legend(fontsize=8)
fig.suptitle('All-country deletion stability; no countries removed from primary models',fontsize=11);fig.tight_layout();figsave(fig,'04_country_deletion_stability')
def table(df,name,note):
 priority=['model','evaluation','income_group','term','estimate','SE_CR2','SE_cluster_CR1','Satterthwaite_df','reference_df','CI95_low','CI95_high','p_value','p_Holm_6_pair_family'] if 'estimate' in df else ['model','n','country_clusters','years','model_rank','residual_df','within_R2','hypothesis','HTZ_F','HTZ_df1','HTZ_df2','HTZ_p']
 df=df[[c for c in priority if c in df]+[c for c in df if c not in priority]]
 return dict(name=name,note=note,columns=list(df.columns),rows=json.loads(df.replace([np.inf,-np.inf],np.nan).to_json(orient='records')))
maincols=['term','estimate','SE_CR2','Satterthwaite_df','CI95_low','CI95_high','p_value'];primarycoef=coef[coef.model.eq('Primary Baseline')][maincols];rqcoef=coef[coef.model.eq('Primary RQ interaction')][maincols];ic=coef[coef.model.eq('Primary Income interaction')][maincols];primarymeta=meta[meta.model.str.startswith('Primary')];allbaseline=coef[coef.block.eq('Baseline')];note='Exploratory conditional associations; frozen IRENA 2025, 2014–2023. No causal interpretation.'
payload={
'descriptive_statistics.xlsx':[table(pd.DataFrame(desc),'Descriptive statistics',note+' Two denominators: primary FE sample and all published financing rows.'),table(finflags,'Annual availability','Published zeros retained; financing missing never replaced by zero.'),table(largest,'Largest commitments','All legitimate extreme observations retained; sample flag shows inclusion.'),table(support[support.model.str.startswith('Primary')],'Primary group support','Historical group counts overlap at country level; do not sum group clusters.')],
'baseline_panel_results.xlsx':[table(primarycoef,'Primary baseline','Country/year FE OLS; country-cluster generalized CR2, coefficient-specific Satterthwaite df; 95% intervals.'),table(primarymeta,'Primary model fit','Within R2 uses country-and-year FE residual outcome; df=N−rank differs from inference df.'),table(allbaseline,'All baseline coefficients','Full planned baseline comparisons, including exact matched and separate native supports.'),table(meta[meta.block.eq('Baseline')],'Baseline sample comparison','Raw index units differ; bank depth is not overall financial development.')],
'financial_institutional_interactions.xlsx':[table(rqcoef,'FD and RQ coefficients','Grand-mean centering before constructing product; include main effects.'),table(contr[contr.model.eq('Primary RQ interaction')],'RQ marginal effects','FD marginal slope evaluated only within observed RQ support; CIs are pointwise.'),table(joint[joint.model.eq('Primary RQ interaction')],'RQ joint test','Single-constraint HTZ equals squared Satterthwaite t.'),table(boot[boot.model.eq('Primary RQ interaction')],'RQ wild bootstrap','Restricted-null Rademacher bootstrap; 9,999 draws, seed 20261016; CR2 studentization.'),table(coef[coef.governance.ne(RQ)],'Alternative governance','Each governance indicator enters separately; no specification selected by significance.')],
'income_group_heterogeneity.xlsx':[table(ic,'Income coefficients','Low income reference; historical group main effects retained when identifiable.'),table(contr[contr.model.eq('Primary Income interaction')],'Group slopes and contrasts','Pairwise six-test Holm adjustment; intervals are pointwise, not simultaneous.'),table(joint[joint.model.eq('Primary Income interaction')],'Income joint test','Generalized CR2, small-sample HTZ joint test.'),table(boot[boot.model.eq('Primary Income interaction')],'Income wild bootstrap','Four planned income hypotheses; raw and five-hypothesis adjusted bootstrap p-values both retained.'),table(coef[coef.model.eq('Fixed 2014 income')|coef.model.eq('Stable income Income interaction')],'Classification sensitivity','Fixed 2014 and stable-income restrictions answer different questions; group main effects may be absorbed.'),table(contr[contr.model.isin(['Fixed 2014 income','Stable income Income interaction'])],'Sensitivity slopes','Independent case support recomputed after each restriction.'),table(pd.DataFrame(vifs),'Primary multicollinearity','VIF after removing country and year FE; matrices keep all identified slopes.')],
'robustness_and_inference.xlsx':[table(meta,'All model specifications','All 60 prespecified OLS fits, rank, sample, fit and FinData provenance; no selective significance filtering.'),table(coef[~coef.model.str.startswith('Primary')],'All sensitivity coefficients','Generalized CR2 and coefficient-specific Satterthwaite inference for every OLS sensitivity.'),table(joint,'All joint tests','HTZ joint inference; not naive cluster-count degrees of freedom.'),table(boot,'All bootstrap tests','9,999 restricted-null draws per hypothesis; Monte Carlo SE reported; seed fixed.'),table(pd.DataFrame(lsummary),'Country deletion summary','All 137 primary countries deleted one at a time for each primary block; parent centers fixed.'),table(loco,'Country deletion detail','Coefficient stability diagnostic; no refitted confidence intervals and no automatic country exclusion.'),table(influence,'Influential observations','Cook distance is classical design/residual diagnostic; DFBETA scaled by CR2 SE, not a clustered hypothesis test.'),table(ppml,'PPML coefficients','Converged IRLS with separation checks; country CR1 t(G−1). Nonlinear CR2/Satterthwaite unavailable.'),table(pd.DataFrame(ppmarg),'PPML marginal effects','Multiplicative conditional mean, distinct from log1p outcome; CR1 inference explicitly labelled.'),table(pd.DataFrame(ppjoint),'PPML joint test','Approximate CR1 F test, not generalized CR2/HTZ.'),table(pd.DataFrame(ppinfluence),'PPML influence sensitivity','Five countries with largest single commitments; distribution-based diagnostic, no primary exclusion.'),table(read('implementation_checks'),'Implementation validation','Every OLS fit checked against full dummy OLS; exact CR2 working-model identity and HC2 limit verified.')]
}
(EV/'workbook_payload.json').write_text(json.dumps(payload,ensure_ascii=False,allow_nan=False),encoding='utf-8')
# Freeze package versions for exact reruns; optional reference reader is not an estimator dependency.
packages=['numpy','pandas','scipy','statsmodels','threadpoolctl','matplotlib'];(OUT/'requirements.txt').write_text('\n'.join(f'{p}=={importlib.metadata.version(p)}' for p in packages)+'\n',encoding='utf-8');(OUT/'requirements-validation.txt').write_text('rdata==1.1.0\nopenpyxl=='+importlib.metadata.version('openpyxl')+'\nPillow=='+importlib.metadata.version('Pillow')+'\n',encoding='utf-8')
print('Prepared',sum(map(len,payload.values())),'workbook sheets; figures',len(list(FIG.glob('*.png'))),flush=True)
