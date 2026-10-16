from planb_setup import *
import pandas as pd,numpy as np,platform,sys
from scipy.linalg import qr
from threadpoolctl import threadpool_limits
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
d=pd.read_csv(SOURCE);co=pd.read_csv(EV/'coefficients.csv');meta=pd.read_csv(EV/'model_plan.csv').merge(pd.read_csv(EV/'R_model_metadata.csv'),on=['model','n','clusters']).merge(pd.read_csv(EV/'Python_model_metadata.csv'),on=['model','n','clusters','rank','residual_df'],suffixes=('_R','_Python'));supp=pd.read_csv(EV/'income_support.csv');win=pd.read_csv(EV/'window_eligibility.csv');ind=pd.read_csv(EV/'independent_validation.csv');loco=pd.read_csv(EV/'country_deletion_diagnostics.csv');ct=pd.read_csv(EV/'within_between_contrasts.csv');jt=pd.read_csv(EV/'R_joint_tests.csv')
G=['Low income','Lower middle income','Upper middle income','High income'];X=['FD','RQ','lnGDP','lnPop','Trade','Inflation'];a=pd.read_csv(EV/'annual_TWFE_absolute_data.csv')
aggcheck=pd.read_csv(EV/'R_independent_aggregation_checks.csv');aggcheck=aggcheck.rename(columns={'absolute_error':'max_absolute_difference'});aggcheck['passed']=aggcheck.max_absolute_difference<=aggcheck.tolerance;ind=pd.concat([ind,aggcheck],ignore_index=True)
coll=[]
with threadpool_limits(limits=1):
 for model in meta.model:
  z=pd.read_csv(EV/(model+'_data.csv'));D=pd.read_csv(EV/(model+'_design.csv'));A=D.to_numpy().copy();r=meta[meta.model==model].iloc[0];ratio=r.tau2/r.sigma2 if r.method=='CRE' else 0
  for code in z.country_code.unique():
   ix=np.flatnonzero(z.country_code==code);c=(1/np.sqrt(1+ratio*len(ix))-1)/len(ix);A[ix]+=c*A[ix].sum(axis=0)
  focus=X if r.method=='TWFE' else ['W_'+x for x in X]+['B_'+x for x in X];idx=[D.columns.get_loc(x) for x in focus];nuis=[i for i in range(D.shape[1]) if i not in idx];Q,_=qr(A[:,nuis],mode='economic');V=A[:,idx]-Q@(Q.T@A[:,idx]);V/=np.linalg.norm(V,axis=0);vif=np.diag(np.linalg.inv(V.T@V))
  for t,v in zip(focus,vif):coll.append(dict(model=model,term=t,VIF_conditional_slopes=float(v),definition='Slope block residualized on included fixed/time/profile nuisance columns, then unit-norm scaled; not raw full-design inverse-Gram diagonal.'))
coll=pd.DataFrame(coll);coll.to_csv(EV/'conditional_slope_collinearity.csv',index=False);meta=meta.rename(columns={'max_VIF':'max_full_design_inverse_Gram_diagonal'}).merge(coll.groupby('model').VIF_conditional_slopes.max().rename('max_conditional_slope_VIF'),on='model');meta['last_calendar_year']=np.where(meta.model.str.contains('window3|matched3'),2022,2023)
d['Yabs']=d.irena2025_pdf;d['Ypc']=d.Yabs*1e6/d.population;d['FD']=d['FD.AST.PRVT.GD.ZS'];d['RQ']=d.regulatory_quality_estimate;d['lnGDP']=np.log(d.gdp_per_capita_constant_2015_usd.where(d.gdp_per_capita_constant_2015_usd>0));d['lnPop']=np.log(d.population.where(d.population>0));d['Trade']=d.trade_pct_gdp;d['Inflation']=d.inflation_cpi_pct
for z in [d,a]:z['log1p_absolute']=np.log1p(z.Yabs);z['log1p_per_capita']=np.log1p(z.Ypc)
variation=[]
for label,z,variables in [('Native complete cases',a,['Yabs','Ypc','log1p_absolute','log1p_per_capita']+X),('All observed financing',d[d.Yabs.notna()],['Yabs','Ypc','log1p_absolute','log1p_per_capita'])]:
 for v in variables:
  zz=z[['country_code','year',v]].dropna();x=zz[v];b=x.groupby(zz.country_code).transform('mean');w=x-b;ssw=float(w@w);ssb=float(((b-x.mean())**2).sum());sst=float(((x-x.mean())**2).sum());assert abs(sst-ssw-ssb)<1e-8*max(sst,1)
  variation.append(dict(sample=label,variable=v,observations=len(zz),countries=zz.country_code.nunique(),overall_SD=x.std(),within_SD=np.sqrt(ssw/(len(zz)-zz.country_code.nunique())),country_mean_SD=zz.groupby('country_code')[v].mean().std(),within_SS_share=ssw/sst,between_SS_share=ssb/sst,SS_identity_error=abs(sst-ssw-ssb)))
variation=pd.DataFrame(variation)
coverage=[];higroupmissing=[];selection=[]
for (year,g),z in d.groupby(['year','group'],dropna=False):
 if g not in G:g='Unclassified'
 observed=z.Yabs.notna();eligible=z.sample_B_Bank_credit.astype(bool)
 coverage.append(dict(year=int(year),income_group=g,economies=len(z),observed_financing=int(observed.sum()),missing_financing=int((~observed).sum()),missing_share=float((~observed).mean()),published_numeric_zero=int(z.financing_published_numeric_zero.sum()),primary_model_observations=int(eligible.sum()),observed_not_primary=int((observed&~eligible).sum())))
 if g=='High income':
  for r in z.itertuples():higroupmissing.append(dict(year=int(year),country_code=r.country_code,country_name=r.country_name,financing_observed=pd.notna(r.Yabs),primary_eligible=r.sample_B_Bank_credit,financing_musd2022=r.Yabs,entity_type=r.entity_type,reason=r.financing_missing_reason))
 for observedflag in [True,False]:
  q=z[observed==observedflag]
  for v in ['FD','RQ','lnGDP','lnPop']:selection.append(dict(year=int(year),income_group=g,financing_observed=observedflag,variable=v,economies=len(q),nonmissing_covariate=int(q[v].notna().sum()),median=q[v].median()))
coverage=pd.DataFrame(coverage);hi=pd.DataFrame(higroupmissing);selection=pd.DataFrame(selection)
assert coverage.economies.sum()==2170 and coverage.observed_financing.sum()==d.Yabs.notna().sum()
hicomp=[]
for code,z in hi.groupby('country_code'):
 r=dict(country_code=code,country_name=z.iloc[0].country_name,historical_high_income_years=len(z),observed_high_income_years=int(z.financing_observed.sum()),primary_FE_high_income_years=int(z.primary_eligible.sum()))
 for model,key in [('balanced_TWFE','balanced_high_income_years'),('window2_TWFE_absolute','two_year_high_income_windows'),('window3_TWFE_absolute','three_year_high_income_windows')]:
  zz=pd.read_csv(EV/(model+'_data.csv'));r[key]=int(((zz.country_code==code)&(zz.group=='High income')).sum())
 hicomp.append(r)
hicomp=pd.DataFrame(hicomp);hicomp.to_csv(EV/'high_income_model_composition.csv',index=False)
concentration=[];annual_peaks=[];countryrows=[]
def conc(label,x):
 x=x.dropna();total=x.sum();s=x/total if total>0 else x*np.nan;sort=s.sort_values(ascending=False)
 concentration.append(dict(scope=label,records=len(x),observed_amount_musd2022=total,largest_share=sort.iloc[0] if len(x) else np.nan,top5_share=sort.head(5).sum(),top10_share=sort.head(10).sum(),HHI=float((s*s).sum()),effective_equal_amount_units=1/float((s*s).sum()) if total>0 else np.nan,definition='Shares of available recorded amounts only; unreported financing excluded, never set zero'))
obs=d[d.Yabs.notna()].copy();conc('All observed country-year records 2014–2023',obs.Yabs);conc('Country observed-period sums (unequal reporting histories)',obs.groupby('country_code').Yabs.sum())
for year,z in obs.groupby('year'):
 conc('Observed recipients in '+str(year),z.Yabs)
 for _,r in z.nlargest(5,'Yabs').iterrows():annual_peaks.append(dict(year=int(year),country_code=r.country_code,country_name=r.country_name,income_group=r.group,amount_musd2022=r.Yabs,share_observed_year_total=r.Yabs/z.Yabs.sum()))
for code,z in obs.groupby('country_code'):
 total=z.Yabs.sum();mean=z.Yabs.mean();countryrows.append(dict(country_code=code,country_name=z.iloc[0].country_name,observed_years=len(z),recorded_numeric_zero=int(z.Yabs.eq(0).sum()),observed_total_musd2022=total,annual_mean_musd2022=mean,largest_year_share=z.Yabs.max()/total if total>0 else np.nan,annual_CV=z.Yabs.std()/mean if mean>0 else np.nan,peak_year=int(z.loc[z.Yabs.idxmax(),'year']),annual_per_capita_mean=z.Ypc.mean(),interpretation='Observed timing proxy, not project-level concentration; incomplete years not assumed zero'))
concentration=pd.DataFrame(concentration);countryrows=pd.DataFrame(countryrows);annual_peaks=pd.DataFrame(annual_peaks)
with threadpool_limits(limits=1):
 D=pd.concat([pd.get_dummies(a.country_code,dtype=float),pd.get_dummies(a.year,drop_first=True,dtype=float)],axis=1).to_numpy();Q,_=qr(D,mode='economic');res=a[X].to_numpy()-Q@(Q.T@a[X].to_numpy())
within_corr=pd.DataFrame(res,columns=X).corr().rename_axis('variable').reset_index();between_corr=a.groupby('country_code')[X].mean().corr().rename_axis('variable').reset_index()
timing=[]
for v in ['log1p_absolute','log1p_per_capita','FD','RQ']:
 z=a[['country_code','year',v]].copy();z['demeaned']=z[v]-z.groupby('country_code')[v].transform('mean');p=z[['country_code','year','demeaned']].copy();p.year+=1;p=p.rename(columns={'demeaned':'previous_demeaned'});j=z.merge(p,on=['country_code','year']);timing.append(dict(variable=v,consecutive_pairs=len(j),countries=j.country_code.nunique(),demeaned_adjacent_correlation=j.demeaned.corr(j.previous_demeaned),definition='Descriptive observed adjacent pairs; country demeaned; no inference of project renewal or reporting continuity'))
timing=pd.DataFrame(timing)
# Native-vs-matched-vs-aggregated comparisons keep the estimand and support explicit.
comparison=[]
for h in [2,3]:
 for out in ['absolute','per_capita']:
  for kind,name in [('Annual native','annual_TWFE_'+out),('Annual matched',f'matched{h}_annual_TWFE_'+out),('Window',f'window{h}_TWFE_'+out)]:
   z=co[(co.model==name)&co.term.isin(['FD','RQ'])].copy();z['horizon']=h;z['comparison']=kind;z['outcome']=out;comparison.extend(z.to_dict('records'))
comparison=pd.DataFrame(comparison)
limits=pd.DataFrame([
 dict(topic='Project concentration',available='Country-year published renewable totals',limit='No transaction IDs, project counts, size, cancellation or lender allocation; cannot calculate actual project concentration or distinguish one large commitment from several.'),
 dict(topic='Provider dependence',available='One frozen selected-public recipient vintage',limit='Exact provider-year roster, reporting breaks and stable institutional perimeter not recovered. Cannot estimate provider fixed effects, leave-provider-out results or institution shares.'),
 dict(topic='Commitment timing',available='Annual country amounts; observed-year peaks and complete windows',limit='Approval/recognition timing can be lumpy; smoothing cannot distinguish reporting revisions, cancellations, refinancing, actual disbursement or private investment.'),
 dict(topic='Currency/price base',available='Specific finance-table headings: million2022USD; unchanged source',limit='Generic publication note conflicts with specific table headings; constant2022 operational convention is provisional, not new institutional confirmation. No price conversion or later-vintage substitution.'),
 dict(topic='Missing financing',available='Published zeros and missing-reason flags',limit='Absence/nonreporting/structural recipient scope not distinguishable from zero commitments; missing never assigned zero; incomplete-window means excluded.'),
 dict(topic='CRE interpretation',available='Within/between X and time-profile means; year FE',limit='Algebraic within equivalence and CR2 do not solve endogenous selection, omitted time-varying factors, reverse causality or financing anticipation.'),
 dict(topic='Income comparisons',available='Historical annual groups, window-start groups, actual country support',limit='No new group interactions or separate income regressions; group-country counts overlap when classifications change, and do not sum to pooled clusters.'),
 dict(topic='Inference scope',available='CR2/Satterthwaite for all19 linear models; HTZ joint tests',limit='No new wild bootstrap, PPML or fsQCA. Small effective-cluster/df warnings retained. Primary8 Holm applies to p-values; displayed 95% CIs are pointwise.'),
 dict(topic='Country deletions',available='All137 annualCRE and119 two-yearTWFE country deletions',limit='CRE diagnostic holds full-sample REML variance ratio fixed; conditional FGLS coefficients, not refitted variance-component inference or a new primary sample.')])
for name,z in [('variation_decomposition',variation),('coverage_by_income_year',coverage),('high_income_observation_status',hi),('observed_missing_selection',selection),('concentration',concentration),('country_timing',countryrows),('annual_largest_recipients',annual_peaks),('within_correlations',within_corr),('between_correlations',between_corr),('adjacent_timing',timing),('annual_window_comparisons',comparison),('measurement_limits',limits)]:z.to_csv(EV/(name+'.csv'),index=False)
plan=json.loads((EV/'preanalysis_plan_frozen.json').read_text());rules=pd.DataFrame([dict(rule=k,definition=('UTC '+v if k.endswith('_utc') else json.dumps(v,ensure_ascii=False) if isinstance(v,(list,dict)) else str(v))) for k,v in plan.items()])
rules=pd.concat([rules,pd.DataFrame([dict(rule='Institutional/source references',definition='IRENA2025 PDF: https://www.irena.org/-/media/Files/IRENA/Agency/Publication/2025/Jul/IRENA_DAT_RE_Statistics_2025.pdf ; WDI private bank credit: https://data.worldbank.org/indicator/FD.AST.PRVT.GD.ZS ; WGI: https://www.worldbank.org/en/publication/worldwide-governance-indicators ; Frozen source definitions and unresolved conflicts: ../reports/final_methodological_clearance.md')])],ignore_index=True)
def rows(z):
 z=z.replace([np.inf,-np.inf],np.nan).astype(object);z=z.where(pd.notna(z),None);return z.to_dict('records')
def sheet(name,z,note):return dict(name=name,note=note,columns=list(z.columns),rows=rows(z))
annual=co[~co.model.str.contains('window|matched')];multi=co[co.model.str.contains('window|matched')]
payload={
 'within_between_mundlak_results.xlsx':[sheet('Primary annual CRE',co[(co.model=='annual_CRE_absolute')],'Within/between components, same1104 observations/137 countries; pointwise CR2 CIs; primary8 Holm p-values where applicable.'),sheet('Annual estimates and sensitivities',annual,'All frozen annual native, balanced-history and lag models; every coefficient retained.'),sheet('Within between contrasts',ct,'Two differences Holm-adjusted as a separate family; no causal interpretation.'),sheet('Mundlak joint tests',jt,'HTZ joint restriction within=between; not a pretest selecting an estimator.'),sheet('Model support and fit',meta[~meta.model.str.contains('window|matched')],'Random-intercept variance, rank, residual df; conditional-slope VIF distinguished from full-design scaling.'),sheet('Conditional slope collinearity',coll,'VIF after nuisance projection; full FE dummy normalization can otherwise exaggerate apparent collinearity.'),sheet('Income support',supp[~supp.model.str.contains('window|matched')],'Historical yearly classifications; group country counts may overlap.'),sheet('Independent verification',ind,'Python whitened GLS/CR2 and concentrated REML profiles independently match R nlme/clubSandwich.'),sheet('Country deletion diagnostics',loco[loco.model=='annual_CRE_absolute'],'Every country omitted once; fixed full-sample variance ratio, conditional FGLS only.')],
 'multi_year_financing_results.xlsx':[sheet('Annual window comparison',comparison,'Compare native, exact matched annual support and prespecified complete windows; outcomes remain distinct.'),sheet('All multi year coefficients',multi,'Two- and three-year windows, two-year CRE and matched annual models; no preferred horizon selected by p-values.'),sheet('Window model support',meta[meta.model.str.contains('window|matched')],'Aggregated counts are country-windows; matched annual counts are country-years.'),sheet('Income support',supp[supp.model.str.contains('window|matched')],'Window group is first-year historical classification; no new group interaction estimation.'),sheet('Complete window eligibility',win,'All required calendar years observed; no incomplete-window means computed.'),sheet('Two year window data',pd.read_csv(EV/'complete_2year_windows.csv'),'Mean annual commitments, not unique-project totals; annual per-capita ratios averaged.'),sheet('Three year window data',pd.read_csv(EV/'complete_3year_windows.csv'),'2014–2016,2017–2019,2020–2022 only; 2023 excluded by frozen design.'),sheet('Country deletion diagnostics',loco[loco.model=='window2_TWFE_absolute'],'All119 country deletions; no trimming or promoted sensitivity.')],
 'financing_measurement_and_selection_audit.xlsx':[sheet('Within between variation',variation,'SS shares exactly decompose total observed variation; country-mean SD uses equal countries.'),sheet('Income year missingness',coverage,'Observed/missing records and actual zeros remain distinct; missing never treated as zero.'),sheet('High income coverage',hi,'All historical high-income economy-years, including nonreporting and ineligible cases.'),sheet('High income composition',hicomp,'Exact recipient-economy retention as complete windows/balanced histories are required.'),sheet('Observed versus missing',selection,'Available-covariate medians only; not evidence that missing financing is zero or ignorable.'),sheet('Recorded concentration',concentration,'Shares of observed amounts; project-level and provider concentration cannot be inferred.'),sheet('Country timing proxies',countryrows,'Peak-year share and annual CV use observed years only; reporting histories differ.'),sheet('Largest annual recipients',annual_peaks,'Five largest recorded recipient amounts per year; records are not projects.'),sheet('Within X correlations',within_corr,'Country/year-residualized X on frozen annual1104 support; no variable selection.'),sheet('Between X correlations',between_corr,'Country-mean X correlations; selected recipients only.'),sheet('Adjacent financing timing',timing,'Observed consecutive annual pairs; descriptive demeaned correlations only.'),sheet('Measurement limits',limits,'Unresolved institutional, project-level, currency-documentation and selection limits.'),sheet('Frozen bounded protocol',rules,'Windows, definitions, models and test families frozen before estimation.')]
}
for sheets in payload.values():
 for s in sheets:
  if s['name']=='Annual estimates and sensitivities':s['name']='Annual model sensitivities'
  assert len(s['name'])<=31,s['name']
(EV/'workbook_payload.json').write_text(json.dumps(payload,allow_nan=False),encoding='utf-8')
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.bbox':'tight'})
def save(fig,name):fig.savefig(FIG/(name+'.png'),dpi=200);fig.savefig(FIG/(name+'.pdf'));plt.close(fig)
fig,axs=plt.subplots(1,2,figsize=(11,4))
for ax,v in zip(axs,['FD','RQ']):
 z=co[(co.model=='annual_CRE_absolute')&co.term.isin(['W_'+v,'B_'+v])];ax.errorbar(z.estimate,[1,0],xerr=[z.estimate-z.CI95_low,z.CI95_high-z.estimate],fmt='o',capsize=4,color='#275d79');ax.axvline(0,color='grey',lw=.8);ax.set(yticks=[1,0],yticklabels=['Within','Between'],xlabel='Coefficient; pointwise95% CR2 interval',title=v+' — selected-public absolute commitments')
fig.suptitle('CRE within/between components — identical annual support');fig.tight_layout();save(fig,'01_within_between')
fig,axs=plt.subplots(1,2,figsize=(12,4.5))
for ax,v in zip(axs,['FD','RQ']):
 z=comparison[(comparison.outcome=='absolute')&(comparison.term==v)];labels=[str(r.horizon)+'yr: '+r.comparison for r in z.itertuples()];pos=np.arange(len(z));ax.errorbar(z.estimate,pos,xerr=[z.estimate-z.CI95_low,z.CI95_high-z.estimate],fmt='o',capsize=3,color='#275d79');ax.axvline(0,color='grey',lw=.8);ax.set(yticks=pos,yticklabels=labels,xlabel='Pointwise95% CR2 interval',title=v+' — absolute financing');ax.invert_yaxis()
fig.tight_layout();save(fig,'02_annual_matched_windows')
fig,ax=plt.subplots(figsize=(8,4.5))
for g in G:
 z=coverage[coverage.income_group==g];ax.plot(z.year,z.missing_share*100,marker='o',label=g)
ax.set(xlabel='Year',ylabel='Financing records missing (%)',ylim=(0,100),title='Historical income-group reporting selection');ax.legend(fontsize=8);save(fig,'03_reporting_selection')
fig,ax=plt.subplots(figsize=(8,4.5));z=variation[(variation['sample']=='Native complete cases')&variation.variable.isin(['Yabs','log1p_absolute','FD','RQ'])];ax.barh(z.variable,z.within_SS_share,label='Within-country',color='#275d79');ax.barh(z.variable,z.between_SS_share,left=z.within_SS_share,label='Between-country',color='#b77b50');ax.set(xlabel='Share of total observed sum of squares',xlim=(0,1),title='Variation decomposition — native1104 support');ax.legend(loc='lower right');save(fig,'04_variation_decomposition')
for n in ['planb_setup.py','planb_prepare.py','planb_estimate.R','planb_run_R.py','planb_verify.py','planb_outputs.py','planb_aggregation_check.R','planb_run_aggregation_check.py']:
 if Path(__file__).parent.resolve()!=SCRIPTS.resolve():shutil.copy2(Path(__file__).parent/n,SCRIPTS/n)
(SCRIPTS/'dependency_versions.json').write_text(json.dumps(dict(Python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,scipy=__import__('scipy').__version__,matplotlib=matplotlib.__version__,R='4.6.1',clubSandwich='0.7.0',nlme='See retained R session information'),indent=2),encoding='utf-8')
print('Outputs prepared',len(payload),'workbooks',sum(map(len,payload.values())),'sheets');print(variation.to_string(index=False));print(concentration.head(2).to_string(index=False))
