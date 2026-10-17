import json,os
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.stats import spearmanr, chi2_contingency
import argparse
parser=argparse.ArgumentParser(description='Explore IRENA record visibility and 2025–2026 data-vintage sensitivity on authorized frozen inputs.')
parser.add_argument('--panel',required=True,help='Private archive input: four_income_group_final_dataset.csv')
parser.add_argument('--vintage',required=True,help='Private archive input: irena_vintage_overlap.csv')
parser.add_argument('--outdir',default='results',help='Output folder for derived tables and receipt')
args=parser.parse_args()
OUT=os.path.abspath(args.outdir);os.makedirs(OUT,exist_ok=True)
panel=pd.read_csv(args.panel)
v=pd.read_csv(args.vintage)
assert len(panel)==2170 and len(v)==1953 and not panel.duplicated(['country_code','year']).any()
assert not v.duplicated(['country_code','year']).any()
panel['reported']=panel['irena2025_pdf'].notna().astype(int)
panel['bank']=pd.to_numeric(panel['FD.AST.PRVT.GD.ZS'],errors='coerce')
panel['rq']=pd.to_numeric(panel['regulatory_quality_estimate'],errors='coerce')
panel['lngdp']=np.log(panel['gdp_per_capita_constant_2015_usd'].where(panel['gdp_per_capita_constant_2015_usd']>0))
panel['lnpop']=np.log(panel.population.where(panel.population>0)
)
panel['trade']=panel.trade_pct_gdp
panel['inflation']=panel.inflation_cpi_pct
panel['year']=panel.year.astype(int)
panel['income']=pd.Categorical(panel.income_group,categories=['Low income','Lower middle income','Upper middle income','High income','Unclassified'])
summary=panel.groupby('income',observed=True).agg(grid=('reported','size'),observed=('reported','sum'),missing=('reported',lambda x:(1-x).sum())).reset_index();summary['observed_pct']=summary.observed/summary.grid*100;summary['missing_pct']=summary.missing/summary.grid*100
summary.to_csv(OUT+'/visibility_coverage.csv',index=False)
from statsmodels.stats.proportion import proportion_confint
summary['coverage_ci_low'],summary['coverage_ci_high']=proportion_confint(summary.observed,summary.grid,alpha=.05,method='wilson')
summary.to_csv(OUT+'/visibility_coverage.csv',index=False)
cc=panel[(panel.income_group!='Unclassified')].dropna(subset=['bank','rq','lngdp','lnpop','trade','inflation']).copy()
print('coverage subset',cc.shape,cc.reported.mean(),cc.groupby('income',observed=True).reported.agg(['size','mean']).to_string())
cc['income']=cc['income'].cat.remove_unused_categories()
cc['bank10']=cc.bank/10;cc['rq01']=cc.rq/0.1
formula='reported ~ bank10 + rq01 + lngdp + lnpop + trade + inflation + C(income, Treatment(reference="Low income")) + C(year)'
model=smf.glm(formula=formula,data=cc,family=sm.families.Binomial()).fit(cov_type='cluster',cov_kwds={'groups':cc.country_code},maxiter=300)
print('GLM convergence',model.converged,model.params[['bank10','rq01']])
res=pd.DataFrame({'term':model.params.index,'beta_logit':model.params.values,'se_country_cluster_CR1':model.bse.values,'p_cluster':model.pvalues.values,'CI_low':model.conf_int()[0].values,'CI_high':model.conf_int()[1].values})
res['odds_ratio']=np.exp(res.beta_logit)
res['OR_low']=np.exp(res.CI_low);res['OR_high']=np.exp(res.CI_high)
res.to_csv(OUT+'/visibility_logit.csv',index=False)
ll_null=smf.glm('reported ~ 1',data=cc,family=sm.families.Binomial()).fit().llf
modelmeta={'N':len(cc),'countries':cc.country_code.nunique(),'observed':int(cc.reported.sum()),'missing':int((1-cc.reported).sum()),'model':'binomial GLM, time FE and income group indicators, country clustered CR1','loglik':float(model.llf),'pseudo_R2_mcfadden':float(1-model.llf/ll_null),'scale_bank':'10 percentage points','scale_RQ':'0.1 WGI estimate point'}
probgroup=cc.assign(pred=model.predict(cc)).groupby('income',observed=True).agg(n=('reported','size'),observed_rate=('reported','mean'),fitted_rate=('pred','mean')).reset_index()
probgroup.to_csv(OUT+'/visibility_predictions.csv',index=False)
switch=cc.groupby('country_code').reported.nunique();print('outcome switch groups',sum(switch>1),'of',len(switch),'countries')
lpm=smf.ols('reported ~ bank10 + rq01 + lngdp + lnpop + trade + inflation + C(country_code) + C(year)',data=cc).fit(cov_type='cluster',cov_kwds={'groups':cc.country_code})
lpmres=pd.DataFrame({'term':lpm.params.index,'coef':lpm.params.values,'se_CR1_country':lpm.bse.values,'p':lpm.pvalues.values,'CI_low':lpm.conf_int()[0].values,'CI_high':lpm.conf_int()[1].values})
lpmres.loc[lpmres.term.isin(['bank10','rq01','lngdp','lnpop','trade','inflation'])].to_csv(OUT+'/visibility_within_lpm.csv',index=False)
modelmeta.update({'lpm_N':len(cc),'lpm_clusters':cc.country_code.nunique(),'switching_countries':int(sum(switch>1))})
v=v.rename(columns={'irena2025_public_musd2022':'finance25','renewable_commitments_million_usd_2023':'finance26'})
v['state']=np.select([v.finance25.notna()&v.finance26.notna(),v.finance25.notna()&v.finance26.isna(),v.finance25.isna()&v.finance26.notna()],['both','2025 only','2026 only'],default='neither')
states=v.groupby(['group','state']).size().unstack(fill_value=0).reset_index()
states['n']=states[[c for c in states if c!='group']].sum(axis=1)
states.to_csv(OUT+'/vintage_coverage.csv',index=False)
trend=v.groupby(['year','state']).size().unstack(fill_value=0).reset_index();trend.to_csv(OUT+'/vintage_by_year.csv',index=False)
both=v.dropna(subset=['finance25','finance26']).copy()
positive=both[(both.finance25>0)&(both.finance26>0)].copy()
cor=spearmanr(positive.finance25,positive.finance26)
print('rank overall',len(positive),cor)
rank_by_year=[]
for year,g in positive.groupby('year'):
 if len(g)>20:
  s=spearmanr(g.finance25,g.finance26)
  rank_by_year.append({'year':year,'n':len(g),'rho':s.statistic,'p':s.pvalue})
pd.DataFrame(rank_by_year).to_csv(OUT+'/vintage_rank_by_year.csv',index=False)
vm=both.merge(panel[['country_code','year','bank','rq','lngdp','lnpop','trade','inflation']],on=['country_code','year'],validate='one_to_one')
vm=vm.dropna(subset=['bank','rq','lngdp','lnpop','trade','inflation']);print('matched regression N',len(vm),'clusters',vm.country_code.nunique())
vm['y25']=np.log1p(vm.finance25);vm['y26']=np.log1p(vm.finance26)
matched=[]
for key in ['y25','y26']:
 m=smf.ols(f'{key} ~ bank + rq + lngdp + lnpop + trade + inflation + C(country_code)+C(year)',data=vm).fit(cov_type='cluster',cov_kwds={'groups':vm.country_code})
 for term in ['bank','rq']:
  c=m.conf_int().loc[term]
  matched.append({'vintage':key,'n':len(vm),'clusters':vm.country_code.nunique(),'term':term,'estimate':m.params[term],'se_CR1_country':m.bse[term],'CI_low':c[0],'CI_high':c[1],'p':m.pvalues[term]})
pd.DataFrame(matched).to_csv(OUT+'/matched_vintage_FE.csv',index=False)
ratio=positive.finance26/positive.finance25
vmeta={'grid_2015_2023':int(len(v)),'2025_observed':int(v.finance25.notna().sum()),'2026_observed':int(v.finance26.notna().sum()),'both':int(len(both)),'2025_only':int((v.state=='2025 only').sum()),'2026_only':int((v.state=='2026 only').sum()),'neither':int((v.state=='neither').sum()),'positive_both':int(len(positive)),'spearman_rank_positive_both':float(cor.statistic),'positive_median_native_base_ratio':float(ratio.median()),'different_price_bases':'2025 constant 2022 USD vs 2026 constant 2023 USD; no conversion performed','2014_excluded_from_2026_overlap':True,'matched_FE_sample':int(len(vm)),'matched_FE_countries':int(vm.country_code.nunique())}
with open(OUT+'/analysis_receipt.json','w') as f:json.dump({'visibility':modelmeta,'vintage':vmeta},f,indent=2)
print('\nMAIN RESULTS LOGIT:');print(res[res.term.str.contains('bank10|rq01|income')].to_string(index=False))
print('\nFE LPM');print(lpmres[lpmres.term.isin(['bank10','rq01'])].to_string(index=False))
print('\nVINTAGE');print(states.to_string(index=False));print(vmeta)
print('\nMATCHED FE');print(pd.DataFrame(matched).to_string(index=False))
