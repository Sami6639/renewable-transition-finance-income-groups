from qca_replication_setup import *
import pandas as pd,numpy as np
from datetime import datetime,timezone
SOURCE=ROOT/'reports/final_model_selection_20261016/four_income_group_final_dataset.csv'
d=pd.read_csv(SOURCE);G=['Low income','Lower middle income','Upper middle income','High income'];FD='FD.AST.PRVT.GD.ZS';RQ='regulatory_quality_estimate';GDP='gdp_per_capita_constant_2015_usd';Y='irena2025_pdf'
protocol=json.loads((EV/'proposed_protocol_frozen.json').read_text(encoding='utf-8'));assert protocol['status']=='PROPOSED — NOT APPROVED';assert len(d)==2170 and d.year.max()==2023 and not d.duplicated(['country_code','year']).any()
rows=[];missing=[]
for code,z in d.groupby('country_code'):
 z=z.set_index('year');b=z.reindex([2018,2019,2020]);o=z.reindex([2021,2022,2023]);first=z.iloc[0]
 valid_fd=bool(np.isfinite(b[FD]).all());valid_rq=bool(np.isfinite(b[RQ]).all());valid_gdp=bool(np.isfinite(b[GDP]).all() and b[GDP].gt(0).all());valid_y=bool(np.isfinite(o[Y]).all() and o[Y].ge(0).all());valid_pop=bool(np.isfinite(o.population).all() and o.population.gt(0).all());known=z.loc[2018,'group'] in G;eligible=valid_fd and valid_rq and valid_gdp and valid_y and valid_pop and known
 reasons=[name for name,ok in [('bank_depth_3_years',valid_fd),('RQ_3_years',valid_rq),('positive_real_GDPpc_3_years',valid_gdp),('nonnegative_financing_3_years',valid_y),('positive_outcome_population_3_years',valid_pop),('classified_income_2018',known)] if not ok]
 r=dict(country_code=code,country_name=first.country_name,entity_type=first.entity_type,sovereign_scope=first.sovereign_scope,income_group_2018=z.loc[2018,'group'],eligible=eligible,exclusion_reasons='; '.join(reasons),baseline_FD_complete=valid_fd,baseline_RQ_complete=valid_rq,baseline_GDPpc_complete=valid_gdp,outcome_financing_complete=valid_y,outcome_population_complete=valid_pop,baseline_income_classified=known,stable_income_2018_2023=z.loc[2018:2023,'group'].nunique(dropna=False)==1)
 for yr in range(2018,2024):r['income_group_'+str(yr)]=z.loc[yr,'group']
 for yr in [2018,2019,2020]:
  for c,label in [(FD,'bank_depth'),(RQ,'RQ'),(GDP,'real_GDP_pc')]:r[f'{label}_{yr}']=z.loc[yr,c]
  r[f'bank_depth_obs_status_{yr}']=z.loc[yr,FD+'_obs_status']
 for yr in [2021,2022,2023]:
  r[f'financing_musd2022_{yr}']=z.loc[yr,Y];r[f'financing_missing_{yr}']=pd.isna(z.loc[yr,Y]);r[f'population_{yr}']=z.loc[yr,'population'];r[f'financing_per_capita_{yr}']=z.loc[yr,Y]*1e6/z.loc[yr,'population'] if pd.notna(z.loc[yr,Y]) and z.loc[yr,'population']>0 else np.nan;r[f'financing_published_zero_{yr}']=bool(z.loc[yr,Y]==0)
 r.update(bank_depth_mean_2018_2020=b[FD].mean() if valid_fd else np.nan,RQ_mean_2018_2020=b[RQ].mean() if valid_rq else np.nan,real_GDP_pc_mean_2018_2020=b[GDP].mean() if valid_gdp else np.nan,financing_mean_musd2022_2021_2023=o[Y].mean() if valid_y else np.nan,financing_cumulative_musd2022_2021_2023=o[Y].sum() if valid_y else np.nan,financing_mean_per_capita_2021_2023=(o[Y]*1e6/o.population).mean() if valid_y and valid_pop else np.nan,financing_ratio_of_sums_per_capita=(o[Y].sum()*1e6/o.population.sum()) if valid_y and valid_pop else np.nan,primary_outcome='mean annual selected-public commitments per capita',source_vintage='IRENA Renewable Energy Statistics 2025, frozen PDF',financing_unit='constant 2022 USD; annual totals expressed in millions',GDP_unit='constant 2015 USD per person',calibration_status='Raw cases only; memberships withheld pending protocol approval',publisher_estimation_status='WDI blank obs_status not proof of unestimated reporting')
 rows.append(r)
 for yr in range(2018,2024):
  variables=[FD,RQ,GDP] if yr<=2020 else [Y,'population']
  for c in variables:missing.append(dict(country_code=code,income_group_2018=z.loc[2018,'group'],year=yr,variable=c,missing=bool(pd.isna(z.loc[yr,c]))))
allcases=pd.DataFrame(rows);cases=allcases[allcases.eligible].copy();assert len(cases)==111 and cases.country_code.nunique()==111 and cases.income_group_2018.isin(G).all()
allcases.to_csv(EV/'all_economy_case_eligibility.csv',index=False);cases.to_csv(OUT/'fsqca_case_dataset.csv',index=False);pd.DataFrame(missing).groupby(['income_group_2018','year','variable'],dropna=False).agg(economies=('country_code','size'),missing=('missing','sum')).reset_index().to_csv(EV/'window_missingness.csv',index=False)
coverage=[]
for g in ['Pooled']+G+['Unclassified']:
 z=allcases if g=='Pooled' else allcases[allcases.income_group_2018.eq(g)];c=z[z.eligible];coverage.append(dict(income_group_2018=g,economies=len(z),eligible_cases=len(c),excluded=len(z)-len(c),baseline_FD_missing=int((~z.baseline_FD_complete).sum()),baseline_RQ_missing=int((~z.baseline_RQ_complete).sum()),baseline_GDPpc_missing_or_nonpositive=int((~z.baseline_GDPpc_complete).sum()),financing_window_incomplete=int((~z.outcome_financing_complete).sum()),population_window_incomplete=int((~z.outcome_population_complete).sum()),stable_income_cases=int(c.stable_income_2018_2023.sum())))
pd.DataFrame(coverage).to_csv(EV/'case_coverage_by_income.csv',index=False)
diagnostics=[];near=[];descriptive=[];alternatives=[]
for row in protocol['anchors']:
 var=row['variable'];v=cases[var]
 for kind,c in [('Primary crossover',row['crossover'])]+[('Alternative crossover',c) for c in row['crossover_alternatives']]:
  alternatives.append(dict(variable=var,scenario=kind,full_out=row['full_out'],crossover=c,full_in=row['full_in'],status='Proposed, not calibrated or optimized'))
  for g in ['Pooled']+G:
   z=cases if g=='Pooled' else cases[cases.income_group_2018.eq(g)];x=z[var];band=row['near_crossover_band'];diagnostics.append(dict(variable=var,scenario=kind,crossover=c,income_group_2018=g,cases=len(z),raw_below_crossover=int(x.lt(c).sum()),raw_equal_crossover=int(np.isclose(x,c,rtol=0,atol=1e-10).sum()),raw_above_crossover=int(x.gt(c).sum()),near_crossover=int(x.sub(c).abs().le(band).sum()),band=band,interpretation='Univariate raw-value review only; not memberships or truth-table configurations'))
  for rr in cases[v.sub(c).abs().le(row['near_crossover_band'])].itertuples():near.append(dict(country_code=rr.country_code,country_name=rr.country_name,income_group_2018=rr.income_group_2018,variable=var,raw_value=getattr(rr,var),crossover=c,scenario=kind))
 for g in ['Pooled']+G:
  x=cases[var] if g=='Pooled' else cases.loc[cases.income_group_2018.eq(g),var];descriptive.append(dict(variable=var,income_group_2018=g,cases=len(x),minimum=x.min(),median=x.median(),maximum=x.max(),raw_at_or_below_full_out=int(x.le(row['full_out']).sum()),raw_at_or_above_full_in=int(x.ge(row['full_in']).sum()),warning='Descriptive support examined only after anchors frozen; no anchors changed'))
for name,rows in [('proposed_anchor_alternatives',alternatives),('univariate_crossover_diagnostics',diagnostics),('cases_near_crossovers',near),('condition_outcome_support',descriptive)]:pd.DataFrame(rows).to_csv(EV/(name+'.csv'),index=False)
corr=cases[[a['variable'] for a in protocol['anchors']]].corr();corr.to_csv(EV/'raw_condition_outcome_correlations.csv')
# Contextual selection diagnostics, without fitting a selection model.
select=[]
for (g,elig),z in allcases.groupby(['income_group_2018','eligible'],dropna=False):select.append(dict(income_group_2018=g,eligible=elig,economies=len(z),FD_mean_median=z.bank_depth_mean_2018_2020.median(),RQ_mean_median=z.RQ_mean_2018_2020.median(),GDPpc_mean_median=z.real_GDP_pc_mean_2018_2020.median(),note='Available raw conditions only; not proof of a missingness mechanism'))
pd.DataFrame(select).to_csv(EV/'eligible_excluded_comparison.csv',index=False)
pd.DataFrame([dict(variable=c,unit=next((r['unit'] for r in protocol['anchors'] if r['variable']==c),'See annual source column / flag'),definition='Three actual years required; no available-year means' if c in [r['variable'] for r in protocol['anchors']] else 'Original source value, contextual history or explicit eligibility/provenance flag; annual financing in million constant 2022 USD, annual population in persons',missing=int(cases[c].isna().sum())) for c in cases]).to_csv(EV/'case_variable_dictionary.csv',index=False)
(EV/'case_preparation_receipt.json').write_text(json.dumps(dict(source=str(SOURCE),source_sha256=sha(SOURCE),protocol_sha256=sha(EV/'proposed_protocol_frozen.json'),prepared_utc=datetime.now(timezone.utc).isoformat(),eligible_cases=len(cases),cases_by_income=cases.income_group_2018.value_counts().to_dict(),memberships_calculated=False,truth_table_created=False,minimization_run=False),indent=2))
print(pd.DataFrame(coverage).to_string(index=False));print('No anchors altered; no memberships or configurations calculated.')
