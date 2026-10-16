from stage2_setup import *
import pandas as pd,numpy as np,math,shutil
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
read=lambda n:pd.read_csv(EV/(n+'.csv'),low_memory=False)
coef=read('coefficients');meta=read('model_metadata');boot=read('wild_bootstrap');support=read('group_support');loco=read('LOCO_coefficients');ppml=read('PPML_results');comparisons=read('matched_comparisons');joint=read('joint_tests')
targets={'A primary':'RISE','B primary':'FDI','C primary':'Deficit'}
from statsmodels.stats.multitest import multipletests
coef['p_Holm_sensitivity_family']=np.nan
for family,term in [('A','RISE'),('B','FDI'),('C','Deficit')]:
 names=[family+' asinh inflation',family+' all-X lag']+(['B asinh signed FDI'] if family=='B' else [])
 mask=coef.model.isin(names)&coef.term.eq(term);coef.loc[mask,'p_Holm_sensitivity_family']=multipletests(coef.loc[mask,'p_value'],method='holm')[1]
coef['p_Holm_institutional2']=np.nan;mask=coef.model.isin(['Voice and accountability benchmark','Political stability benchmark'])&coef.term.eq('Governance');coef.loc[mask,'p_Holm_institutional2']=multipletests(coef.loc[mask,'p_value'],method='holm')[1]
ppml['p_Holm_PPML_primary3']=np.nan;mask=np.array([r.model in targets and targets[r.model]==r.term for r in ppml.itertuples()]);ppml.loc[mask,'p_Holm_PPML_primary3']=multipletests(ppml.loc[mask,'p_value'],method='holm')[1]
coef.to_csv(EV/'coefficients.csv',index=False);ppml.to_csv(EV/'PPML_results.csv',index=False)
priority=coef[np.array([r.model in targets and targets[r.model]==r.term for r in coef.itertuples()])].merge(meta[['model','n','clusters','rank','residual_df','within_R2']],on='model',validate='many_to_one').merge(boot[['model','p_wild_restricted','p_wild_Holm_primary3']],on='model',validate='one_to_one')
priority=priority.merge(support[support.income_group.eq('High income')][['model','n','clusters']].rename(columns={'n':'high_income_n','clusters':'high_income_clusters'}),on='model',validate='one_to_one');priority.to_csv(EV/'priority_summary.csv',index=False)
# Independent R verification: not another invocation of the Python engine.
r=read('R_coefficients');ci=read('R_confidence_intervals');rc=r.merge(ci[['model','term','CI_L','CI_U']],on=['model','term'],validate='one_to_one');py=coef[coef.model.isin(list(targets)+['A readiness interaction'])];rep=py.merge(rc,on=['model','term'],validate='one_to_one');assert len(rep)==30
checks=[]
for pcol,rcol,tol in [('estimate','beta',5e-9),('SE_CR2','SE',5e-9),('Satterthwaite_df','df_Satt',5e-6),('p_value','p_Satt',5e-8),('CI95_low','CI_L',5e-8),('CI95_high','CI_U',5e-8)]:
 rep['difference_'+pcol]=rep[rcol]-rep[pcol];checks.append(dict(quantity=pcol,terms=len(rep),max_absolute_difference=rep['difference_'+pcol].abs().max(),absolute_tolerance=tol,passed=bool(rep['difference_'+pcol].abs().le(tol).all())))
jr=joint.merge(read('R_joint_tests'),on='model',validate='one_to_one')
for pcol,rcol,tol in [('HTZ_F','Fstat',5e-8),('HTZ_df2','df_denom',5e-6),('HTZ_p','p_val',5e-8)]:
 jr['difference_'+pcol]=jr[rcol]-jr[pcol];checks.append(dict(quantity=pcol,terms=len(jr),max_absolute_difference=jr['difference_'+pcol].abs().max(),absolute_tolerance=tol,passed=bool(jr['difference_'+pcol].abs().le(tol).all())))
rm=read('R_metadata');mm=rm.merge(meta[['model','n','clusters','rank','residual_df']],on='model',suffixes=('_R','_Python'),validate='one_to_one')
for c in ['n','clusters','rank','residual_df']:assert mm[c+'_R'].eq(mm[c+'_Python']).all()
pm=read('sample_membership');rr=read('R_sample_membership')
for label in targets.keys()|{'A readiness interaction'}:
 a=pm[pm.model.eq(label)];b=rr[rr.model.eq(label)];assert set(zip(a.country_code,a.year))==set(zip(b.country_code,b.year))
rep.to_csv(EV/'independent_coefficient_comparison.csv',index=False);jr.to_csv(EV/'independent_joint_comparison.csv',index=False);checks=pd.DataFrame(checks);checks.to_csv(EV/'replication_discrepancy_summary.csv',index=False);assert checks.passed.all(),checks
mm.to_csv(EV/'independent_design_comparison.csv',index=False)
def table(df,cols):
 def f(v):
  if pd.isna(v):return 'not applicable'
  return f'{v:.6g}' if isinstance(v,(float,np.floating)) else str(v).replace('|','/')
 return '| '+' | '.join(cols)+' |\n|'+'|'.join(['---']*len(cols))+'|\n'+'\n'.join('| '+' | '.join(f(row[c]) for c in cols)+' |' for row in df.to_dict('records'))
loco_summary=[]
for label,term in list(targets.items())+[('A readiness interaction','RQ_x_RISE')]:
 z=loco[loco.model.eq(label)&loco.term.eq(term)];largest=z.loc[z.change.abs().idxmax()];loco_summary.append(dict(model=label,term=term,deletions=len(z),minimum_estimate=z.estimate.min(),maximum_estimate=z.estimate.max(),largest_change_country=largest.omitted_country,largest_absolute_change=abs(largest.change),parent_estimate=z.full_sample_estimate.iloc[0]))
ls=pd.DataFrame(loco_summary);ls.to_csv(EV/'LOCO_summary.csv',index=False)
sens=coef[coef.term.isin(['RISE','FDI','Deficit','RQ_x_RISE'])].copy();sens.to_csv(EV/'mechanism_sensitivities.csv',index=False)
# Coefficient-specific scales: none of these outcomes is an ordinary percentage change in commitments.
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':200})
def finish(fig,name):
 for ax in fig.axes:ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
 if name=='matched_depth_governance_comparison':
  handles,labels=fig.axes[0].get_legend_handles_labels()
  for ax in fig.axes:
   if ax.get_legend() is not None:ax.get_legend().remove()
  fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,1.08),ncol=2,frameon=False)
 fig.tight_layout();fig.savefig(FIG/(name+'.png'),bbox_inches='tight');fig.savefig(FIG/(name+'.pdf'),bbox_inches='tight');plt.close(fig)
fig,axs=plt.subplots(1,3,figsize=(12,3.8))
for ax,(label,term) in zip(axs,targets.items()):
 z=priority[priority.model.eq(label)].iloc[0];ax.errorbar(z.estimate,0,xerr=[[z.estimate-z.CI95_low],[z.CI95_high-z.estimate]],fmt='o',color='#245577',capsize=5);ax.axvline(0,color='grey',lw=1);ax.set_yticks([]);ax.set_title(label+' — '+term);ax.set_xlabel('Coefficient; CR2 95% interval');ax.text(.03,.1,f'n={z.n:,}; clusters={z.clusters}\nCR2 df={z.Satterthwaite_df:.1f}; Holm p={z.p_Holm_primary3:.2f}',transform=ax.transAxes)
fig.suptitle('Primary exploratory mechanism estimates',y=1.03);finish(fig,'primary_mechanism_intervals')
fig,axs=plt.subplots(1,2,figsize=(11,4))
for ax,term in zip(axs,['FD','Governance']):
 for i,family in enumerate(['A','B','C']):
  for shift,model,color in [(-.09,family+' matched baseline','#888888'),(.09,family+' primary','#245577')]:
   z=coef[coef.model.eq(model)&coef.term.eq(term)].iloc[0];ax.errorbar(z.estimate,i+shift,xerr=[[z.estimate-z.CI95_low],[z.CI95_high-z.estimate]],fmt='o',capsize=3,color=color,label='Matched baseline' if i==0 and shift<0 else 'Extension' if i==0 else None)
 ax.axvline(0,color='grey',lw=1);ax.set_yticks(range(3),['A support','B support','C support']);ax.set_title(term);ax.set_xlabel('Coefficient; CR2 95% interval');ax.legend(frameon=False)
finish(fig,'matched_depth_governance_comparison')
fig,axs=plt.subplots(1,3,figsize=(13,4.5))
for ax,(label,term) in zip(axs,targets.items()):
 family=label[0]
 for i,(m,short) in enumerate([(family+' primary','Main'),(family+' asinh inflation','asinh inflation'),(family+' all-X lag','All X at t−1')]):
  z=coef[coef.model.eq(m)&coef.term.eq(term)].iloc[0];ax.errorbar(z.estimate,i,xerr=[[z.estimate-z.CI95_low],[z.CI95_high-z.estimate]],fmt='o',capsize=3,color='#245577')
 ax.axvline(0,color='grey',lw=1);ax.set_yticks(range(3),['Main','asinh inflation','All X at t−1']);ax.set_title(family+' — '+term);ax.set_xlabel('Same target units; CR2 95% interval')
finish(fig,'prespecified_timing_inflation_sensitivities')
fig,axs=plt.subplots(1,3,figsize=(13,4))
for ax,(label,term) in zip(axs,targets.items()):
 z=loco[loco.model.eq(label)&loco.term.eq(term)].sort_values('estimate').reset_index(drop=True);ax.scatter(z.estimate,np.arange(len(z)),s=11,alpha=.55,color='#245577');ax.axvline(z.full_sample_estimate.iloc[0],color='#c17b38',label='Parent estimate');ax.axvline(0,color='grey',lw=.8);ax.set_title(label+' — '+term);ax.set_yticks([]);ax.set_xlabel('Leave-one-country-out estimate')
 ax.text(.02,.92,'Largest changes: '+', '.join(z.loc[z.change.abs().nlargest(3).index,'omitted_country']),transform=ax.transAxes,fontsize=8)
fig.text(.5,-.025,'Orange: full-sample estimate; grey: zero. Each dot omits one country.',ha='center',fontsize=9)
axs[0].set_ylabel('Country deletions, ordered by estimate');finish(fig,'country_deletion_sensitivity')
rep_report=f'''# Independent second-stage replication

Result: **PASSED — no material discrepancies.** R `stats::lm` and `clubSandwich` independently reproduced A, B and C primary models and the secondary A interaction. R reads the frozen candidate grid, independently reconstructs the 1,104-row benchmark, applies the new complete-case checks and iterated singleton removal, builds country/year effects and centered interaction itself. It does not import Python model matrices or coefficients.

CR2 uses identity working covariance, country clusters, coefficient-specific Satterthwaite tests/95% intervals, and HTZ tests. All **30 non-FE coefficients**, all four priority/secondary joint tests, exact country-year keys, ranks and residual df agree within recorded tolerances.

{table(checks,['quantity','terms','max_absolute_difference','absolute_tolerance','passed'])}

R session/package information is retained in `audit_evidence/R_session_info.txt`; local R 4.6.1 and clubSandwich 0.7.0 are reused without changing the original replication. R ran successfully through the Python stdout-capture launcher, which avoids this desktop environment's earlier R file-export visibility issue. Captured stdout, stderr and process status are preserved. Any locale warnings in stderr did not prevent successful fits.

Scope: this independent check covers the three new priority OLS models and secondary A interaction, not all 22 OLS sensitivities, the custom restricted bootstrap, or nonlinear PPML inference. Python additionally checked every OLS fit against a separate full-dummy statsmodels fit and the CR2 working-covariance identity. Independent numerical agreement verifies implementation under the chosen assumptions; it does not establish causal identification or completeness of the financing outcome.

Reproduction: run `scripts/stage2_run_R.py`, then `scripts/stage2_outputs.py`. The R script is `scripts/stage2_independent.R`; it uses the previously installed project-local library. The original 60-model estimates, bootstrap and independent replication are preserved. See [clubSandwich CR2 documentation](https://jepusto.github.io/clubSandwich/reference/vcovCR.html).
'''
(OUT/'independent_replication_second_stage.md').write_text(rep_report,encoding='utf-8')
primarycols=['model','term','estimate','SE_CR2','Satterthwaite_df','CI95_low','CI95_high','p_value','p_Holm_primary3','p_wild_restricted','p_wild_Holm_primary3','n','clusters']
fitrows=[]
for model in targets:
 z=comparisons[comparisons.extension.eq(model)].iloc[0];fitrows.append(dict(model=model,baseline_within_R2=z.baseline_within_R2,extension_within_R2=z.extension_within_R2,change=z.within_R2_change,partial_R2=z.partial_R2))
fitrows=pd.DataFrame(fitrows);fitrows.to_csv(EV/'fit_comparison_summary.csv',index=False)
rq=coef[coef.model.isin(targets)&coef.term.eq('Governance')]
report=f'''# Second-stage econometric results

## Research decision

**All three primary mechanism tests remain inconclusive.** Keep the original 2014–2023 model primary, report the full exploratory second stage transparently, and do not select replacement specifications for significance. This stage does not confirm the previously inconclusive FD, FD×RQ or income-heterogeneity hypotheses. No fsQCA membership, calibration or minimization was run.

The authorized execution completed **22 OLS specifications**, **four restricted-null wild-cluster tests with 9,999 valid draws each**, all-country deletion coefficient diagnostics for four parent models, and **six converged PPML specifications**. No planned fit failed. Independent R replication passed for three priority models plus the secondary interaction. All files are in the new `second_stage_econometric_results` directory; no earlier result is overwritten.

## Frozen data and source clearance

The outcome is `log1p(irena2025_pdf)`, with untransformed commitments measured in million constant-2022 USD. Maintain frozen IRENA 2025 recipient commitments, 2014–2023, bank private credit/GDP, Regulatory Quality, log real GDP/person, log population, trade/GDP, inflation, country and calendar-year FE. A adds lagged RISE; B adds lagged net FDI/GDP and current real GDP growth; C adds lagged access deficit, 100 minus access percentage. Missing commitments are never zeros. These are selected public-finance commitments, not disbursements or total public/private renewable investment. Published zero amounts may be rounded. The unresolved financing-release differences and high-income source selection from earlier reports remain.

The [official RISE methodology](https://rise.esmap.org/methodology) explicitly describes reconstructing past scores under revised indicators. Its renewable pillar uses five sectoral indicators. All **1,960 country-year pillar records** in the downloaded series equal their five-component arithmetic mean, maximum discrepancy 2.84e−14. Annual retrospective 2010–2023 scores are therefore consistent with the documented framework; they are not independent annual surveys or real-time policy information. Legacy Data360 descriptions (three pillars, 25/75 bands) conflict with the current framework (four pillars, 33/67 bands). Those descriptions are retained as metadata discrepancies; no traffic-light thresholds are used in estimation or calibration. The data gate supports provisional exploratory interpretation, not a causal-policy design.

## Primary results and inference

{table(priority,primarycols)}

Each displayed 95% interval is an unadjusted coefficient-specific CR2/Satterthwaite interval; adjusted p-values apply Holm across exactly the three declared mechanism targets. The bootstrap p-values have a separate three-test Holm family. They use country-level Rademacher weights, seed 20261016, restricted OLS null and re-estimated CR2 Wald studentization at every draw, with the plus-one p-value correction. They are not independently bootstrapped in R. Main controls' p-values are not part of the three-mechanism Holm family.

A does not provide evidence that greater formal policy readiness increases these selected public commitments. A near-zero estimate is compatible with countervailing bankability and development-finance targeting channels; it is not evidence that energy policies are irrelevant. B does not distinguish complementarity from substitution: its negative conditional point estimate is imprecise, with only about six Satterthwaite df and seven effective design clusters. C does not establish that access need predicts commitments conditional on the benchmark. High-income access is near universal, and high-income rows contribute only about 3% of deficit-slope design information; do not generalize a high-income need slope.

Coefficients are per one point of the policy score, one GDP-percentage-point of lagged FDI, or one population-percentage-point access deficit. They act on log(1 + commitments in millions), not directly on percentage changes in financing. Do not compare coefficient magnitudes across these distinct units as mechanism strength.

## Matched-support fit and coefficient comparisons

A uses **914 rows/109 clusters**, B **1,094/135**, C **1,094/136**. High-income support is respectively **149/26**, **178/34**, **178/34**. All four historical groups remain represented. Baselines are independently refitted on exactly the same retained keys as each extension. Changes in FD/RQ and their intervals are in `matched_baseline_comparisons.xlsx`; differences from the original 1,104-row baseline must not be attributed to a mechanism without matching support.

{table(fitrows,list(fitrows.columns))}

Within R² here is explanatory fit after removing both country and year effects; it is not out-of-sample predictive performance. Partial R² is the reduction in residual SSE relative to the exact matched baseline. B's fit increment belongs jointly to FDI and growth, not FDI alone. The small increments in A and C provide little additional within-panel explanatory fit. None is a causal mediation decomposition.

## Regulatory Quality and financial depth

{table(rq,['model','estimate','SE_CR2','Satterthwaite_df','CI95_low','CI95_high','p_value'])}

RQ remains positively associated in the primary OLS extensions. The retained coefficient is a conditional association in the observed selected-public-finance sample, not a causal governance effect. Financial depth remains imprecise in the primary OLS extensions. The six previously estimated GE/RL/CC models were not rerun; the three earlier baseline result sets are referenced directly in the robustness workbook. VA and Political Stability were separately examined as predeclared comparisons of participatory legitimacy and political-risk conditions, not as significance-selected replacements for RQ.

## Interaction, sensitivity and concentration

The secondary centered RQ×lagged RISE estimate is **−0.011902**, CR2 95% interval **[−0.028110, 0.004307]**, Satterthwaite df **10.17**, p **0.1331**; restricted bootstrap p **0.0901**. Both main effects are retained, with sample grand-mean centering. This weak/inference-sensitive result does not establish complementarity, substitution or validate the earlier FD×RQ hypothesis; the interaction is between different constructs. Sweden supplies about 23.6% of its design information.

Pre-frozen inflation sensitivity replaces inflation percentage points with `asinh(inflation / 1 percentage point)` in all observations and refits its matching baseline. No extreme-inflation country is deleted. The three mechanism conclusions remain inconclusive. All-X timing alternatives lag baseline controls one year while keeping already-lagged mechanisms at t−1; growth is also lagged in B. Their samples are A **910/109**, B **1,088/135**, C **1,089/136**, with matched lagged baselines. No lag skips a calendar gap. These alternatives remain inconclusive.

Signed FDI can legitimately be negative or extreme. The top 20 actual lagged FDI rows are disclosed, not deleted. B's asinh signed-FDI sensitivity retains all 1,094 observations and 135 clusters. Its slope is −0.012507, p 0.8773, in transformed units; it must not be compared one-for-one with the raw percentage-point slope. The transformation increases effective design clusters to about 27, but does not create evidence of a mechanism. Raw FDI information is concentrated: Netherlands 28.9%, top five 68.5%. Predeclared top contributors and all countries are retained in parent fits.

{table(ls,list(ls.columns))}

Country-deletion estimates are diagnostics on heavily overlapping samples, not independent replications or automatic exclusion rules. A and C point estimates cross zero across deletions. B remains negative across individual deletions, which does not make its wide parent interval significant. The secondary interaction remains negative in these deletions but still has imprecise cluster inference. No observation was removed from a primary model to strengthen a result. Original inflation is dominated by extreme-inflation country information; the asinh sensitivity is reported in full rather than used to select a preferred estimate.

## PPML robustness and limits

All six extension/matched-baseline PPML fits passed the positive-outcome design-rank separation screen and converged, with scaled score checks saved. Country/year FE use nonnegative observed commitment levels; there are no artificial zeros. Mechanism PPML intervals all include zero. RQ intervals also include zero in these level-mean models, unlike the OLS log1p models: governance evidence is therefore conditional on functional form and weighting, not universally robust across estimators.

PPML uses **country-cluster CR1 with t(G−1)**, explicitly not nonlinear CR2/Satterthwaite. Nonlinear CR2 and PPML wild bootstrap were **not implemented**. The four completed wild-bootstrap tests apply to OLS only. PPML targets a different conditional mean and weights large commitment amounts differently. It cannot remedy selective source coverage, missing reporting, endogenous allocations or unresolved vintage differences. A fit passing the separation screen is not a guarantee of causal identification or adequate group-specific inference.

## Verification, reproducibility and interpretation

Independent R/clubSandwich replication passed all 30 identified non-FE coefficients and four HTZ tests, including SEs, Satterthwaite df, intervals, exact support and ranks. `independent_replication_second_stage.md` gives numerical discrepancies and scope. Python additionally checked every OLS fit against full-dummy statsmodels estimates and the CR2 identity. Execution scripts, the frozen pre-estimation protocol, source hashes, R session, all model keys, covariance matrices, fit/inference diagnostics and plotted data are retained. The original 60 OLS specifications, their bootstrap results, independent replication, 111-country fsQCA preparation and unapproved anchors are preserved.

Country/year FE and lagging do not remove endogenous lending allocation, reverse policy responses, project pipelines, donor eligibility rules or correlated shocks across countries. Cluster df and information concentration, particularly for FDI and high-income contrasts, matter more than nominal country counts alone. No new income heterogeneity is established or re-estimated here. Results apply to observed recipient-country commitments within the institutional perimeter, not the complete global renewable-finance market.

**Final recommendation:** present the three unsuccessful primary mechanism tests and predeclared sensitivities as exploratory additional evidence, retain the original primary specification, qualify RQ as an OLS conditional relationship with functional-form sensitivity, and make no new claim of causal policy effects, FDI complementarity/substitution or robust income heterogeneity. Do not seek replacement mechanisms by statistical significance. Any separate fsQCA work remains behind the existing protocol-review boundary and is not validation of these inconclusive panel results.
'''
(OUT/'second_stage_econometric_report.md').write_text(report,encoding='utf-8')
def sheet(name,df,note):
 df=df.replace([np.inf,-np.inf],np.nan).astype(object);df=df.where(pd.notna(df),None)
 return dict(name=name,note=note,columns=list(df.columns),rows=df.to_dict(orient='records'))
primary_models=list(targets);rawdict=pd.DataFrame([dict(variable=v,definition=de,unit=u) for v,de,u in [('outcome','log1p(frozen IRENA2025 selected public-finance commitments)','Input million constant2022USD'),('FD','Domestic bank private credit/GDP','% GDP'),('Governance','WGI Regulatory Quality latent estimate, centered','approx -2.5 to2.5'),('RISE','Official retrospective renewable pillar, exact t−1, centered','0–100 score points'),('FDI','Net inward FDI/GDP, exact t−1, signed','% GDP'),('growth','Frozen contemporaneous real GDP growth, B only','annual %'),('Deficit','100 minus publisher electricity access%, exact t−1','population percentage points'),('asinh sensitivities','Researcher transforms explicitly labelled, no imputation','asinh per1 percentage point'),('within_R2','Two-way-FE residualized fit; not prediction','fraction')]])
payload={'second_stage_regression_results.xlsx':[sheet('Primary mechanism summary',priority,'All primary intervals include zero. Holm family has exactly three mechanism targets; outcomes and estimates unchanged by selection.'),sheet('Primary coefficients',coef[coef.model.isin(primary_models)],'Country-cluster CR2; coefficient-specific Satterthwaite intervals and df.'),sheet('Secondary policy interaction',coef[coef.model.eq('A readiness interaction')],'Centered RQ × lagged RISE, with both main effects; secondary exploratory hypothesis.'),sheet('OLS model metadata',meta,'All 22 planned OLS fits; full FE rank; residual df is not cluster-test df.'),sheet('Historical group support',support,'All four historical income groups; country counts can overlap across groups over time.'),sheet('Variable dictionary',rawdict,'Definitions remain distinct; capacity or generation never substitutes for financing.')], 'matched_baseline_comparisons.xlsx':[sheet('Primary fit comparison',fitrows,'Fit increments use identical retained keys; B increment belongs to both FDI and growth.'),sheet('Coefficient and interval shifts',comparisons,'Matched baseline refitted independently; no causal mediation interpretation.'),sheet('Matched baseline coefficients',coef[coef.model.str.contains('baseline')],'Includes corresponding transformed-inflation and all-X-lag baseline definitions.'),sheet('Primary sample keys',pm[pm.model.isin(primary_models+[x[0]+' matched baseline' for x in primary_models])],'Exact priority keys; all22 model keys retained in audit_evidence/sample_membership.csv.')], 'second_stage_robustness_and_inference.xlsx':[sheet('Restricted wild bootstrap',boot,'9,999 valid draws, seed20261016; OLS null imposed; re-studentized CR2 Wald; separate Holm primary3 family.'),sheet('Prespecified sensitivities',sens,'Target units differ for signed-asinh FDI; all variants disclosed.'),sheet('Country deletion summary',ls,'All-country deletion diagnostics; fixed parent centering; no deletion changes a parent model.'),sheet('Country deletion coefficients',loco,'Each country omitted once from four parent designs; overlapping diagnostic samples, not independent evidence.'),sheet('PPML results',ppml,'Levels, country/year FE; country CR1 only. No nonlinear CR2 or PPML bootstrap.'),sheet('HTZ joint tests',joint,'Four one-restriction tests; same target nulls as declared priority/secondary hypotheses.'),sheet('Original governance reused',read('reused_original_governance_results'),'Original GE/RL/CC baseline rows reused unchanged; all earlier specifications remain in their archive.'),sheet('VA and stability definitions',coef[coef.model.str.contains('benchmark')],'Conceptually distinct participation/legitimacy and political-risk comparisons, not p-value-selected controls.'),sheet('Largest lagged FDI rows',read('largest_FDI_observations'),'Actual signed observations; no winsorization, exclusion or zero replacement.'),sheet('Independent R discrepancies',checks,'All30 non-FE coefficients and four HTZ tests; exact sample/rank checks passed.'),sheet('Implementation checks',read('implementation_checks'),'Separate full-dummy OLS estimates/residuals and CR2 working-covariance identity.'),sheet('Execution status',read('execution_status'),'Actual failures, if any; nonlinear CR2/PPML bootstrap intentionally unavailable.') ]}
(EV/'workbook_payload.json').write_text(json.dumps(payload,ensure_ascii=False,allow_nan=False),encoding='utf-8')
for src in Path(__file__).parent.glob('stage2_*'):
 if src.suffix not in ['.py','.R']:continue
 dst=SCRIPTS/src.name
 if src.resolve()!=dst.resolve():shutil.copy2(src,dst)
(OUT/'execution_instructions.md').write_text('''# Reproduction\n\nUse Python 3.13 with numpy, scipy, pandas, statsmodels, matplotlib and threadpoolctl. From scripts run stage2_setup.py, stage2_estimate.py, stage2_run_R.py, stage2_outputs.py, node build_second_stage_results.mjs, then stage2_verify.py. The R launcher reuses project-local R4.6.1 and clubSandwich0.7.0 and captures stdout to avoid desktop file-export visibility issues. Existing R outputs are cached; remove only the new R_stdout.txt if deliberately repeating that check.\n\nThe original econometric engine was copied unchanged into scripts/engine_frozen_copy.py. The new adapter constructs the authorized designs while leaving the earlier engine/results untouched. The frozen pre-estimation protocol is audit_evidence/estimation_protocol_frozen.json. Wild bootstrap: 9999 draws, seed20261016, Rademacher country weights, restricted null, CR2 studentization. No fsQCA code is called.\n\nThe workbook builder uses @oai/artifact-tool from the bundled runtime. Exports are static audit artifacts: rerun the scripts to refresh, rather than editing numeric cells. All model keys/covariances, independent comparisons, original-governance references, raw FDI diagnostics and plotted data are retained. Source paths are explicit for this project. Permission is required where the research folder is outside the shell workspace.\n''',encoding='utf-8')
(EV/'replication_validation.json').write_text(json.dumps(dict(status='PASSED',terms=30,joint_tests=4,maxima=checks.set_index('quantity').max_absolute_difference.to_dict(),sample_keys_identical=True,rank_df_identical=True),indent=2),encoding='utf-8')
print('REPORTS FIGURES PAYLOAD COMPLETE',flush=True);print(checks.to_string(index=False),flush=True)
