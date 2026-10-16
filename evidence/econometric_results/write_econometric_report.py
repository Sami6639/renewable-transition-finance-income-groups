from pathlib import Path
import json,pandas as pd,numpy as np,hashlib
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'econometric_results';EV=OUT/'audit_evidence';read=lambda n:pd.read_csv(EV/(n+'.csv'))
def md(df,cols):
 def fmt(x):return '' if pd.isna(x) else f'{x:.6g}' if isinstance(x,(float,np.floating)) else str(x)
 return '| '+' | '.join(cols)+' |\n| '+' | '.join(['---']*len(cols))+' |\n'+'\n'.join('| '+' | '.join(fmt(r[c]) for c in cols)+' |' for r in df.to_dict('records'))
c=read('coefficients');m=read('model_metadata');co=read('marginal_effects_and_contrasts');j=read('joint_tests');b=read('wild_bootstrap');l=read('LOCO_summary');pp=read('PPML_results');ppj=read('PPML_joint_tests');s=read('group_support');des=read('descriptive_statistics');ref=json.loads((EV/'independent_CR2_reference_validation.json').read_text());receipt=json.loads((EV/'run_receipt.json').read_text())
main=c[c.model.eq('Primary Baseline')];rq=c[c.model.eq('Primary RQ interaction')];inc=co[co.model.eq('Primary Income interaction')&co.evaluation.eq('Group FD slope')];base=c[c.block.eq('Baseline')&c.term.eq('FD')];top=read('largest_published_commitments').iloc[0];mdrows=m[m.model.str.startswith('Primary')];groups=s[s.model.eq('Primary Baseline')];d=pd.read_csv(ROOT/'reports/final_model_selection_20261016/four_income_group_final_dataset.csv');a=d[d.sample_B_Bank_credit];dsc=des[des.scope.eq('Primary FE sample')&des.income_group.eq('Pooled')];selcov=m[m.model.isin(['FinData native Baseline','IMF native Baseline'])]
text=f'''# Exploratory econometric results

**Financing the Renewable Energy Transition: Financial Development, Institutional Quality, and Income-Group Heterogeneity**

The approved 2014–2023 IRENA 2025 analysis has been executed. Sixty prespecified OLS specifications, five restricted-null wild cluster bootstraps with 9,999 draws each, all-country deletion checks for the three primary blocks, and three converged PPML specifications are reported. No fsQCA was estimated. Primary models and legitimate extremes were retained regardless of signs or p-values.

The main results are **inconclusive about a contemporaneous financial-depth association, financial–institutional complementarity or systematic four-group slope heterogeneity**. Some lagged and PPML sensitivities show signals; these are disclosed below and are not promoted into the primary specification. Inconclusive tests do not prove zero effects or group equality. The evidence concerns conditional associations in selected reporting institutions' published public-finance commitments to recipients, not total renewable investment, private finance or causal effects.

## 1. Frozen data, scope and descriptive evidence

Source: `reports/final_model_selection_20261016/four_income_group_final_dataset.csv`, SHA256 `{receipt['source_sha256']}`. All 2,170 grid rows remain unchanged. No IRENA 2026 observations or 2024 outcome were appended. The outcome is the frozen PDF-derived published renewable total, million constant 2022 USD. Its earlier API/provider-year and transaction-revision limitations remain; the prior provisional methodological clearance permits exploration, not a complete national financing census.

The six slopes are bank private credit/GDP, Regulatory Quality, log real GDP per capita (constant 2015 USD), log population, trade/GDP and CPI inflation. Country and calendar-year effects enter every primary model. All macrofinancial values are official source values; no assistant interpolation, missing-to-zero or index splicing occurs. FinData publisher imputations are flagged rather than redescribed as observed inputs.

There are 1,114 complete rows before FE singleton removal and **1,104 estimation rows/137 country clusters/10 years** afterwards. Fourteen published numeric financing zeros are retained; all unavailable financing stays missing. A published 0.00 may reflect rounding and is not established as a structural absence of finance. The largest published record is {top.country_name}, {int(top.year)}: {top.irena2025_pdf:,.2f} million constant 2022 USD, retained when sample-eligible. The distribution figure documents the heavy tail; no winsorization is applied.

{md(groups,['income_group','rows','clusters','years'])}

Country counts across historical income groups overlap because countries switch groups; they cannot be summed as independent clusters. High income has 178 pooled rows across 34 clusters, with selective recipient coverage. The separate-group FE audit's 175/31 is a different support definition; no separate-group OLS regressions are substituted for the approved pooled interactions.

{md(dsc,['variable','nonmissing','mean','SD','median','maximum','skewness','numeric_zeros'])}

The initial dependent variable is **ln(1 + financing measured in million USD)**. The additive one is unit dependent: changing to USD or billion USD changes the model, especially around zero and small commitments. These unit sensitivities are disclosed without changing the prespecified primary transformation. A coefficient times ten is a change in conditional log-plus-one financing per ten percentage points of bank credit/GDP. It is not automatically a percentage change in mean untransformed financing; retransformation and Jensen's inequality preclude that shortcut. Positive-only log and level PPML are separate estimands.

## 2. Estimation and inference

OLS uses the Frisch–Waugh–Lovell projection against the full country/year FE span, with rank-aware QR and slope scaling for numerical stability. FE and identified main effects are retained. FD and governance are grand-mean centered **before** their product; the common baseline centering does not change slopes. Interaction blocks are separate in primary models; the combined model is an explicitly labelled sensitivity.

Generalized **CR2** uses an identity working covariance and country clusters. With full-model hat matrix H and cluster selector E_g, the adjustment is the symmetric inverse square root of (I−H)_gg on its positive residual subspace. Null directions introduced by country effects are handled by the Moore–Penrose inverse. Each contrast's residual-linear score vectors provide its Satterthwaite moment degrees of freedom. The full FE span is included in leverage adjustment; ordinary standard errors or naive G−1 degrees of freedom were not substituted for OLS. The multi-constraint joint test uses the moment-matched **HTZ** correction. Confidence intervals are pointwise 95%; the six income pairwise p-values additionally receive Holm correction, not simultaneous CI coverage.

Implementation checks: every OLS fit reproduces an independent full FE-dummy `statsmodels.OLS` fit; the CR2 working-model covariance identity is verified on every real design; the singleton-cluster limit matches independent HC2. A separate, authentic clubSandwich **MortalityRates** benchmark reproduces published unweighted two-way-FE estimates, CR2 SEs and Satterthwaite df: legal-drinking-age coefficient 7.587708, SE 2.513082, df 24.578519; beer-tax coefficient 3.818671, SE 5.265016, df 5.768415. Those benchmark observations are used solely for method validation and never merged into the research panel. R is unavailable; the empirical estimates use the documented Python implementation, not an unreported R run.

For the five declared interaction hypotheses, unrestricted-model CR2 studentization is recomputed in every **restricted-null wild cluster bootstrap** draw. Constrained OLS fitted values and restricted residuals generate clusterwise Rademacher sign flips; 9,999 draws/hypothesis, seed **20261016**, two-sided squared-Wald statistics, plus-one finite-draw p-values. The original model matrix and samples are held fixed. Bootstrap Monte Carlo SE and a supplementary Holm adjustment across all five declared bootstrap tests are retained. No simulated empirical observations or fabricated findings are introduced.

References: [CR2 and working models](https://jepusto.github.io/clubSandwich/reference/vcovCR.html), [panel FE reference example](https://jepusto.github.io/clubSandwich/articles/panel-data-CRVE.html), [joint HTZ implementation](https://github.com/jepusto/clubSandwich/blob/master/R/Wald_test.R), [restricted wild-bootstrap methods](https://www.stata.com/manuals/rwildbootstrap.pdf). The code implements its labelled CR2-studentized variant and does not claim byte-for-byte equivalence to a particular software default.

## 3. Primary baseline

{md(main,['term','estimate','SE_CR2','Satterthwaite_df','CI95_low','CI95_high','p_value'])}

Bank-depth coefficient **−0.010955**, CR2 SE **0.007017**, Satterthwaite df **14.23**, 95% CI **[−0.025983, 0.004073]**, p **0.1404**. The point estimate is negative, but the data do not distinguish it precisely from zero under the specified inference. This is not evidence of a positive financing channel and does not establish crowding out. Regulatory Quality is also reported with its own uncertainty; no control is omitted because of insignificance.

Regulatory Quality has a positive conditional association in the baseline: coefficient **2.302846**, CR2 SE **0.577764**, Satterthwaite df **39.11**, 95% CI **[1.134312, 3.471381]**, p **0.000285**. This institutional main-effect result must be reported alongside the inconclusive depth and complementarity results. It does not establish that institutional reform causes higher financing or that the financial-depth slope changes with governance.

The baseline within R² is **0.03417**, defined relative to the country-and-year-FE-only residual outcome; overall R² is **0.58988**. The latter largely includes country/year differences and is not evidence that the six slopes explain most financing variation. Primary baseline rank 152 leaves 952 algebraic residual df; these are distinct from the 14.23 inferential df for bank depth.

## 4. Financial–institutional complementarity

{md(rq,['term','estimate','SE_CR2','Satterthwaite_df','CI95_low','CI95_high','p_value'])}

FD×RQ is **−0.006519**, SE **0.009298**, df **17.95**, CI **[−0.026057, 0.013018]**, Satterthwaite p **0.4922**; restricted-bootstrap p **0.4801** (Monte Carlo SE about 0.0050). This is inconclusive, not support for a positive complementarity claim.

The financial-depth marginal slope is β_FD + β_interaction × (RQ−mean RQ). Observed RQ in this sample spans **{a.regulatory_quality_estimate.min():.6g} to {a.regulatory_quality_estimate.max():.6g}**; center **{a.regulatory_quality_estimate.mean():.6g}**. Tables show observed quantiles and extremes; a 61-point curve and pointwise contrast-specific CR2/Satterthwaite intervals stay within that range. Tail estimates are supported by fewer actual cases, visible in the plot's rug. No causal marginal effect or out-of-support extrapolation is asserted.

## 5. Historical income-group heterogeneity

Group main effects and three FD interactions accompany main terms, with low income as reference. Country effects absorb time-invariant group effects; historical group main effects are identified through actual switchers. Main group coefficients concern differences at the FD centering value and are not slopes.

{md(inc,['income_group','estimate','SE_CR2','Satterthwaite_df','CI95_low','CI95_high','p_value','effective_design_clusters','largest_information_country','largest_information_share'])}

The three interactions jointly yield **HTZ F(3, 9.679)=0.2955, p=0.8279**; restricted-bootstrap p **0.7592**. All six primary pairwise Holm-adjusted p-values are **1.0**. All four pointwise own-slope intervals include zero. This is a weakly informative heterogeneity result, not proof that slopes are equal.

The high-income own slope has only **8.80 Satterthwaite df** and approximately **9.67 effective design clusters**, despite 34 countries represented in high-income rows. Greece contributes about 28% of its design information; top five about 52%. The lower-middle-minus-low contrast has about 3.74 effective design clusters, about half its design information from Nepal, and Satterthwaite df only **3.25**. These design concentration measures differ from realized influential-country diagnostics and from nominal cluster counts. They reinforce the prior warning that four-group comparisons are exploratory.

Fixed-2014 income recoding and the stable-income sample are both estimated and reported, with freshly recalculated rank and support. Stable income means no classification changes across the entire 2014–2023 grid, not merely observed financing years. Fixed-baseline models absorb group main effects; the stable restriction has 765 rows/96 clusters and changes the target population. The combined interaction model is reported separately; its FD×RQ curve is the **low-income-reference** marginal slope, and its group slopes are evaluated at mean RQ. These conditional quantities must not be confused with the common-slope RQ model.

## 6. Robustness, influence and alternative estimands

All planned baseline variants:

{md(base,['model','estimate','SE_CR2','Satterthwaite_df','CI95_low','CI95_high','p_value'])}

Exact all-X lags use authentic prior-calendar-year values, including 2013 inputs for 2014 outcomes. The lagged native baseline has 1,114 rows/138 clusters, β_FD **−0.015276**, p **0.0360**. On the exact current/lag intersection (1,099 rows/137 clusters), current β **−0.011284**, p **0.1304**, versus lagged β **−0.015050**, p **0.0409**. This negative lagged association is a disclosed sensitivity signal; the many tested sensitivities, temporal confounding and outcome selection prevent turning it into a new confirmatory or causal result. No lead-lag model was selected by significance.

Government Effectiveness, Rule of Law and Control of Corruption replace RQ **one at a time** in baseline and FD×governance models. Their full coefficients, marginal slopes and tests remain in evidence. No preferred governance indicator is chosen after results. The sovereign-scope restriction follows the existing flag, removing Kosovo's ten observations and leaving 1,094/136; it is a scope sensitivity, not a new legal judgment about sovereignty. The primary sample already contains no dependent-territory rows; observer-state rows remain in the prior scope flag.

Positive-only log uses 1,090 rows/137 clusters and answers an intensive-margin question conditional on a published strictly positive amount. Its baseline depth estimate is −0.017534, p 0.0889; it does not replace log-plus-one. The USD and billion-USD log-plus-one sensitivities demonstrate scale dependence and are both disclosed. They do not establish a preferred monetary unit or independent robustness of the zero-handling convention.

{md(l,['model','term','full_estimate','minimum_leave_one_out','maximum_leave_one_out','sign_reversals','largest_change_country','countries_deleted'])}

Each of 137 primary countries is omitted separately in each of the three primary blocks. Parent-sample FD/RQ centers are fixed for comparable deletion coefficients. The common depth estimate ranges from −0.013725 to −0.008220 and the FD×RQ estimate from −0.009807 to −0.001518; no sign reversal does not make either coefficient statistically established. The low-income reference slope changes sign under 12 deletions; income differences can approach zero after omitting Angola. Full deletion tables are retained. Refitted deletion CIs were not computed, so no claim is made about deletion-specific significance.

Observation leverage, residuals, classical Cook distance and CR2-SE-scaled DFBETAs are supplied. Cook distance and scaled DFBETAs are descriptive influence diagnostics, not clustered p-values or grounds for automatic removal. All legitimate large commitments, including Angola, Brazil and the United Kingdom, remain in primary estimates. The prior UAE-2022 release discrepancy remains unresolved; its 2025 value is retained, not replaced by a 2026 value.

### PPML

The three level-outcome PPML models converge in 11 IRLS iterations each. The strictly-positive-outcome design has the full model column rank in each specification, which rules out a nonzero linear separating direction after the FE-core screen; all-zero FE groups/singletons and nonnegative outcomes are checked. Scaled score errors are below 10⁻⁶. These certificates and convergence diagnostics are saved, rather than assuming that an optimizer's flag alone proves existence. No unsupported row is created to represent a national extensive margin.

For nonlinear PPML, generalized CR2/Satterthwaite and a nonlinear restricted wild bootstrap were **not implemented**. The explicitly labelled available inference is the country-cluster **CR1** sandwich with t(G−1) intervals and an approximate CR1 F test; it is not presented as small-sample CR2. This is a methodological limitation of these robustness models.

PPML common-depth β **−0.002248**, CR1 p **0.8228**; FD×RQ β **−0.001972**, CR1 p **0.8551**. The PPML income block gives an approximate CR1 joint p **0.0366**, unlike primary log-plus-one HTZ/bootstrap results; some raw interaction coefficients are also statistically distinguishable under CR1. All PPML group slopes, six Holm-adjusted pairwise contrasts and intervals are reported. These signals cannot be suppressed, but they do not validate a robust universal four-group result: the level conditional-mean estimand places different weight on large commitments and the nonlinear inference lacks the primary small-sample correction. An additional distribution-based diagnostic deletes each of the five countries with the largest single commitments (Angola, Brazil, India, the UK and Egypt), irrespective of coefficient significance, without dropping them from primary PPML. Fractional monetary outcomes are permitted under a Poisson pseudo-mean model; a count distribution is not asserted.

## 7. Alternative financial-development measures

Identical country-year supports were established **before estimation** and FE eligibility recalculated: bank–FinData 922/114 (2014–2023); FinData–IMF 631/108; bank–IMF 740/128; all three measures 613/106 (2014–2020). Each comparison estimates baseline, RQ and income-interaction blocks on identical keys. Broader native FinData 956/116 and IMF 763/131 are reported separately and cannot isolate measure changes from sample composition.

Bank credit is banking **depth**, FinData is a 0–100 overall depth/access/efficiency composite, IMF FD is a different 0–1 composite. Their raw coefficient magnitudes are not directly comparable. Do not splice them or treat a unit change as equal economic variation. Baseline native FinData p **0.8460** and IMF p **0.9423** are imprecise; exact-matched specifications are likewise fully disclosed. Sign changes across supports and measures are not evidence that one index measures the opposite causal channel.

**944 of 956** FinData native estimation rows have at least one known publisher-imputed input component. All 29 component provenance flags and raw/imputed/unresolved counts remain in the source dataset; missing unresolved components are distinct from filled components. No assistant-generated values are introduced. Blank WDI observation status likewise does not certify nationally unestimated measurements. No IMF values after 2020 were imputed; the preserved older 2015–2020 benchmark remains untouched.

## 8. Scientific interpretation and reporting decision

These are exploratory within-country associations conditional on observed selected-public commitments and controls. Country effects absorb stable differences, year effects common annual shocks; neither removes time-varying omitted factors, reverse causality, provider/recipient selection or project-timing/reporting revisions. Cross-country dependence may persist despite year effects; ten annual clusters do not support casual reliance on ordinary two-way clustering. A wild bootstrap does not repair missing institutional coverage, concentrated contrast information or endogenous finance demand.

Missing financing can reflect structural source scope, unavailable reporting or absence of qualifying records; it is never automatically zero. High-income observed recipients are a selected subset, with financing instruments/providers unlike a national census. External validity is limited to the included reporting recipients and this financial-depth construct. No missing-at-random assumption, selection correction or all-investment causal claim is justified.

Retain the original primary model and report the imprecise primary evidence, all sensitivity signals and inference differences. A Q1 submission should emphasize measurement limits, uncertainty, transparent sample selection and the narrow selected-public-finance estimand. The current results do **not** justify a headline that financial development causes renewable investment or that complementarity/income heterogeneity has been established. Do not change the period, transformation, group coding or covariates to pursue desirable signs or p-values. The separate fsQCA design remains unestimated pending review.

## Outputs and execution

Five requested workbooks, four PNG/vector-PDF figures, exact model membership, covariance matrices, bootstrap receipts, full deletion/influence tables, method-validation evidence, scripts and pinned requirements are in `econometric_results`. Workbook cells are static exports of the reproducible Python calculations, not independent spreadsheet estimators. See `execution_instructions.md`. No fitting procedure failed in this run; unavailable nonlinear small-sample inference, lack of deletion-specific CIs, broad primary uncertainty and limited external validity are explicitly recorded rather than silently substituted.
'''
(OUT/'econometric_results_report.md').write_text(text,encoding='utf-8')
instructions='''# Reproduce the exploratory panel results

Use Python 3.13.9 and the pinned packages in requirements.txt. Keep the existing project folder hierarchy and the frozen source CSV at reports/final_model_selection_20261016/four_income_group_final_dataset.csv. The scripts currently point explicitly to this project's absolute Windows path; update ROOT in scripts if relocating the project. Do not substitute a 2026 financing file.

From the econometric_results folder, use this sequence (an explicit Python executable may be needed on Windows):

```powershell
python -m pip install -r requirements.txt
python econometric_analysis.py
python prepare_econometric_outputs.py
python write_econometric_report.py
node build_econometric_workbooks.mjs
python verify_econometric_results.py
```

The numerical pipeline estimates the matched samples first, then native and sensitivity models, five 9,999-draw null-restricted Rademacher bootstraps (seed 20261016), all-country deletion diagnostics and three PPML fits. prepare_econometric_outputs.py creates extra PPML contrasts and figures. Runtime depends on BLAS; numerical linear algebra uses one thread for reproducibility. Logs print every planned OLS fit, with explicit procedure failures if any. Coefficients, model metadata, membership and covariances are in audit_evidence CSVs. Never treat a missing CSV financing value as zero.

Excel authoring uses the bundled Node @oai/artifact-tool runtime through the local node_modules junction. It is not an undeclared Python Excel writer. If that runtime is unavailable elsewhere, numerical results and figures can still be reproduced; use workbook_payload.json as the complete workbook data specification rather than silently changing estimator calculations. Preserve raw inputs and create a new results directory for a changed analysis. The builder exports five workbooks and renders sheet previews; unchanged previews are reused.

Optional independent CR2 benchmark validation requires rdata==1.1.0 and the retained MortalityRates.rda fixture. check_cr2_reference.py uses the task-local R-data-reader path; adjust that path if needed. It reproduces the authors' published real-data reference, not a simulation and not a replacement research dataset. requirements-validation.txt supplies optional validation dependencies. Final verification uses openpyxl read-only, never for workbook authoring. No fsQCA or 2024 financing belongs in this run.
'''
(OUT/'execution_instructions.md').write_text(instructions,encoding='utf-8')
(EV/'method_availability.json').write_text(json.dumps(dict(OLS_CR2_Satterthwaite='Implemented; independently validated against published real-data reference',joint_CR2_HTZ='Implemented; q=1 identity checked; full formulas documented in source',restricted_wild_bootstrap='5 declared hypotheses, 9999 valid draws each, CR2 studentization',PPML_CR2_Satterthwaite='UNAVAILABLE: nonlinear extension not implemented; explicitly labelled CR1 supplied',PPML_restricted_bootstrap='Not implemented; primary predeclared bootstrap hypotheses are OLS interactions',LOCO_refitted_CIs='Not calculated; coefficient stability only',failed_estimation_procedures=receipt['failed_procedures'],fsQCA='Not run'),indent=2),encoding='utf-8')
print('Scientific report and execution instructions saved')
