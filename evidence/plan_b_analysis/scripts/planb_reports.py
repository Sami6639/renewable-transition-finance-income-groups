from planb_setup import *
import pandas as pd,numpy as np,re
co=pd.read_csv(EV/'coefficients.csv');meta=pd.read_csv(EV/'R_model_metadata.csv').merge(pd.read_csv(EV/'model_plan.csv'),on=['model','n','clusters']);supp=pd.read_csv(EV/'income_support.csv');var=pd.read_csv(EV/'variation_decomposition.csv');con=pd.read_csv(EV/'concentration.csv');coverage=pd.read_csv(EV/'coverage_by_income_year.csv');loco=pd.read_csv(EV/'country_deletion_diagnostics.csv');ct=pd.read_csv(EV/'within_between_contrasts.csv');jt=pd.read_csv(EV/'R_joint_tests.csv');checks=pd.read_csv(EV/'independent_validation.csv');agg=pd.read_csv(EV/'R_independent_aggregation_checks.csv');plan=json.loads((EV/'preanalysis_plan_frozen.json').read_text());hist=pd.read_csv(EV/'country_timing.csv');design=pd.read_csv(EV/'design_metadata.csv')
def table(x,cols=None):return (x[cols] if cols else x).to_markdown(index=False,floatfmt='.5g')
primary=co[(co.model=='annual_CRE_absolute')&co.term.isin(['W_FD','B_FD','W_RQ','B_RQ'])];pcols=['model','term','estimate','SE_CR2','Satterthwaite_df','p_value','p_Holm_primary8','CI95_low','CI95_high'];annual=co[(co.model.isin(['annual_TWFE_absolute','annual_MundlakOLS_absolute','annual_TWFE_per_capita','balanced_TWFE','lagged_TWFE']))&co.term.isin(['FD','RQ','W_FD','W_RQ','B_FD','B_RQ'])]
multi=pd.read_csv(EV/'annual_window_comparisons.csv');gsummary=coverage.groupby('income_group')[['economies','observed_financing','missing_financing','primary_model_observations']].sum().reset_index();gsummary['missing_share']=gsummary.missing_financing/gsummary.economies
focusmodels=['annual_TWFE_absolute','annual_CRE_absolute','balanced_TWFE','lagged_TWFE','window2_TWFE_absolute','window2_CRE_absolute','matched2_annual_TWFE_absolute','window3_TWFE_absolute','matched3_annual_TWFE_absolute'];supportfocus=supp[supp.model.isin(['annual_TWFE_absolute','balanced_TWFE','window2_TWFE_absolute','window3_TWFE_absolute'])]
vs=var[(var['sample']=='Native complete cases')&var.variable.isin(['FD','RQ','Yabs','Ypc','log1p_absolute','log1p_per_capita'])];full=hist[hist.observed_years==10];ctshort=ct[['term','estimate','SE_CR2','Satterthwaite_df','p_value','p_Holm_two_differences','CI95_low','CI95_high']]
lsummary=loco.groupby(['model','term']).agg(min_deleted_estimate=('estimate','min'),max_deleted_estimate=('estimate','max'),max_absolute_change=('change',lambda x:abs(x).max())).reset_index()
report=f'''# Bounded exploratory Plan B: within–between variation and financing timing

## Executive assessment

Plan B does **not** establish a positive banking-depth mechanism that was hidden in the original TWFE results. Banking depth remains inconclusive both within and between selected recipient economies in the primary annual CRE model, and in the complete two-year TWFE models. RQ is positively associated with absolute recorded commitments in both components. These are associational findings, not evidence that better regulation causes financing or that cross-sectional significance confirms a within-country mechanism. The original panel and approved fsQCA findings remain unchanged.

There is substantial between-country variation in financial depth/RQ, limited within-country information and lumpy observed financing. Yet no diagnostic isolates these as the cause of the original inconclusive estimates. Complete windows reduce sample coverage, especially high income. A prespecified three-year sensitivity produces negative depth estimates with nominal p<.05, alongside inconclusive two-year and matched annual depth results. That sensitivity is fully disclosed and is not selected as a new primary finding. No unrestricted variable, window, threshold or calibration search was run.

## Frozen design and source

The bounded plan was frozen at {plan['frozen_utc']} before any Plan B estimation, followed by a separate full-design hash freeze. Input is the unchanged 2,170-row,217-economy 2014–2023 CSV `{SOURCE.name}`, SHA256 `{sha(SOURCE)}`. The annual primary sample exactly matches the original1,104 observations/137 country clusters. Each observation is a country-year; windows are country-windows, always clustered by country. No new data source or financing vintage was added.

The outcome remains the IRENA2025 PDF-derived **published recipient-country public renewable-energy financing commitments from selected institutions**, operationally million2022USD. It is not total renewable investment, domestic/private finance, disbursements or financing needs. The specific financing-table headings supply the provisional2022-price convention; the conflicting generic publication note remains unresolved. This run does not provide institutional confirmation of the price basis, commitment-recognition rules or a stable provider-year perimeter. The archived2025 API is still unavailable and 2025/2026 release differences, including UAE2022, remain unreconciled. No later-vintage values or price conversions were used. Prior clearance authorizes bounded exploratory estimation of this published construct, not broader national-finance claims.

Covariates are the original bank private credit/GDP, WGI RQ, ln(real GDP per capita, constant2015USD), ln(population), trade/GDP and CPI inflation. Depth is not overall financial development, efficiency or stability. WGI estimates have uncertainty. No IMF/FinData splice, new interaction, income-specific regression or extra mechanism was introduced.

## Within–between CRE/Mundlak specification

For every actual estimation sample, each of the six X is decomposed into x_it−mean_i(x) and mean_i(x). Means use observed estimation rows, not imputed full histories. The model is y_it = intercept + beta_W*(x_it−xbar_i) + beta_B*xbar_i + year effects + coefficients*country means of year indicators + random country intercept + error. Country means of year indicators account for the unbalanced calendar observation profile. They are not selection instruments and do not make missing reporting ignorable. Exact redundant profile columns are dropped, particularly in balanced histories.

Gaussian random-intercept CRE is estimated by nlme REML, with country-clustered CR2/Satterthwaite inference from clubSandwich. Augmented pooled Mundlak OLS is a reference, and matched TWFE is fitted on identical observations. CR2 uses identity working covariance for lm and the fitted random-intercept working covariance for lme. In this fully augmented random-intercept design, both the within coefficient and its CR2 inference reproduce matched TWFE numerically. This is an algebraic re-expression, **not an independent replication sample or new within-country evidence**. RE variance assumptions affect efficiency and between estimates; cluster robustness does not repair misspecified conditional means or endogenous selection. Standard Mundlak restriction beta_B=beta_W is tested without using that p-value to select the estimator.

Primary annual results (unadjusted95% pointwise intervals; primary-family Holm p-values shown separately):

{table(primary,pcols)}

Within-minus-between contrasts:

{table(ctshort)}

The RQ contrast has nominal p=.0468 but Holm2 p=.0936; the FD contrast is inconclusive. The six-variable HTZ joint restriction has F={jt.iloc[0].HTZ_F:.4f}, denominator df={jt.iloc[0].HTZ_df2:.4f}, p={jt.iloc[0].p_value:.4f}. These results do not robustly establish that within and between RQ slopes differ, nor justify a causal between-country story. Between RQ can reflect institutional reputation, recipient eligibility, donor/lender targeting, historical development and omitted country characteristics.

All annual reference and prespecified sensitivity results:

{table(annual,['model','term','estimate','SE_CR2','Satterthwaite_df','p_value','CI95_low','CI95_high'])}

Balanced histories retain610 observations/61 countries but only25 high-income country-years across3 ever-high-income countries; this selected restriction is not an improved four-group primary design. All-X one-year lags use frozen2014–2022 X for2015–2023 outcomes (1,003/136); no2013 values or gap-bridging interpolation were added. These timing changes do not supply instruments or eliminate anticipation.

## Multi-year definitions and matched comparisons

Main windows are2014–15,2016–17,2018–19,2020–21,2022–23. The prespecified sensitivity uses2014–16,2017–19,2020–22;2023 is omitted by design rather than padded or averaged into an unequal final window. Each horizon is analysed separately, never stacked as independent observations.

Every financing, population and covariate year must be observed and valid throughout a window. No partial-window means were computed. FD/RQ/trade/inflation use complete arithmetic annual means; logged GDPpc/population use means of annual logs. Absolute y is log1p(mean annual commitments in million2022USD). Intensity y is log1p(mean annual commitments_million*1e6/that year's population), constant2022USD/person/year. It is not a ratio of sums, a cumulative unique-project total or log of a per-window cumulative amount. Population is retained as the same essential control in both outcome models. Mean annual normalization keeps monetary units comparable across horizons. Published numeric zeros remain recorded values with rounding cautions; absent records remain missing.

Country and calendar-window fixed effects replace country and individual-year effects in aggregated models; individual-year effects cannot be separately identified from a single complete-window mean. Iterative FE eligibility requires at least2 windows per country and2 countries per window. Window context is the historical first-year income group, with changes flagged. No confirmatory group interaction is fitted.

Support, model rank and conventional residual df (not the denominator of CR2 t tests):

{table(meta[meta.model.isin(focusmodels)],['model','n','clusters','periods','rank','residual_df','tau2','sigma2','ICC'])}

Two-year data have487 country-windows/119 clusters, based on974 exact contributing country-years. Three-year data have251 country-windows/93 clusters, based on753 country-years. The full native annual1104 sample, matched annual samples, and corresponding window estimates are all reported so composition losses are not mistaken for a timing mechanism:

{table(multi,['horizon','outcome','comparison','term','estimate','SE_CR2','Satterthwaite_df','p_value','CI95_low','CI95_high'])}

Two-year depth remains inconclusive for both volume and intensity. RQ's absolute two-year association survives the frozen primary8 Holm family (adjusted p≈.0365), whereas the intensity result does not (≈.2372). Three-year depth is negative with unadjusted p≈.0218(volume)/.0229(intensity); the matched annual depth estimates remain inconclusive. RQ's three-year intervals include zero. These differences show sensitivity to aggregation/remaining within variation and functional form, not proof that timing explains any causal effect. log1p(mean amount) is not mean(log1p amount), so matched support alone cannot isolate a pure smoothing effect. The three-year design has only3 periods and8 high-income clusters, and its sensitivity tests are explicitly unadjusted exploratory results, not two independent replications of a preferred hypothesis.

## Measurement, timing and selection audit

Observed-variation decomposition:

{table(vs,['variable','observations','countries','within_SD','country_mean_SD','within_SS_share','between_SS_share'])}

Only3.71% of depth and1.51% of RQ total observed SS is within-country in the native sample; within-year identification is more restricted still. This helps explain the difference between abundant pooled variation and effective within information, but does not demonstrate attenuation, measurement error or a hidden true positive effect. Conditional-slope VIFs project out the FE/time/profile nuisance columns before normalizing; raw full-design inverse-Gram diagonals are separately labelled and must not be read as ordinary within-slope VIFs. Both within-residualized and between-country correlations are retained without choosing variables from them.

All1,399 published country-year financing records span179 economies;771 source-grid records are missing. Historical group accounting is explicit, including9 unclassified country-years (3 observed), rather than silently reclassifying them:

{table(gsummary)}

High-income588/791 historical economy-years (74.34%) have missing financing; only203 are recorded and178 enter the original complete-case FE sample. Missingness is recipient/source selection, structural scope or unavailable reporting of unknown proportions; it is not verified zero commitments. The all217-economy frame includes territories/SARs. Available-covariate medians for observed versus missing economies are descriptive and cannot correct missing-not-at-random selection. Actual high-income recipients and every observation status are listed in the audit workbook.

Income representation by specification:

{table(supportfocus)}

Group country counts overlap because historical categories change, so they do not sum to pooled country clusters. All four groups remain represented, but high-income support falls from34 ever-high-income countries in the annual benchmark to20 in complete two-year windows and8 in three-year windows. Preserving category labels does not make heterogeneity credible; no new four-group inferential claim is made.

Concentration among observed amounts:

{table(con.head(2),['scope','records','observed_amount_musd2022','largest_share','top5_share','top10_share','HHI','effective_equal_amount_units'])}

Across observed-country totals, the largest recipient contributes10.25%, top5 recipients28.02% and top10 recipients42.46%. These sums have unequal reporting histories and cover only observed amounts, not global total investment. Among{len(full)} economies with all10 financing years recorded, the median largest annual share of the observed decade is{full.largest_year_share.median():.2%}; for countries with at least5 years it is{hist.loc[hist.observed_years>=5,'largest_year_share'].median():.2%}. This is evidence of lumpy country-year totals, **not project-level concentration**. The dataset has no unique transaction IDs, individual project sizes, cancellations or provider allocation. Actual project/lender concentration, reporting-institution dependence and leave-provider-out tests cannot be calculated from these tables. Peak timing may combine project recognition, revisions and reporting, which cannot be separated here.

## Inference, influence and multiplicity

All19 planned linear models completed. Country-clustered CR2/Satterthwaite inference is available for every model; Satterthwaite df, effective design clusters and largest/top5 country information shares accompany every slope. The annual within-FD contrast has df14.23/effective15.08 clusters despite137 nominal clusters; Cambodia accounts for20.25% of its design information. Two-year FD has df11.39/effective12.19 clusters, with Cambodia23.85%. Between FD is most exposed to China (14.68% of design information). Conventional residual df1,073 in annual CRE is not1,073 independent countries. Inflation's annual within df≈1.33 warrants particular caution; all macro results are retained without trimming.

Eight primary coefficient tests were frozen as one Holm family: annual CRE W/B FD/RQ plus two-year TWFE FD/RQ for each outcome. Two within–between differences use a separate Holm2 family. Displayed95% intervals remain pointwise, not simultaneous. All other sensitivities are disclosed as unadjusted exploratory comparisons; no p-value-selected horizon or covariate replaces the baseline.

All137 annualCRE and119 two-yearTWFE country deletions were evaluated. CRE deletions hold the full-sample REML variance ratio fixed and are labelled conditional FGLS diagnostics, not fresh variance-component inference. No country was removed from the primary analysis:

{table(lsummary)}

No wild bootstrap, new PPML, nonlinear CRE or additional fsQCA was specified or run in this bounded stage; their completion is not claimed. Original PPML sensitivity and all fsQCA conclusions remain active. CR2 assumes independence across country clusters for its reference approximation and does not solve shared-donor/regional shocks, omitted effects or endogenous selection.

## Independent verification and preservation

R4.6.1 with nlme3.1-169 and clubSandwich0.7.0 is the primary implementation. Independent Python GLS, whitened CR2 geometry, Satterthwaite and HTZ calculations match all156 coefficient estimates/inferences across19 models. Independent one-dimensional concentrated REML profiles reproduce all6 fitted random-intercept variance ratios; maximum relative discrepancy is8.4e-8. The frozen original baseline coefficient/CR2 results and matched within–TWFE identities pass. R independently reconstructs annual raw-source eligibility and every738 eligible complete window; maximum aggregation discrepancy is2.3e-13. Validation tolerances, actual errors, full logs, covariances and failures are preserved. R locale-startup warnings did not prevent fitting. No fit failures were observed.

All earlier panel/fsQCA outputs, source versions, frozen datasets and protocols are preserved. No researcher-imputed observations, missing-as-zero conversion, finance-vintage splice or new calibration was introduced. `independent_verification_report.md`, `audit_evidence/final_validation.json` and `scripts/README.md` provide the numerical, preservation, workbook-readability and execution receipts. The three Excel workbooks contain all results, support/rank/df, eligible windows, missingness, concentration and diagnostics; figures are in `figures/`.

## Recommendation

Retain the original2014–2023 annual primary regression and approved fsQCA as frozen exploratory evidence. Add Plan B as a bounded explanation of variation, outcome timing and selection limits, not a replacement discovery. Do not promote the nominally negative three-year depth estimate, between-RQ significance or any calibration-sensitive absolute fsQCA pathway as validation of the original financial-development mechanism. For a broad claim about financing the whole renewable transition or credible four-income-group comparisons, prioritize a verified comparable financing dataset with project/provider-year reporting and true-zero definitions, or substantially narrow the research question to the allocation of selected published public commitments. Further specification searches on the present outcome are not the recommended next step.

Method sources: [clubSandwich lme CR2 documentation](https://jepusto.github.io/clubSandwich/reference/vcovCR.lme.html), [R nlme lme documentation](https://stat.ethz.ch/R-manual/R-devel/library/nlme/html/lme.html), and [Wooldridge's institutional paper on CRE/time-dummy means](https://harris.uchicago.edu/sites/default/files/wooldridge_ppe_seminar_paper_6-1-22_0.pdf). Source interpretation and unresolved documentation are retained in `../reports/final_methodological_clearance.md`; the [original IRENA2025 PDF](https://www.irena.org/-/media/Files/IRENA/Agency/Publication/2025/Jul/IRENA_DAT_RE_Statistics_2025.pdf) remains the defined outcome source.
'''
report+='''\n## Clarification of test scope and recorded inputs\n\nThe six-variable HTZ within-equals-between restriction is conditional on the retained calendar observation-profile terms. It is not an omnibus random-intercept exogeneity test, a test that all profile coefficients are zero, or evidence that missing financing is ignorable.\n\nComplete observed-year requirements refer to retained published values. They do not certify absence of estimation by the source agency: original WDI observation-status flags, including not-disclosed status, remain in the frozen input and annual snapshots. No new researcher-imputed value was introduced.\n'''
(OUT/'plan_b_econometric_report.md').write_text(report,encoding='utf-8')
decision='''# Research contribution and redesign decision

**Decision: preserve the existing exploratory findings; do not replace them with a preferred Plan B specification. A new verified financing dataset or a substantially narrower allocation question is necessary before making broad claims about national renewable-transition financing and income-group mechanisms.**

Plan B makes a useful diagnostic contribution. Depth and RQ vary mainly between recipients (96.29% and98.49% of observed SS), whereas the TWFE question relies on much smaller within variation. Yet the fully augmented CRE within estimates are algebraically identical to matched TWFE, not new corroboration. Banking depth is inconclusive in both components. Between-country RQ is positive, but it can reflect omitted country attributes and provider/recipient selection. The within-minus-between RQ difference does not survive its prespecified Holm2 adjustment. Neither between significance nor a nonrejected Mundlak joint restriction establishes a causal mechanism or selects the correct estimator.

The complete two-year analysis does not reveal a robust banking-depth association. Three-year sensitivities yield negative depth coefficients at nominal p≈.022, alongside inconclusive matched annual depth estimates; these are unadjusted, horizon-sensitive exploratory results. Reduced periods, nonlinear averaging and sparse high-income support prevent attributing the difference solely to project timing. Promote neither the negative three-year result nor the positive between-RQ estimate as a substitute main finding. Annual log1p million2022USD volume, per-capita intensity and approved fsQCA absolute volume remain distinct constructs and transformations.

High-income representation deteriorates: annual178 observations/34 ever-high-income countries, two-year57 windows/20, three-year18 windows/8. Strict balanced histories retain only25 high-income years across3 countries. All four labels remain present, but this is not broad four-group external validity. Historical group cluster counts overlap after income transitions. High-income missingness is588/791 economy-years; missing is not zero and cannot be repaired by an estimator.

The measurement audit supports a timing concern: among81 economies with complete10-year financing reporting, the median peak annual share is about37.2%. It does not identify concentration in individual projects or reporting institutions. Published totals lack transaction identifiers, provider-year rosters and disbursement/cancellation histories. A single frozen vintage avoids cross-vintage mixing but does not verify stable reporting within that vintage. Existing price-note inconsistencies, API unavailability and unreconciled releases remain unresolved. Robust standard errors cannot address these construct limits.

Two defensible paths remain:

1. **Narrow the present article to selected public-commitment allocation.** Make the recipient/perimeter limitation explicit in the research question, abstract, outcome definition and external-validity section. Retain the independently replicated annual benchmark, mixed/inconclusive depth and interaction findings, RQ's conditional association with PPML sensitivity, the primary fsQCA null sufficiency finding, and all predefined alternative/sensitivity results. Use Plan B as a diagnostic appendix or a tightly framed exploratory component. Avoid claims about total investment, private mobilization, disbursements or causal credit-to-finance mechanisms. A narrower title may be appropriate, subject to the authors' decision; no project title or old manuscript was changed.

2. **Acquire a financing dataset suited to the original broad question.** Require an explicit eligible-recipient frame and verified zero/nonreporting states; stable provider-year coverage; renewable-only transaction identifiers and duplicate rules; approval/recognition/disbursement dates; original and revised/cancelled commitments; currency and price-base documentation; instrument and public/private/domestic/international distinctions. Retain separate recipient/provider roles where comparability cannot be established. Reconcile overlapping records and vintages before combining amounts. A financing database must not be replaced with capacity, generation or estimated financing needs. No new database search was initiated here.

If neither a comparable financing outcome nor an honest narrower allocation construct can be defended, stop the broad financial-development/transition-finance causal or universal four-group claim. The current data can support a transparent selected-recipient associational study, not a guaranteed Q1 contribution. Credible novelty must come from a clear theory and measurement contribution with external-validity limits, not additional regressions or calibration tuning. Future causal claims would require a separately justified identification design; neither Mundlak decomposition nor multi-year averaging supplies one.

No additional fsQCA memberships, truth tables or solutions were run; the111-country approved results are unchanged. No new income interactions, PPML or threshold search was performed. All19 Plan B models and their null/mixed/sensitive results are retained, with CR2 inference, primary-family Holm adjustments, matched-support comparisons and independent numerical verification.
'''
(OUT/'research_contribution_and_redesign_decision.md').write_text(decision,encoding='utf-8')
independent=f'''# Independent numerical verification — bounded Plan B

Primary estimates: R4.6.1, nlme3.1-169 Gaussian random-intercept REML and clubSandwich0.7.0 CR2/Satterthwaite; HTZ restriction tests. Independent implementation: Python NumPy/SciPy whitened GLS and residual-space CR2, reusing the unchanged previously benchmark-validated inference formulas through a read-only import. This does not re-use R's covariance/adjustment matrices. FullR logs, design columns and source hashes are preserved.

All19 planned models completed, with156 slope coefficient comparisons. For random intercept covariance sigma²I+tau²J, Python uses the analytic symmetric inverse-square-root transformation within each cluster; global scale cancels from GLS/CR2. An independent scalar profile REML optimization recomputes the variance ratio, rather than treating the nlme optimizer as independently verified merely by recycling its parameters. The coefficient/covariance verification uses R's fitted ratio to test identical estimators, with the profile check separately assessing variance fitting.

{table(checks)}

R independently reconstructs annual source eligibility and aggregate means from the original raw-source columns, rather than merely reading Python's final outcomes:

{table(agg)}

Maximum coefficient discrepancy is1.7e-12, SE discrepancy1.3e-9, Satterthwaite df discrepancy4.0e-8, and CI endpoint discrepancy2.9e-9. Joint tests agree within2.8e-12 for denominator df and8.1e-15 for p. Independent profile REML variance ratios differ by at most8.4e-8 relative to max(1,Rratio). Frozen-original baseline coefficient/CR2 comparisons and matched within/TWFE identities pass. No material discrepancy was silently adjusted.

Country-deletion CRE uses fixed full-sample variance ratios and is a conditional FGLS diagnostic, not an independent nlme refit for every deletion. Aggregation source checks and coefficient comparisons do not verify that the underlying institutional records are a complete financing census. CR2's working covariance and country-independence assumptions remain explicit. No bootstrap or nonlinear PPML claim is made. Locale startup warnings are recorded and did not affect successful verification.

The final separate file-validation receipt checks every XLSX data cell, readable figures/scripts, freeze hashes and protected prior outputs. It is not a claim of native Excel application execution or independent source-transaction verification.
'''
(OUT/'independent_verification_report.md').write_text(independent,encoding='utf-8')
readme='''# Reproduce the bounded Plan B analysis

All outputs live in Renewable_Energy_Finance_Q1/plan_b_analysis. Originals are read only. The frozen source is ../reports/final_model_selection_20261016/four_income_group_final_dataset.csv. The exact bounded model/window plan and full-design hashes precede estimation in audit_evidence. No financing download, vintage splice, new variable, fsQCA or threshold search is part of this run.

Run the following scripts from this directory with Python3.13.9 (see dependency_versions.json), preferably with BLAS threads set to1. Use the Python executable from your own environment. Scripts explicitly configure the project root; change it only for an intentionally relocated project copy.

1. planb_prepare.py — hash prior outputs, freeze bounded plan, derive native/lag/balanced support and complete nonoverlap windows, and freeze model designs. Does not estimate coefficients.
2. planb_run_R.py — executes planb_estimate.R using the preserved R4.6.1 runtime and existing nlme3.1-169/clubSandwich0.7.0 libraries. All19 fits and full logs are captured by Python into real files. No packages need downloading.
3. planb_verify.py — independently profiles REML and reproduces GLS/OLS coefficients, CR2, Satterthwaite and HTZ; checks original benchmark and matched TWFE; applies frozen Holm families and country-deletion diagnostics. Requires all checks to pass.
4. planb_run_aggregation_check.py — executes planb_aggregation_check.R for independent raw-source annual eligibility and complete-window mean checks, capturing the CSV and full log automatically.
5. planb_outputs.py — descriptive diagnostics, full-precision JSON workbook payload and figures; no models estimated.
6. planb_reports.py — source-backed main report, redesign decision and independent report.
7. node build_planb_workbooks.mjs — authors three workbooks with the bundled @oai/artifact-tool runtime, scans formula errors, renders every worksheet and exports. node_modules is a junction to the immutable bundled Node dependencies; no third-party Excel writer is used.
8. planb_validate.py — read-only XLSX cell roundtrip, numerical/source/window/freeze validation, protected-file SHA256 checks, readable PNG/PDF and script syntax checks. See audit_evidence/final_validation.json.

Original Gaussian models: y=log1p(million2022USD annual commitments); six original X. CRE splits all X into within and between means and adds country means of calendar indicators to year effects for unbalanced profiles. lme working CR2 target is its fitted covariance; lm target identity. Full conditional slope VIFs project out nuisance columns before normalization. Raw full-design inverse-Gram diagonals are separately labelled and are not ordinary slope VIFs.

Two-year windows2014–15/16–17/18–19/20–21/22–23; three-year sensitivity2014–16/17–19/20–22. Every financing/population/X year is required. Window amounts are means in the original annual monetary unit; per-capita outcome averages annual ratios before log1p. No missing/partial window becomes zero or a complete mean. No partial-window-mean sensitivity was computed.

Read completed files without rerunning the pipeline. Intentional reproduction may re-export this new results folder; old datasets, panel/fsQCA results and protocols must stay unchanged. No wild bootstrap/PPML/additional fsQCA is claimed. The last receipt records all completed checks and methodological limits.
'''
(SCRIPTS/'README.md').write_text(readme,encoding='utf-8')
if Path(__file__).resolve()!=(SCRIPTS/'planb_reports.py').resolve():shutil.copy2(Path(__file__),SCRIPTS/'planb_reports.py')
versions=json.loads((SCRIPTS/'dependency_versions.json').read_text());versions['nlme']='3.1-169';(SCRIPTS/'dependency_versions.json').write_text(json.dumps(versions,indent=2),encoding='utf-8')
print('Main report, redesign decision, independent report and execution instructions saved',flush=True)
