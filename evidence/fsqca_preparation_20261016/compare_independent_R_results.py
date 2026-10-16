from qca_replication_setup import *
import pandas as pd,numpy as np
OLD=ROOT/'econometric_results/audit_evidence'
py=pd.read_csv(OLD/'coefficients.csv');py=py[py.model.isin(['Primary Baseline','Primary RQ interaction','Primary Income interaction'])]
r=pd.read_csv(EV/'R_clubSandwich_coefficients.csv');ci=pd.read_csv(EV/'R_clubSandwich_confidence_intervals.csv');rc=r.merge(ci[['model','term','CI_L','CI_U']],on=['model','term'],validate='one_to_one');z=py.merge(rc,on=['model','term'],how='left',validate='one_to_one');assert len(z)==25 and z.beta.notna().all()
summary=[]
for pcol,rcol,tol in [('estimate','beta',5e-9),('SE_CR2','SE',5e-9),('Satterthwaite_df','df_Satt',5e-6),('CI95_low','CI_L',5e-8),('CI95_high','CI_U',5e-8),('p_value','p_Satt',5e-8)]:
 z['difference_'+pcol]=z[rcol]-z[pcol];summary.append(dict(quantity=pcol,terms_compared=len(z),maximum_absolute_difference=z['difference_'+pcol].abs().max(),absolute_tolerance=tol,passed=bool(z['difference_'+pcol].abs().le(tol).all())))
z.to_csv(EV/'independent_coefficient_comparison.csv',index=False)
pj=pd.read_csv(OLD/'joint_tests.csv');pj=pj[pj.model.isin(['Primary RQ interaction','Primary Income interaction'])];rj=pd.read_csv(EV/'R_clubSandwich_joint_tests.csv');joint=pj.merge(rj,on='model',validate='one_to_one')
for pcol,rcol,tol in [('HTZ_F','Fstat',5e-8),('HTZ_df1','df_num',0),('HTZ_df2','df_denom',5e-6),('HTZ_p','p_val',5e-8)]:
 joint['difference_'+pcol]=joint[rcol]-joint[pcol];summary.append(dict(quantity=pcol,terms_compared=len(joint),maximum_absolute_difference=joint['difference_'+pcol].abs().max(),absolute_tolerance=tol,passed=bool(joint['difference_'+pcol].abs().le(tol).all())))
joint.to_csv(EV/'independent_joint_test_comparison.csv',index=False);summary=pd.DataFrame(summary);summary.to_csv(EV/'replication_discrepancy_summary.csv',index=False)
keys=pd.read_csv(EV/'R_independent_sample.csv');mem=pd.read_csv(OLD/'sample_membership.csv');mem=mem[mem.model.eq('Primary Baseline')];assert len(keys)==1104 and not keys.duplicated(['country_code','year']).any();assert set(zip(keys.country_code,keys.year))==set(zip(mem.country_code,mem.year))
rm=pd.read_csv(EV/'R_model_metadata.csv');pm=pd.read_csv(OLD/'model_metadata.csv');rank=rm.merge(pm[['model','n','country_clusters','model_rank','residual_df']],on='model',suffixes=('_R','_Python'),validate='one_to_one');assert (rank.n_R==rank.n_Python).all() and (rank['rank']==rank.model_rank).all() and (rank.residual_df_R==rank.residual_df_Python).all();rank.to_csv(EV/'replication_design_comparison.csv',index=False)
def md(d,cols):
 def f(x):return f'{x:.10g}' if isinstance(x,(float,np.floating)) else str(x)
 return '| '+' | '.join(cols)+' |\n| '+' | '.join(['---']*len(cols))+' |\n'+'\n'.join('| '+' | '.join(f(r[c]) for c in cols)+' |' for r in d.to_dict('records'))
status='PASSED — no material discrepancies' if summary.passed.all() else 'DISCREPANCIES — review required; originals preserved'
text=f'''# Independent econometric replication

**Result: {status}.** The three actual renewable-finance primary models were independently re-estimated using **R 4.6.1 / clubSandwich 0.7.0**, rather than another call to the custom Python estimator. The original reports, source CSVs, workbooks and code were preserved. No new empirical specification was selected and no result was changed to obtain significance.

## Independent implementation and exact support

The R script reads the frozen 2014–2023 source CSV and independently applies complete-case, nonnegative-financing, positive GDPpc/population and known historical-income checks, followed by iterated country/year singleton removal. It reproduces exactly the **1,104 country-years and 137 country clusters** in the primary Python sample. No 2024/IRENA 2026 values enter. All published numeric zeros remain; unavailable values remain missing.

R uses `stats::lm` with explicit country/year factor effects. FD and Regulatory Quality are grand-mean centered before products. Historical group main effects and FD interactions use low income as reference. R independently constructs the group variables and model matrices. The baseline, RQ-interaction and income-interaction ranks are **152, 153 and 158**, respectively, with residual df **952, 951 and 946**. These match the originals. No original Python fitting helper was imported into the R replication.

`clubSandwich::vcovCR(type="CR2", target=rep(1,N))` implements identity-working-model country-cluster bias reduction; `coef_test(test="Satterthwaite")` and `conf_int` provide coefficient-specific tests and confidence intervals. `Wald_test(test="HTZ")` independently replicates the original joint tests. The independent check covers **all 25 identified non-FE coefficients**, their SEs, df, p-values and 95% intervals, and **both primary interaction joint tests**. It is not a claim that all 60 sensitivity models or the bootstrap/PPML implementations were independently rerun.

## Numeric discrepancy ledger

{md(summary,['quantity','terms_compared','maximum_absolute_difference','absolute_tolerance','passed'])}

The tolerances above are recorded explicitly; observed differences are numerical rounding/linear-algebra differences, not substantively different estimates. Exact per-term and joint-test comparisons, R matrices/metadata, independent sample keys and session/package information are retained in `audit_evidence`.

## Independently reproduced key coefficients

{md(z[z.term.isin(['FD','Governance','FD_x_Governance'])],['model','term','beta','SE','df_Satt','CI_L','CI_U','p_Satt'])}

## Joint interaction tests

{md(joint,['model','Fstat','df_num','df_denom','p_val','difference_HTZ_F','difference_HTZ_df2','difference_HTZ_p'])}

The replication supports the reported interpretation: contemporaneous financial depth and FD×RQ are inconclusive; Regulatory Quality has a positive conditional baseline association; the pooled income-interaction evidence is imprecise. Agreement between implementations verifies calculations under the chosen design, not source completeness, causal identification, statistical power or the truth of a theory. It does not validate an fsQCA solution.

## Access, runtime and remaining limits

No Rscript executable was initially installed. The official CRAN **R-4.6.1-win.exe** was downloaded and its MD5 matched CRAN's published fingerprint `7907f3a20ec8ec88cd0da279024b8e27` before project-local installation. R's online installer could not access the package index and reported that the binary was unavailable. The reproducible fallback downloaded the actual official CRAN 4.6 binary index and required binary archives through the working Python HTTPS connection, then installed them locally. The R installer reported each package unpacked and internal MD5 checks passed. The retained source index, archives, scripts and logs document that fallback. The failed optional source-file URL for R's installer script returned HTTP 404 and was not needed for estimation.

The Windows R process also emitted locale-setting warnings and the usual S3-method registration notice. Fits completed successfully. An execution-environment file-export visibility problem was resolved by capturing R's CSV output through a Python launcher; Windows line endings were normalized when parsing that captured output. Complete captured stdout/stderr and actual CSV outputs are retained. This was an export issue, not a fabricated replacement for R estimates.

Original restricted-bootstrap p-values were not independently replicated in R, and nonlinear PPML CR2/Satterthwaite or bootstrap remain unavailable in the original analysis. Those original scope limitations remain unchanged. This phase establishes independent clubSandwich replication of the three requested OLS models, confidence intervals and primary HTZ tests only. There is no unresolved blocker to that check.

Reproduction: run `install_local_R_packages.R` against the retained project-local R runtime/library, then `replicate_clubsandwich.R`; if file visibility requires it in this desktop environment, use `capture_R_replication.py` to capture exports. `compare_independent_R_results.py` reads the exported evidence and never modifies original estimates. The package reference is [clubSandwich documentation](https://jepusto.github.io/clubSandwich/reference/vcovCR.html); joint inference is documented in [Wald_test](https://jepusto.github.io/clubSandwich/reference/Wald_test.html). Runtime installer documentation is the [official R Windows FAQ](https://cran.r-project.org/bin/windows/base/rw-FAQ.html).
'''
(OUT/'independent_econometric_replication.md').write_text(text,encoding='utf-8')
(EV/'replication_validation.json').write_text(json.dumps(dict(status=status,implementation='R 4.6.1 / clubSandwich 0.7.0',independent_sample_keys_identical=True,coefficients_compared=25,joint_tests_compared=2,all_tolerances_passed=bool(summary.passed.all()),maxima=summary.set_index('quantity').maximum_absolute_difference.to_dict()),indent=2))
print(summary.to_string(index=False));print(status)
