# Renewable Energy Public Finance: Financial Development, Institutions, and Income Groups

## Research extension for the revised manuscript (October 2026)

The authoritative, final manuscript title is **Renewable Energy Public Finance: Financial Development, Institutions, and Income Groups**. Earlier titles in preserved historical audit records are working titles and do not supersede this title. In addition to the frozen 2014–2023 analyses below, its new reporting-visibility and IRENA 2025/2026-vintage comparisons are publicly documented in [extensions/reporting_visibility_vintages/](extensions/reporting_visibility_vintages/).

The new extension includes the exact regression program, eight aggregate results CSV files, an analysis receipt, a code-to-table mapping, frozen input SHA256 identifiers, aggregate verifier, input-hash verifier and figure regeneration script. **Authorized underlying country-year files remain excluded**, so full end-to-end reproduction still requires obtaining the documented exact source inputs. Run:

```sh
python verify_public_release.py
python extensions/reporting_visibility_vintages/verify_aggregate_results.py
```

These two public checks verify published aggregates and the original public release manifest, **not** the correctness of institutional source transactions or the source-dependent estimation itself.


This repository supports the final manuscript named above, retaining the original 2014–2023 core panel and the exploratory reporting-visibility / data-vintage extensions. Historical artifacts may contain earlier working titles.

## What is public here

This repository contains research Python/R code, aggregate model coefficients and covariance matrices, model-fit and inference tables, selected aggregate figures, methodological/calibration plans, dependency records and verification utilities. No new project-wide license is granted. Source-data terms remain applicable.

Country-year source/derived panels, country-level fuzzy memberships, estimation membership lists, model design matrices, case-support records and raw institutional downloads are intentionally omitted while their redistribution terms are reviewed. No manuscripts, personal paths, private approval text, editorial task reports or third-party source reports are included. The code is therefore not a self-contained reproduction environment.

The frozen 2014–2023 outcome is selected recipient-country public-finance commitments in **million USD at 2022 prices and exchange rates**, using donor-country deflators, according to the contemporaneous [IRENASTAT 2025 finance-table metadata](https://pxweb.irena.org/pxweb/en/IRENASTAT/IRENASTAT__Finance/PUBFIN_2025_H2_PX.px/) (updated 23 July 2025; matrix PUBFIN_2025_H2_PX). The yearbook's generic note instead says *current prices*: this contradictory boilerplate is disclosed, but the dated series-specific metadata is used for the unit. The source combines OECD-DAC 2024 statistics and IRENA Public Finance Database 2025; the yearbook describes 16 major development-finance institutions in addition to OECD-DAC. The detailed provider-by-year roster is **not** verified, and the data do not constitute total national investment. Missing records are not zeros. The separate 2026 release uses **2023-price USD** and is not spliced into the 2025 sample. See [source unit and reporting-boundary record](DATA_SOURCE_SCOPE_2025.md).

## Results and organization

- evidence/econometric_results/audit_evidence/: baseline and interaction coefficients, CR2 covariances, bootstrap tests, PPML checks, descriptive/group aggregates and fit metadata.
- evidence/second_stage_econometric_results/audit_evidence/: secondary-model aggregate estimates and comparisons.
- evidence/plan_b_analysis/audit_evidence/: within/between and time-window aggregate estimates and diagnostics.
- evidence/fsqca_results_20261016/audit_evidence/: necessity and solution metrics, without country memberships or case-support tables.
- evidence/fsqca_preparation_20261016/: case-variable definitions and aggregate selection counts; the country-case dataset is omitted.

Manuscript baseline FD and RQ estimates correspond to model='Primary Baseline' in econometric_results/audit_evidence/coefficients.csv. Primary/alternative QCA solution metrics correspond to primary_Ypc and primary_Yabs in Python_solutions.csv. Model metadata and inference files should accompany coefficients; values are not causal effects.

All published numerical outputs are copied unchanged from the verified source archive. Selected code has path/privacy adaptations only, recorded in PORTABILITY_CHANGELOG.json. Original scripts can write into their results directories; use a separate working copy for intentional reproduction.

## Verification

Run the standard-library-only public check from the repository root:

```sh
python verify_public_release.py
```

It validates MANIFEST.csv hashes and selected published-result identities. It does not re-estimate models or validate omitted source data. A successful check certifies only those bounded checks.

The separately included verify_with_local_inputs.py performs partial independent sample/OLS/calibration checks only after the required original frozen inputs are placed at their documented paths. It requires numpy, pandas and threadpoolctl. See ACQUIRE_INPUTS.md. It is expected to stop with FileNotFoundError while these inputs are absent.

Earlier private preparation verified the original source archives and reconstructed the 1,104-row/137-country sample. A **separate local independent R execution on 17 October 2026** subsequently reported 114/114 passing comparison groups for selected CR2/Satterthwaite, HTZ, Mundlak/timing, and 275 fsQCA scenarios. Those original R receipts are held locally and have **not** been released here. PPML, bootstrap, raw-to-panel construction and reporting-visibility/vintage estimators were not independently rerun in R. This public checkout does not contain the private input panels necessary to execute the entire study from raw transactions; the historical code portability adaptations were not retested end-to-end.

## Execution limitations

Python project roots resolve relative to the scripts. Original numeric environment: Python 3.13.9 and recorded package versions. R pipelines require separately installed recorded QCA/nlme/clubSandwich versions; set RSCRIPT if Rscript is not on PATH and REPLICATION_EVIDENCE_ROOT for the secondary independent R script. No external R runtime, workbook builders, source downloads or historical protection manifests are bundled.

Missing-input and scientific hash checks were not bypassed. Historical setup/protection routines may fail on missing files or sanitized-script hash differences; these must be resolved explicitly in an intentional reproduction. Original administrative steps copying private chat-approval text were removed; no estimator, sample rule, threshold or numerical result was changed.
