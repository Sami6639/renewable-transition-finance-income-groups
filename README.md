# Renewable transition finance and income groups

Code and aggregate research outputs for **Financing the Renewable Energy Transition: Financial Development, Institutional Quality, and Income-Group Heterogeneity**.

## What is public here

This repository contains research Python/R code, aggregate model coefficients and covariance matrices, model-fit and inference tables, selected aggregate figures, methodological/calibration plans, dependency records and verification utilities. No new project-wide license is granted. Source-data terms remain applicable.

Country-year source/derived panels, country-level fuzzy memberships, estimation membership lists, model design matrices, case-support records and raw institutional downloads are intentionally omitted while their redistribution terms are reviewed. No manuscripts, personal paths, private approval text, editorial task reports or third-party source reports are included. The code is therefore not a self-contained reproduction environment.

The study concerns selected public renewable-energy commitments in the frozen IRENA 2025 vintage for 2014–2023, provisionally in million constant-2022 USD. Missing financing records must not be recoded as zero; no 2026 financing values should be substituted.

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

Earlier private preparation verified the original source archives, reconstructed the 1,104-row/137-country sample and independently reproduced baseline OLS coefficients and five fsQCA calibrations. That work is not a fresh end-to-end run of this public package. CR2, bootstrap, PPML, all minimizations, second-stage and Plan B pipelines have not been rerun after portability changes.

## Execution limitations

Python project roots resolve relative to the scripts. Original numeric environment: Python 3.13.9 and recorded package versions. R pipelines require separately installed recorded QCA/nlme/clubSandwich versions; set RSCRIPT if Rscript is not on PATH and REPLICATION_EVIDENCE_ROOT for the secondary independent R script. No external R runtime, workbook builders, source downloads or historical protection manifests are bundled.

Missing-input and scientific hash checks were not bypassed. Historical setup/protection routines may fail on missing files or sanitized-script hash differences; these must be resolved explicitly in an intentional reproduction. Original administrative steps copying private chat-approval text were removed; no estimator, sample rule, threshold or numerical result was changed.
