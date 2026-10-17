# Reporting visibility and IRENA 2025–2026 vintage sensitivity

**Article:** *The Measurement Boundaries of Renewable-Energy Public Finance: Reporting Visibility and Institutional Associations*

**Scientific status:** Exploratory, post-primary extension. These analyses do not estimate causal effects or national renewable investment. Missing published records are **not** national zero-finance amounts; observed numeric zeros remain observations.

## What is deposited
- `analyze_reporting_and_vintages.py`: original extension analysis program. It takes two authorized local inputs and writes eight aggregate result tables plus `analysis_receipt.json`.
- `RESULTS_SHA256.csv`: SHA256 and byte-length manifest for all nine aggregate results/receipt files, validated by `verify_aggregate_results.py`.
- `render_aggregate_figures.py`: regenerates manuscript Figures 6 and 7 directly from published aggregate tables.
- `results/`: published **aggregate model results** for the two new research margins: presence of a 2025 IRENA recipient record and continuity between IRENA 2025/2026 source vintages.
- `SOURCE_AND_CLAIM_MAP.md`: identification of each input, corresponding manuscript claims and limits.
- `requirements.txt`: environment used for the frozen extension.

## Required inputs (not distributed)
1. `four_income_group_final_dataset.csv` — SHA256 `5dfe182725b16192ca00d93d7c34e900471e8b11855542fd6319cf08b8b958e6`, 2,170 economy-year grid, 2014–2023. From the private evidence archive at `evidence/reports/final_model_selection_20261016/`.
2. `irena_vintage_overlap.csv` — SHA256 `fdd711c71b851b50754d4bde3f266bd497532f60db922a122c49b43b1b40b5af`, 1,953 matched economy-year cells, 2015–2023. From `evidence/reports/alternative_evidence/` in the private archive.

**Run after verifying source rights and exact hashes**:
```sh
python -m pip install -r requirements.txt
python verify_frozen_inputs.py /authorized/four_income_group_final_dataset.csv /authorized/irena_vintage_overlap.csv
python analyze_reporting_and_vintages.py --panel /authorized/four_income_group_final_dataset.csv --vintage /authorized/irena_vintage_overlap.csv --outdir reproduced
python verify_aggregate_results.py
python render_aggregate_figures.py
```
`verify_aggregate_results.py` validates deposited aggregate identities only; it does **not** recalculate missing licensed source data. A genuine analysis rerun requires the two inputs above.

## Scope and interpretation
- Pooled record-visibility binomial GLM: 1,522 complete cases and 163 country clusters. Logit associations are not replicated in country-FE LPM.
- Matched-vintage 2015–2023 reference: 1,953 cells, 1,121 in both, 142 only in 2025, 13 only in 2026, 677 in neither.
- Two separate log(1+financing) FE sensitivity models on 863 matched cases and 116 country clusters. They use the vintages' **different price bases** (2022-price vs 2023-price USD), so coefficients are not estimates of a common dollar change.
- Published rank correlation among 1,099 jointly positive records does not remedy missingness in the target population.
- Provider-by-year rosters, transaction-level revision history and complete harmonized amounts are unavailable. Do not impute missing records or infer causal provider churn.

## Availability / versioning caveat
This extension is not a complete one-click raw-data replication archive. The GitHub repository contains aggregate outputs and executable code only; institutional source extracts are omitted while redistribution rights are reviewed. The original analysis and its independent partial checks were performed against privately held frozen files. No agency raw files or country-level memberships were uploaded by this change.
