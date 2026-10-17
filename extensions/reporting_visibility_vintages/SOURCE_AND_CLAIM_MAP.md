# Source-to-claim and manuscript-to-code map

The public repository complements the original 2014–2023 panel and the separately reported fsQCA; it does not supersede the frozen 2025 core estimand.

| Manuscript item | Public artifact | Source and limitation |
|---|---|---|
| Table 6, source-record presence regression | `results/visibility_logit.csv`, `results/visibility_within_lpm.csv`, `results/visibility_predictions.csv` | Conditional logit N=1,522/163; same-support country/year FE LPM. Neither identifies provider choice or causal effects. |
| Source visibility descriptive discussion | `results/visibility_coverage.csv` | 2,170-cell reference grid, income-group counts and Wilson limits. |
| Table 7, vintage source-record continuity | `results/vintage_coverage.csv`, `results/vintage_by_year.csv` | 1,953 matched cells, 2015–2023. Interpret absence as absence in that release, not no finance. |
| Figure 7 | `results/vintage_coverage.csv` | Group-level coverage rates are calculable from aggregates. |
| Table 8, matched FE models | `results/matched_vintage_FE.csv` | 863 matched cells; CR1; vintage-specific native price bases; no exchange-rate/deflator bridge. |
| Vintage rank statement | `results/vintage_rank_by_year.csv`, `results/analysis_receipt.json` | Spearman on 1,099 positive common records only. |
| Model denominators and design | `results/analysis_receipt.json` | Original numerical receipt (not a fresh raw-source audit). |

**Frozen inputs:** `four_income_group_final_dataset.csv` SHA256 `5dfe182725b16192ca00d93d7c34e900471e8b11855542fd6319cf08b8b958e6`; `irena_vintage_overlap.csv` SHA256 `fdd711c71b851b50754d4bde3f266bd497532f60db922a122c49b43b1b40b5af`. Neither is publicly redistributed here.

For reproducibility, `analyze_reporting_and_vintages.py` generates the result tables from exact source input versions. The baseline coefficient and fsQCA public outputs remain in the existing `evidence/` directories.

**Independence boundary:** The previously reported reruns reproduce coefficients from the author-frozen inputs. They do not independently verify the institutional transaction population or underlying vintage revision causes. An external researcher needs lawful access to the frozen datasets or a comparably documented independently reconstructed panel.
