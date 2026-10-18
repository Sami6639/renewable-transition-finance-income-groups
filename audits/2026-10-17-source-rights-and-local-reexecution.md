# Source licenses and frozen-input re-analysis boundary (17 October 2026)

**Final manuscript title (18 October 2026):** Renewable Energy Public Finance: Financial Development, Institutions, and Income Groups. The former title (*The Measurement Boundaries of Renewable-Energy Public Finance: Reporting Visibility and Institutional Associations*) was a working title at the date of this historical audit.

## Re-execution verified on the author's frozen research archive

A local rerun (Python 3.13, frozen October 2026 country-year inputs) reproduced:
- Original econometrics: 60 OLS/CR2 models, 5 Rademacher null-restricted wild-cluster bootstrap tests (9,999 draws), 3 PPML/CR1 models.
- Exploratory mechanisms: 22 OLS/CR2 models, 4 bootstrap tests, and 6 PPML/CR1 specifications.
- Configuration: 275 Python fsQCA calibration/sufficiency sensitivity scenarios; 550 solution output rows.
- New reporting-visibility / 2025–2026 vintage sensitivity: 8 aggregate published tables reproduced from the frozen inputs.
- Comparative numerical ledgers: primary OLS 482 coefficient rows, secondary OLS 148 coefficient rows, fsQCA 550 solution rows. Maximum absolute differences for inspected estimands below 1.1e-10 (and exact matching in bootstrap p-values and fuzzy solution metrics).

**Scope:** These checks reproduce estimators and outputs *from the identical frozen derived panels*, **not** fresh acquisition and independent transaction-level reconstruction from institution databases. No new independent R run was conducted in this runtime (Rscript unavailable); previously archived R checks remain historical verification, not a newly rerun R audit. The statistical rerun does not verify the source-vintage provider roster or any transaction-level inclusion decision. **Metadata update 18 October 2026:** the archived IRENASTAT 2025 matrix PUBFIN_2025_H2_PX explicitly resolves the operational source-unit identification as 2022-price and exchange-rate USD using donor-country deflators; the yearbook generic note saying current prices remains an unresolved textual inconsistency.

## Copyright, terms, and dataset segmentation

| Source | Published terms relevant to sharing | Release position |
| --- | --- | --- |
| IRENA *Renewable Energy Statistics 2025* | Publication copyright page says material can be freely used, shared, copied or reproduced if clearly attributed and carrying the © IRENA 2025 notice; third-party material may have separate rights. [Publication PDF](https://www.irena.org/-/media/Files/IRENA/Agency/Publication/2025/Jul/IRENA_DAT_RE_Statistics_2025.pdf), [website terms](https://www.irena.org/terms-and-conditions). | Derived aggregates and code are published with attribution. Public bulk release of the mixed source-compiled financing panel is withheld pending third-party-rights and source-field review. |
| World Bank WDI | [Dataset terms](https://data.worldbank.org/summary-terms-of-use): default CC BY 4.0 except indicator-specific restrictions. | Eligible open data may be shared with proper attribution and change notation; verify each frozen indicator's metadata. |
| World Bank Worldwide Governance Indicators | [Data Catalog](https://datacatalog.worldbank.org/search/dataset/0038026/worldwide-governance-indicators) lists CC BY 4.0. | Attribute World Bank, WGI vintage and revisions. |
| World Bank RISE | [Data Catalog](https://datacatalog.worldbank.org/search/dataset/0038387/regulatory-indicators-for-sustainable-energy-rise-application) lists CC BY 4.0. | Attribute World Bank/ESMAP and exact data vintage. |
| IMF statistics | [IMF Copyright and Usage](https://www.imf.org/en/about/copyright-and-terms) permits statistical-data distribution subject to data-specific attribution, accuracy and use conditions; additional commercial permissions may apply. | Do not assume all IMF website content falls under the special statistical data clause. |
| OECD-DAC component (if incorporated in IRENA financing tables) | [OECD Data Terms](https://www.oecd.org/en/about/terms-conditions.html) permit reuse of OECD-owned data with attribution, subject to metadata and third-party restrictions. | Attribute relevant provider and obtain clarification if the incorporated finance series' contribution/ownership is unclear. |

These are public source terms, not an assertion that a combined derivative panel is universally licensed or that an irrevocable license is granted by the repository owner.

## Open items for institutional confirmation

The IRENA 2025 PDF financing tables are headed in 2022 USD, but its generic finance methodological note describes current-price USD figures and rounding conventions inconsistent with the displayed decimals. The country-level funding source identifies several institutional reporters but lacks transaction-level cross-vintage reconciliation. Requests should seek (i) the authoritative price basis for country public-finance tables, (ii) the financing-data provider scope and ownership/redistribution eligibility, and (iii) availability of a frozen historical machine-readable export for reproducing the 2025 publication.

Contact identified in the IRENA 2025 publication: **statistics@irena.org** (no correspondence has been sent by this repository).

## Repository boundary

- Codes and aggregate outputs are public; country-year institutional downloads and mixed derived panels remain local while license questions are unresolved.
- `verify_public_release.py` validates existing public release files and core published identities; `extensions/reporting_visibility_vintages/verify_aggregate_results.py` validates new extension outputs and hashes.
- Running these public checks is **not** equivalent to rerunning all models.
- Authorized researchers with access to the original frozen datasets may execute the original scripts and record new code/data/environment hashes.
- Scientific interpretations remain associational; no evidence is fabricated to fill nonreported financing.

**Attribution for finance data:** IRENA (2025), *Renewable Energy Statistics 2025*, International Renewable Energy Agency, Abu Dhabi. © IRENA 2025. Data attribution does not imply IRENA endorsement.
