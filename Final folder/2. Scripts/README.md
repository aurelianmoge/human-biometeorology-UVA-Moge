# Script reading order

Run `run.py` to execute the following stages in order. Every R script must run from the clean project root; the runner sets that directory automatically.

| Theme | Script | Purpose |
| --- | --- | --- |
| Data | `01_collect_sources.py` | Initial immutable input/evidence collection and hash manifests; only stage needing Exploration |
| Data | `02_extract_workbooks.py` | Preserve original Excel cell provenance and corroborated date repairs |
| Data | `03_assemble_panels.R` | Continuous calendars, duplicate/missingness policy, physical weather checks, fixed station proxies |
| Data | `04_census_populations.py` | Correct July year codes, infant handling, population partitions and Bedford harmonization |
| Descriptions | `01_describe_data.R` | Regional/annual/seasonal summaries, cloud plots, correlations and county descriptions |
| Models | `01_fit_models.R` | Explicit 56-request count-model registry, cached fits and diagnostics |
| Models | `02_minima_and_thresholds.R` | Finite-domain minima, conditional minimum draws and fixed published-event rules |
| Models | `03_blocked_validation.R` | Buffered internal prediction holdouts, with per-day scoring |
| Models | `04_segmented_thresholds.R` | Descriptive same-day hinge profiles; no invented selection-adjusted tests |
| Research | `01_direct_comparisons.R` | Paired demographic/metric contrasts, historical arithmetic reproduction and temporal covariance sensitivity |
| Research | `02_classification_accounting.R` | Complete pooled/own qualifying sets and signed date changes |
| Research | `03_county_inequalities.R` | Annual ecological log-rate regressions with county-clustered covariance |
| Research | `04_gam_interactions.R` | Genuine penalized same-day additive/tensor GAMs and supported surfaces |
| Research | `05_population_allocations.R` | Corrected Census-weighted allocation scenarios, not verified catchment rates |
| Research | `06_research_atlas.R` | Subgroup/cause/specification/lag plots, including uncertain and null results |
| Shared | `verify_results.R` | Numerical identities and independent historical checkpoints |
| Shared | `render_tables.R` | Readable paginated PNGs for analytical tables up to 500 rows |
| Shared | `build_reports.py` | Detailed reports, equations and offline visual browser |

The shared engine deliberately uses readable explicit basis operations. It does not install packages or source legacy scripts. The `model_registry.csv`, `model_status.csv`, `verification_checks.csv`, manifests and stage logs document what ran and what stopped. Large daily and curve CSVs are plotted rather than turned into thousands of table-image pages.

Specialized inherited branches retain their output/report snapshots and scientific limits. The clean runner does not claim to rerun every exploratory search or reconstruct missing original random draws.
