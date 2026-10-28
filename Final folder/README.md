# Virginia human biometeorology — clean project

Start with **[the general report](4.%20Reports/GENERAL_REPORT.html)** ([PDF](4.%20Reports/GENERAL_REPORT.pdf)) or **[the visual results browser](3.%20Outputs/index.html)**. The four thematic reports explain the data, descriptions, models and original research. A [mathematical appendix](4.%20Reports/MATHEMATICAL_APPENDIX.html) develops the intuition, equations and assumptions. The [historical evidence review](4.%20Reports/HISTORICAL_EVIDENCE_REVIEW.html) distinguishes reproduced results from preserved specialized branches and later corrections.

## Four folders

| Folder | Contents |
| --- | --- |
| `1. Data` | Immutable raw copies, derived panels, corrected Census cells, source manifests, variable/crosswalk dictionaries and frozen numerical references |
| `2. Scripts` | Numbered, commented English scripts, shared statistical functions and a single runner |
| `3. Outputs` | Four matching themes, each with CSV, PNG, paginated table PNGs, model objects, logs and clearly labeled historical evidence |
| `4. Reports` | Four thematic reports, general report, mathematical appendix, historical evidence review; Markdown, HTML and PDF |

The four themes are data description, descriptive statistics, modelization and original research. New analytical results are immediately inside each theme. Earlier reports/figures/tables retain package names under `historical_evidence`; use the visual browser's search to find them. Source packages stay together so relative report/image links remain useful. Some inherited documents reference old model files or scripts outside the snapshots; their links are historical, not a claim that this clean runner executes them.

## Reproduce

From any directory:

```sh
python3 '/absolute/path/to/Clean folder/2. Scripts/run.py'
```

From this folder:

```sh
python3 '2. Scripts/run.py'
python3 '2. Scripts/run.py' --verify-only
python3 '2. Scripts/run.py' --reports-only
```

Use `--import-exploration` only to refresh the source/evidence snapshots intentionally. It is not required on ordinary reruns after collection. The clean analytical stages read their own copies and frozen checkpoints. The original Exploration files are preserved. Raw data never change during analysis, missing deaths never become zero, and generated models never read `.RData` or `.Rhistory`.

Requirements: Python 3.9+ standard library; R 4.6.1 was used locally, with dplyr, tidyr, splines, ggplot2, gridExtra, mgcv and MASS. Pandoc/Typst render reports; the runner discovers PATH versions or RStudio's bundled macOS tools. Rendering uses Georgia. Installed package versions are recorded in session logs. No package installation, upgrade or network access occurs during reproduction.

Count-model caches use data and fitting-code hashes. A changed fitting dependency invalidates those caches. All stages produce logs under `3. Outputs/01_data_description/logs/pipeline`; requested optional model failures remain in the model registry. Unresolved negative-binomial alternating optimizers are stopped rather than reported as successful fits. Rebuilding plots, table images and PDFs can take several minutes.

## Scientific scope

The rewritten core independently reconstructs source cells and the corrected regional panel, fits distributed-lag count models, computes paired demographic contrasts, explores minima/hinges and classification rules, rebuilds Census population allocations, fits annual county gradients, and produces genuine penalized GAM interaction surfaces. Primary response checkpoints, all 12 historical demographic contrasts and the county poverty checkpoint are independently checked.

Specialized historical regional architectures, Stage 3 latent-allocation models, large calibration experiments and morbidity branches are preserved with their reports and outputs. **They are not all rewritten or freshly refitted by this runner.** The result concordance makes that distinction explicit. Historical corrections and later undercoverage findings take precedence over an older report's stronger wording. The central numerical and interpretive checks do not certify missing source provenance, original catchments, statistical calibration or causal identification.

The project supports conditional aggregate associations and methodological comparisons. Marginal race/sex contrasts are not age-standardized susceptibility. Annual county poverty gradients are not daily poverty-by-temperature interactions. Population-days are allocation scenarios under unverified catchments, not unique people harmed or attributable deaths. See the reports for precise limitations and equations.

## Reading and editing scripts

`00_shared/setup.R` declares common paths, groups, dates and graphical conventions. `distributed_lag.R` documents the basis, fitting and paired covariance mathematics. Numbered theme scripts follow the order in `run.py`, using explicit intermediate objects and explanatory English section comments inspired by the user's ENSAE script organization. No personal machine path is required for the analytical stages.

Change analytical choices in the source scripts deliberately and rerun the full pipeline. Keep model IDs, footprints, controls, references, windows, calendars and uncertainty labels together when comparing CSV rows. Generated report text lives in `00_shared/build_reports.py`; edit that source instead of editing rendered HTML/PDF files that a rerun replaces.

## Git and possible release

The local `.gitignore` excludes raw/derived/reference data, extracted documents, generated evidence, models, logs and outputs. This supports a compact code-first commit without accidentally publishing supplied mortality extracts or gigabytes of historical snapshots. Source manifests and numerical references remain available locally for auditing. Generated HTML/PDF/PNG deliverables can be distributed separately after deciding which data and figure permissions apply.

No license is imposed, and no Git commit, push or public release has been performed. Choose a license and document the input-data access procedure before a public code release. The pipeline can run offline with the documented input/reference bundle; a code-only checkout does not contain those excluded inputs.
