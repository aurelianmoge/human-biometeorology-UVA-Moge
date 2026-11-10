# Temperature, DTR and mortality inequalities in Virginia, 2005–2020

Reproducible analysis pipeline for the Human Biometeorology project (UVA, 2026–2027).
Question: **do sex, race and age groups differ in how cold, heat and the diurnal temperature
range (DTR) affect their daily mortality, and do group-specific minimum-mortality temperatures
(Cao et al. 2026) change who is counted as "at risk"?**

Data: about 983,000 resident deaths from the Virginia Department of Health, in 12 climate regions,
each represented by one weather station (Davis et al. 2023). Methods: two-stage distributed lag
non-linear models, the published statewide design of DeGuzman et al. (2026), time-stratified
conditional Poisson regression, attributable burden, and the Cao et al. (2026) days-at-risk comparison.
See [`4. Reports/DATA_AND_METHODS.md`](4.%20Reports/DATA_AND_METHODS.md).

## Layout

| Folder | Content | In git |
|---|---|---|
| `1. Data/raw/` | 12 regional workbooks (VDH deaths + station weather), Census county estimates, GHCN-Daily donor stations | no (restricted / large) |
| `1. Data/reference/` | frozen reference panel, station and county metadata (Davis et al. 2025 roster) | no |
| `1. Data/derived/` | everything rebuilt by stage 01 | no |
| `2. Scripts/` | the pipeline (below) | **yes** |
| `3. Outputs/` | tables (`.csv`), figures (`.png`), fitted models (`.rds`), logs | no (rebuilt) |
| `4. Reports/` | Data & Methods document | yes |

## Running

From this folder (`Final folder`), with R ≥ 4.3 (`data.table`, `dlnm`, `mixmeta`, `fixest`, `ggplot2`, `dplyr`, `tidyr`) and Python 3 (standard library only):

```sh
Rscript "2. Scripts/run_all.R"                    # about 25 min on 8 cores
Rscript "2. Scripts/run_all.R" --with-simulation  # adds the 300-replicate null simulation (about 100 min)
Rscript "2. Scripts/run_all.R" --from 05          # restart from a stage
```

The input data are not public. The pipeline needs the `1. Data/raw` and `1. Data/reference` folders described above.

## Pipeline

| Stage | Script | What it does |
|---|---|---|
| 00 setup | `config.R`, `dlnm_tools.R` | Paths, constants, figure style; two-stage DLNM toolkit (stage 1 fits, `mixmeta` pooling, curve tests) |
| 01 data | `01_extract_workbooks.py` | Original Excel cells, corroborated date repairs |
| | `02_census_populations.py` | Census county × age × race × sex populations (July-year mapping, infant bin, Bedford) |
| | `03_assemble_panel.R` | 12 × 5,844 region-day calendar; duplicate dates set to missing; weather quality rules |
| | `04_fetch_ghcn.py` | Public NOAA GHCN-Daily Tmin/Tmax for the 12 stations and nearby donors |
| | `05_fill_weather.R` | Fills the 544 missing station-days (0.8%): same GHCN station, nearby station or other Virginia stations, chosen by masked validation |
| | `06_county_crosswalk.R` | County → region crosswalk (Davis 2025 roster + Davis 2023 rules) and regional populations |
| | `07_analysis_panel.R` | Region-day analysis panel |
| 02 descriptives | `01_descriptives.R` | Deaths, rates, standardised mortality ratios, regional climate, trends, seasonality, raw association |
| | `02_exposure_choice.R` | Quasi-AIC of Tmin vs Tmean vs Tmax vs afternoon temperature |
| 03 temperature | `01_two_stage_dlnm.R` | Main model: region × margin DLNM on Tmin, multivariate meta-analysis, MMT bootstrap |
| | `02_group_contrasts.R` | Whole-curve tests, ratios of RR, age composition (race), reference dependence, power |
| | `03_sensitivity_grid.R` | 14 one-factor-at-a-time specifications |
| 04 DTR | `01_two_stage_dtr.R` | DTR net of temperature level; seasonal, humidity, anomaly and lag variants; lag-response by sex |
| 05 group differences | `01_published_design.R` | DeGuzman et al. (2026) statewide design for every margin, Tmin and DTR, with direct tests |
| | `02_conditional_poisson.R` | Time-stratified conditional Poisson slopes and equality tests |
| | `03_attributable_burden.R` | Cold/heat-attributable deaths per 100,000; age-standardised race comparison |
| 06 days at risk | `01_days_at_risk.R` | Cao et al. (2026) own-curve vs overall-curve days and person-days at risk, with noise benchmark |
| 07 methods supplement | `01_null_simulation.R`, `02_real_vs_null.R` | What the methods report when no group differs; calibrated power |
| 08 figures | `01_main_figures.R` | Figures 1–8 |

Earlier exploratory work (Codex pipeline, Claude exploration) is archived outside this folder, in `../Exploration/` and `../Claude/exploration/`.
