"""Build thematic scientific reports and an offline visual results browser.

All numbers in the newly fitted result tables come from clean CSV exports.
Historical findings have explicit source links and are not counted as new fits.
"""
from pathlib import Path
from urllib.parse import quote
from collections import Counter
import csv
import html
import json
import os
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[2]
THEMES = ["01_data_description", "02_descriptive_statistics", "03_modelization", "04_original_research"]
TITLES = ["Data description and audit", "Descriptive statistics", "Temperature–mortality models and thresholds", "Inequalities and methodological research"]
REPORTS = ROOT / "4. Reports"


def read(relative):
    with (ROOT / relative).open(encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def number(value, digits=3):
    try:
        return f"{float(value):.{digits}f}"
    except (TypeError, ValueError):
        return "unavailable"


def markdown_table(rows, columns, limit=30):
    result = ["| " + " | ".join(label for field, label in columns) + " |", "| " + " | ".join("---" for _ in columns) + " |"]
    for row in rows[:limit]:
        values = []
        for field, label in columns:
            value = row.get(field, "")
            if field in {"RR", "Lower", "Upper", "SE", "Ratio", "Coverage", "Mean_deviance", "MAE", "Minimum", "First_upper_RR105", "Change_percentage_points"}:
                value = number(value)
            else:
                try:
                    numeric = float(value)
                    value = format(numeric, ".4g")
                except (TypeError, ValueError):
                    pass
            values.append(str(value).replace("|", ";").replace("_", " "))
        result.append("| " + " | ".join(values) + " |")
    return "\n".join(result) + "\n"


def url(path, parent):
    return quote(os.path.relpath(ROOT / path, parent), safe="/._-")


def link(path, label, parent):
    return f"[{label}]({url(path, parent)})"


def figure(theme, name, caption, parent):
    return f"\n![{caption}]({url('3. Outputs/' + theme + '/png/' + name + '.png', parent)})\n"


def historical(source, label, parent):
    manifest = read("1. Data/metadata/evidence_manifest.csv")
    entry = next((row for row in manifest if row["Source"] == source), None)
    return link(entry["Clean_path"], label, parent) if entry else f"{label} (source not collected: {source})"


def render(path):
    bundled = Path("/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64")
    pandoc = shutil.which("pandoc") or str(bundled / "pandoc")
    typst = shutil.which("typst") or str(bundled / "typst")
    if not Path(pandoc).exists():
        raise RuntimeError("Pandoc is required for HTML/PDF reports. Markdown reports have been written.")
    common = [pandoc, str(path), "--from=markdown+tex_math_dollars", "--standalone", "--toc", "--resource-path=" + str(path.parent)]
    subprocess.run(common + ["--embed-resources", "--css=" + str(REPORTS / "report.css"), "-o", str(path.with_suffix(".html"))], check=True, capture_output=True, text=True)
    subprocess.run(common + ["--pdf-engine=" + typst, "--pdf-engine-opt=--root", "--pdf-engine-opt=/", "-V", "margin-left:20mm", "-V", "margin-right:20mm", "-V", "mainfont=Georgia", "-V", "fontsize=10pt", "-V", "papersize=a4", "-o", str(path.with_suffix(".pdf"))], check=True, capture_output=True, text=True)


def write_report(relative, title, body):
    path = REPORTS / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    content = f'---\ntitle: "{title}"\nsubtitle: "Virginia 2005–2020 · Clean project · 28 October 2026"\nlang: en\n---\n\n' + body
    path.write_text(content, encoding="utf-8")
    print("Rendering", relative, flush=True)
    try:
        render(path)
    except subprocess.CalledProcessError as error:
        raise RuntimeError(f"Report rendering failed for {path}: {error.stderr}") from error


def main():
    REPORTS.mkdir(exist_ok=True)
    (REPORTS / "report.css").write_text("body{max-width:1050px;margin:40px auto;padding:0 25px;font:17px/1.6 Georgia,serif;color:#203343}h1,h2,h3{font-family:system-ui,sans-serif;line-height:1.2}h1{margin-top:2em}img{max-width:100%;height:auto}table{border-collapse:collapse;width:100%;font-size:13px;display:block;overflow:auto}td,th{padding:8px;border-bottom:1px solid #cbd5dc;text-align:left}th{background:#e8eff4}a{color:#196583}pre{white-space:pre-wrap;background:#eef3f5;padding:15px}blockquote{border-left:4px solid #b44d30;padding-left:20px}nav{background:#eef3f5;padding:20px}")
    evidence = read("1. Data/metadata/evidence_manifest.csv")
    status = read("3. Outputs/03_modelization/csv/model_status.csv")
    checks = read("1. Data/metadata/verification_checks.csv")
    counts = Counter(row["Status"] for row in status)
    estimates = read("3. Outputs/03_modelization/csv/response_estimates.csv")
    direct = read("3. Outputs/04_original_research/csv/direct_demographic_contrasts.csv")
    county = read("3. Outputs/04_original_research/csv/county_social_gradient_regressions.csv")
    minima = read("3. Outputs/03_modelization/csv/fitted_minima_and_crossings.csv")
    classification = read("3. Outputs/04_original_research/csv/pooled_vs_own_date_classification.csv")
    calibration = read("3. Outputs/04_original_research/csv/archived_calibration_recomputed.csv")
    validation = read("3. Outputs/03_modelization/csv/blocked_validation_scores.csv")
    gam = read("3. Outputs/04_original_research/csv/GAM_model_diagnostics.csv")

    parent = REPORTS / THEMES[0]
    body = r'''# What the data can support

The core observations are twelve supplied Excel workbooks, with daily weather and separate mortality margins over 1 January 2005–31 December 2020. There are 5,844 dates and 70,128 source rows before resolving duplicate dates. The observational unit is a workbook-region date. It is not an individual or a county-date. A region label is an assigned station/workbook label; the original residential catchment and possible overlap between extracts are not verified.

The clean scripts independently extract the original XML cells, reconstruct continuous regional calendars, and reproduce the earlier corrected panel. They never read a saved R workspace. Public county mortality and Census demographic populations are separate inputs. Annual county observations cannot be expanded into invented daily county deaths. Population age-by-race-by-sex cells do not supply crossed death cells.

# Source structure and variable dictionary

`1. Data/raw` contains byte-preserved copies of 52 usable source files. The unfinished ERA5 download remains in Exploration and is excluded from the input manifest. The source hashes and original paths are in `raw_manifest.csv`. An exhaustive inherited workbook dictionary records each field's original name, position, unit, missingness, range and candidate role in `1. Data/reference/metadata/variables.csv`. It describes source fields; it is not independent validation of their meaning.

| Quantity | Clean variable | Definition and interpretation |
| --- | --- | --- |
| Daily all-cause deaths | All | Supplied Mort column; nonnegative count, never interpolated |
| Race margins | Black, White | Recorded separate margins; residual includes other/unknown classifications |
| Recorded sex | Female, Male | Supplied categories; not a measure of gender identity |
| Age | Under50, Age50_74, Age75plus | Sums of supplied age bands; all bands retained in reconciliation |
| Causes | Respiratory, Circulatory | Supplied MDC.04/MDC.05 labels; exact extract coding remains conditional |
| Daily minimum/maximum | Tmin, Tmax | Same-row extrema in Celsius, after basic physical QC |
| Temperature midrange | Tmean_midrange | (Tmax + Tmin)/2; not a measured 24-hour mean |
| Daily range | DTR | Tmax − Tmin computed within station before averaging |
| Relative humidity | RH07 | Supplied morning humidity only when physically supported |
| Pollution | PM, Ozone | Local RAW anomalies; not measured concentration increments |

Repeated source labels and cause categories remain documented in the original dictionary. Only explicitly selected, interpretable variables enter the clean analytical panels. Afternoon observations, wind, pressure, pollution and demographic margins remain available for descriptions or source-bounded extensions.

# Dates, duplicate records and missing observations

Excel serial dates can be fractional. The reader preserves serials, source row numbers, original floor dates and year/month/day fields. Fourteen dates are corrected only when the corroborating internal checks agree: the target is one day later; the serial is near midnight; the chronological row implies that date; the target is absent; the previous date is duplicated; and adjacent year/month/day fields support the target. Fifteen unresolved date-field conflicts remain unresolved in the historical audit. No blanket one-day shift is applied.

For duplicate region-date records, mortality is missing rather than summed, averaged or set to zero. Available duplicate weather is averaged under the inherited continuity convention. All absent dates are then inserted. This is essential: deleting missing days before constructing lags would turn a nominal seven-day lag into a variable-duration lag.

Primary weather histories remain observed-only. Death counts are never filled. An observed zero is retained as zero. The inherited short-gap sensitivity and masked-weather validation are preserved in the evidence library; that sensitivity is not the primary new fit. Missingness is reported by variable, station and year. Missing counts include inserted absent dates as well as blank/invalid observations.

# Weather checks and exposure proxies

Inverted extrema and broad physically implausible extrema are excluded in derived data. Morning RH must lie in [0,100], have dew point no more than 0.5°C above synchronized temperature, and agree within 15 percentage points with the Magnus-formula RH reconstructed from morning temperature/dew point. This is a declared screening rule, not proof of instrument accuracy or synchronization.

The supplied PHF DTR column disagrees with its daily extrema on 4,236 rows. This is independently recovered from the original cells. New analyses always use Tmax−Tmin. Other source-unit disagreements remain in the audit; a consistency check cannot decide which conflicting measurement is correct. Warm-season zero or subzero readings are suspicious but do not automatically establish an error.

The temperature proxy is the equal-weight mean of twelve fixed stations. An aggregate day is missing if any required station is missing. Deaths are summed over the same twelve supplied extracts, conditional on unknown catchments. DTR models use eleven fixed stations with more than 1,000 physically supported RH07 observations; VJI is excluded. The station footprint is fixed before fitting. We never change weights according to available stations or daily mortality.

DTR is derived at station level first. With a fully matched fixed footprint, averaging station DTR equals mean Tmax minus mean Tmin. The explicit order still prevents accidentally mixing different station subsets. Temperature and DTR primary results use different footprints, controls and eligible dates, so their comparison does not isolate a biological metric effect.

# County-to-region correspondence

The 2025 article supplement documents seven article catchments, comprising 70 historical localities or 69 current county-equivalents. The inherited audit estimates 77.2% of Virginia's 2020 population in that published roster. Those assignments are distinct from the unknown memberships of the twelve supplied workbooks. Roanoke County is not silently added to the paper roster. Nearest-station assignments are exploratory geographic scenarios rather than verified residence linkages.

The complete crosswalk distinguishes Reference_region, Project_region, status, locality FIPS, boundary vintage and nearest-station sensitivity. Blank Project_region is an unknown assignment, not an excluded county. Consequently the clean acute models do not attach county poverty to aggregate daily deaths, invent regional population offsets, or estimate neighborhood vulnerability.

# Census years and demographic populations

The older upload uses AGEGRP99 for totals and AGEGRP0 for infants; infants must be added to the age bins. The later ALLDATA6 upload uses YEAR3–12 for July 2010–2019, YEAR13 for April 2020 and YEAR14 for July 2020. The clean code independently reconstructs 2,133 source county-year age-sum identities and 38,304 harmonized age×race×sex population cells. Historical Bedford city (51515) is combined with Bedford County (51019) for consistent county presentation. State totals are unaffected by that merge.

These are race-alone Census populations, including Hispanic residents. Compatibility with death race coding is not proven. Daily population values in allocation scenarios interpolate between July 1 annual estimates and hold endpoints constant. Those daily values are analytical assumptions, not official daily population estimates. Older population-dependent results retain their historical values and should not be mixed with the corrected clean allocation denominator.

# Identification and practical consequences

Without joint daily age×race×sex deaths, marginal racial and sex responses cannot be interpreted as age-standardized susceptibility. Without verified catchments, no local population mortality rate can be constructed. Without physical pollution concentrations and monitor metadata, a local anomaly effect is not a concentration inequality. Without verified observation clocks, the supplied afternoon field cannot support a precise local-standard-time exposure claim.

The missing inputs are specific: daily crossed mortality counts by residence geography, original extract-to-locality membership and overlap, cause/race coding definitions, weather clock/observation-window metadata, and monitor-level pollution metadata. These are limitations of the available design. No provider request is sent by this workflow.
'''
    body += figure(THEMES[0], "missingness_heatmap", "Missingness on complete regional calendars; blanks/invalid values remain missing.", parent)
    body += "\n# Files to inspect\n\n" + "\n\n".join([
        link("1. Data/metadata/raw_manifest.csv", "Immutable input manifest", parent),
        link("1. Data/reference/metadata/variables.csv", "Full inherited workbook dictionary", parent),
        link("1. Data/reference/metadata/county_region_crosswalk.csv", "Verified article roster versus unverified workbook mapping", parent),
        link("3. Outputs/01_data_description/csv/independent_panel_reproduction.csv", "Independent regional-panel comparison", parent),
        link("3. Outputs/01_data_description/csv/supplied_DTR_identity_audit.csv", "Independent DTR identity audit", parent)])
    body += "\n\nThe year mapping is checked against the [official Census layout](https://www2.census.gov/programs-surveys/popest/technical-documentation/file-layouts/2010-2020/cc-est2020-alldata6.pdf)."
    write_report(THEMES[0] + "/REPORT.md", TITLES[0], body)

    parent = REPORTS / THEMES[1]
    body = r'''# Reading descriptive results correctly

The plots describe observed data, before causal interpretation. A raw mortality–temperature curve mixes seasonal infections, long-term population trends, changing composition, reporting coverage and delayed weather histories. A dense cloud at moderate temperature can reflect the abundance of moderate days, not a precise health optimum. The regression chapters address some of these sources of structure but do not make the raw cloud causal.

# Regional weather and missingness

Regional summary tables report observed N, mean, standard deviation, P01/P10/median/P90/P99 for four weather metrics and ten outcome margins. Quantiles use observed finite rows. Different variables can therefore have different denominators. Use the annual missingness tables alongside any regional comparison. Geography affects both weather and the quality of observation support.

The seasonal profiles average available values within region-month. Temperature minima/maxima and midrange show the annual seasonal cycle; DTR describes the daily span rather than temperature level. A region with large DTR is not necessarily hotter. Temperature and range are mechanically related through the same extrema, which motivates controlling temperature when studying DTR.

# Mortality over time

Annual death plots sum observed region-date counts. They are not annual mortality rates, and totals with incomplete observation calendars are not directly comparable with complete-year totals. The supporting CSV gives the number of observed all-cause dates for each region-year. The aggregates may change with population growth, age structure, pandemic mortality, coding or catchment composition. Excluding 2020 later tests one consequential source of change.

Age, race and sex are separate margins. The same death can occur in an age group, a race group and a sex group, so summing across all displayed demographic categories double-counts deaths. Female+Male need not equal All exactly because a residual recorded category can remain. Black+White does not exhaust recorded race. The data chapter reconciles these margins without assigning residual deaths to a preferred category.

# Mortality versus Tmin, Tmax, midrange and DTR

The four cloud plots show daily supplied aggregate deaths against the equal-station exposure values. The smooth line is a descriptive local regression, with no lag or demographic adjustment and no inferential ribbon. Points with missing outcome or missing exposure are omitted. A U-shape may be seasonal; a weak hot-end trend may reflect sparse observations. Neither the visual minimum nor the line's slope is an estimated causal threshold.

The correlation heatmap pools observed regional rows and uses pairwise-complete correlations. Every cell can use a different set of days. Between-region differences can also induce associations even if within-region changes differ. Correlation is useful for identifying collinearity and choosing interpretable model specifications, not for deciding that correlated weather variables are independent causes. Tmin, Tmax and DTR satisfy an exact arithmetic identity and must not enter together as three unrestricted independent predictors.

# Demographic population descriptions

Corrected Census trajectories describe annual statewide resident populations for recorded sex, race-alone categories and broad age groups. They supply context for count-versus-rate interpretation. They are not validated denominators for unknown workbook catchments. Separate group trajectories overlap across dimensions. A racial mortality count difference may partly reflect different population sizes and age structures; neither counts nor population-only adjustment recovers the missing joint death cells.

# County socioeconomic and geographic context

The annual county panel contains 396 rows for 132 localities in 2010, 2013 and 2016, with public crude and age-adjusted mortality rates, poverty, household income, racial composition, older-age share, population and coordinates. These are annual ecological observations. A poverty–mortality cloud indicates a between-place gradient, not a daily poverty-specific temperature response.

Public age-adjusted rates reduce confounding by broad age composition in the annual outcome, but remain standardized indices. They cannot be fed into a Poisson model as if they were observed death counts. Repeated observations of the same locality are dependent, motivating county-clustered covariance in the original-research chapter. Longitude, latitude and independent-city status are descriptive geographic controls, not a complete spatial causal model.

# How to use the tables

Read the observed N first, then the central tendency and spread, and compare equivalent units and footprints. Regional summary CSVs preserve full precision. Selected compact analytical tables have paginated PNG companions. Large daily panels and prediction grids remain CSV because thousands of rows are more useful as plot data than as table images. This distinction is documented in the table-image register.
'''
    for name, caption in [("seasonal_weather", "Monthly regional weather profiles in degrees Celsius."), ("annual_deaths", "Annual deaths on observed days; inspect coverage before comparing totals."), ("mortality_temperature_clouds", "Unadjusted mortality–weather cloud plots."), ("correlation_heatmap", "Pairwise-complete pooled regional correlations."), ("population_demographic_trajectories", "Corrected July-year demographic population trajectories."), ("county_poverty_mortality_cloud", "Annual county poverty and age-adjusted mortality.")]:
        body += figure(THEMES[1], name, caption, parent)
    write_report(THEMES[1] + "/REPORT.md", TITLES[1], body)

    parent = REPORTS / THEMES[2]
    primary = [row for row in estimates if row["Outcome"] == "All" and row["Specification"] == "primary" and row["Window"] != "Days_3_7"]
    primary_minima = [row for row in minima if row["Outcome"] == "All"]
    body = f'''# Model program and execution\n\nThe clean registry requests 56 count-model fits: primary analyses across four exposures and nine outcome margins, temporal-adjustment and pandemic sensitivities, negative-binomial alternatives, sequential DTR controls and a temperature–DTR interaction. Execution status: {dict(counts)}. Every requested model remains in the registry, including stopped fits. A convergence warning is not silently treated as a successful sensitivity model.\n\n'''
    body += r'''# Why use a count model with distributed lags?

Mortality is a nonnegative count. A quasi-Poisson log model estimates a positive conditional mean while allowing variance greater than the mean. A one-degree temperature change today may have a different association from the same change two weeks earlier. A distributed-lag nonlinear model represents both the exposure shape and the lag pattern:

$$\log\mu_t=\alpha+s(t)+\gamma_{\mathrm{weekday}(t)}+\sum_{\ell=0}^{21}f(x_{t-\ell},\ell)+z_t^\top\eta.$$

Here $s(t)$ has 112 total natural-spline degrees of freedom over sixteen years (7 per year), rather than 112 per year. Calendar adjustment and weekday are nuisance controls. DTR models add the distributed lag of midrange temperature and a three-degree-of-freedom morning RH spline. Temperature models use the common twelve-station proxy. DTR uses its fixed eleven-station proxy. There is no population offset because the workbook catchments are unverified.

The clean engine constructs a quadratic B-spline in exposure with two equally spaced interior knots across the observed exposure range. The lag basis is a natural spline with an intercept and two log-spaced interior knots from day 1 to day 21. Their convolution is built on the complete calendar before selecting fit rows. The code checks several matrix rows against explicit convolution. Exposure references are weather-only medians, common across groups and sensitivities. No mortality outcome selects the reference or spline knots.

The underlying cross-basis method follows [Gasparrini, Armstrong and Kenward (2010)](https://pmc.ncbi.nlm.nih.gov/articles/PMC2998707/). This implementation uses base-R splines because dlnm is not installed; it reproduces the project's audited mathematical specification rather than claiming exact source-paper software equivalence.

# Relative responses, windows and controls

A prediction contrasts a specified exposure with the common median and sums the requested lag contribution. Acute means lag 0–2; cumulative means lag 0–21; days 3–7 are a separate timing summary. A cumulative RR is not the actual mortality prediction for a specific observed day with its mixed weather history. It is a response-profile contrast under the model.

Temperature anchors are P10 and P90. DTR anchors are P95 and P99. These are metric-specific physical values, not equivalent doses. A cold association consistent across models is different from proof of causality. A cumulative RR below one at a warm anchor is not evidence that heat is harmless or protective; temporal adjustment, selection and lag aggregation can contribute.

The DTR sequence holds the eligible calendar constant while changing temporal-only, same-day temperature, distributed-temperature and humidity controls. This separates control effects from simple sample changes. It does not solve residual confounding or measurement error. The product-interaction model estimates a same-day DTR×midrange term in addition to the lagged main effects; its DTR curve evaluated with interaction set to zero refers to median midrange temperature.

# Primary newly reproduced estimates

The following are model-conditional pointwise summaries. Nominal uncertainty is retained for scientific continuity; later calibration findings rule out interpreting every nominal interval as a validated 95% procedure.

'''
    body += markdown_table(primary, [("Exposure", "Metric"), ("Window", "Lag window"), ("Exposure_value", "°C"), ("RR", "RR"), ("Lower", "Lower"), ("Upper", "Upper")])
    body += figure(THEMES[2], "cumulative_response_curves", "Cumulative exposure-response curves relative to weather medians. Ribbons are nominal pointwise intervals.", parent)
    body += figure(THEMES[2], "lag_response_profiles", "Single-lag profiles reveal timing without equating acute and cumulative effects.", parent)
    body += r'''# Minimum mortality temperature and threshold alternatives

A minimum is the grid point minimizing a fitted cumulative log-response on P01–P99. It is a fitted optimum, not the temperature where harm begins. Boundary minima and broad/competing troughs make a single number particularly fragile. The clean workflow retains the full curve, a boundary flag and conditional coefficient-draw minimum locations for all-cause models. Draw ranges are not calibrated minimum confidence intervals.

An RR1.05 crossing asks where the response exceeds the median-reference response by 5%. A missing crossing means not reached within the domain, never zero degrees or evidence of safety. Classification uses the entire qualifying set, including disconnected parts, not every value beyond a first crossing. A percentile cutoff is a weather-rarity rule; an MMT split is an optimum-location rule; a support set is an evidence-and-precision rule. They need not identify the same dates.

'''
    body += markdown_table(primary_minima, [("Exposure", "Metric"), ("Minimum", "Grid minimum °C"), ("Boundary_minimum", "Boundary?"), ("First_upper_RR105", "First upper RR1.05 °C")])
    body += r'''# Segmented threshold exploration

The clean hinge model is continuous with one slope below a candidate threshold and another above it, adjusted for the same temporal and weather nuisance terms. The candidate profile spans P15–P85 in 2.5-percentile steps. It uses the cumulative-model eligibility mask for consistency, but the hinge term is same-day. The vertical profile is Poisson deviance above the best candidate; quasi-Poisson likelihood AIC is not invented.

The best candidate is an in-sample descriptive optimum. A naive interval conditional on that selected cutoff ignores selection. No new breakpoint p-value or biological onset claim is reported. The older segmented-search, extended-support and held-year validation packages are preserved separately, including their null findings. Their larger discovery families are not replaced with a smaller family after inspecting results.

'''
    body += figure(THEMES[2], "segmented_threshold_profiles", "Exploratory same-day breakpoint profiles; not selection-adjusted confidence sets.", parent)
    body += r'''# Article replication and source compatibility

The seven-region DTR benchmark follows the archived Davis et al. 2025 article and supplement (DOI 10.1007/s00484-025-02850-6). The clean fixed-rule rerun applies location-specific physical DTR cutoffs to audited workbook extrema. It counts exposure-event dates separately from the union of forward-lag outcome dates. IAD includes a finite DTR interval; ROA has five DTR-dependent lag bands. Each outcome date is counted once even when multiple events lead to it.

Missing-event lower and upper bounds count unknown weather classifications. They are not sampling confidence intervals and do not include uncertainty in author-estimated thresholds. Published counts remain alongside the workbook counts; disagreement is not tuned away. The inherited NOAA benchmark, matched paper-basis/current-region approximations, complete curve atlas and all 21 composition odds ratios are retained with source and assumption labels. Exact reproduction remains limited by author source snapshot, catchments, clocks and basis conventions.

The earlier Tmin consolidation and article replication use distinct anchors, bases and calendars. Their figures are historical evidence rather than extra independent replications of this clean fit. A number such as RR1.088 at the historical cold Tmin anchor must not be substituted for a new Tmax percentile response.

# Sensitivity, diagnostics and validation

Calendar splines with 48/80/112 total degrees of freedom show sensitivity to seasonal control. The no-2020 models use 105 total temporal degrees of freedom for the fifteen-year fitting sample. Exposure bases and references remain fixed. A negative-binomial sensitivity estimates extra dispersion parametrically; unresolved alternating-optimizer fits are stopped, retained in the registry and excluded from the interpretation.

Residual diagnostics correlate pairs separated by real calendar days. HAC bandwidths 14/28/56 are exported. HAC addresses a specific covariance approximation; it does not validate model choice or remove source error. Later known-truth experiments showed undercoverage for some confidence-region and classification procedures. New nominal results remain exploratory diagnostics in light of that evidence.

Four 28-day interior validation blocks have 42-day buffers. Temporal spline knots are rebuilt from training dates; weather bases use only exposure data. A fold with fewer than five usable test dates or deficient rank is explicitly omitted. Daily held-out predictions and Poisson deviance contributions are exported. Scores below weight each evaluated day equally, not each fold equally. These few internal blocks assess prediction on observed support, not external validation or causal performance.

'''
    body += markdown_table(validation, [("Exposure", "Metric"), ("Candidate", "Candidate"), ("Test_observations", "Held-out rows"), ("Mean_deviance", "Mean deviance"), ("MAE", "MAE")])
    body += "\n# Historical benchmark entry points\n\n" + historical("article_replication_benchmark_v1/reports/CLEAR_RESULTS.md", "Seven-region benchmark and all conventional intervals", parent) + "\n\n" + historical("3. Results/Segmented threshold validation/README.md", "Historical calibrated segmented-threshold validation", parent)
    write_report(THEMES[2] + "/REPORT.md", TITLES[2], body)

    parent = REPORTS / THEMES[3]
    body = r'''# Research question and estimands

The project asks how social inequalities affect temperature–mortality relationships, and how methodological choices affect apparent vulnerability. The available design can compare marginal demographic responses, describe annual county social gradients, and study instability in threshold/classification summaries. It cannot identify individual susceptibility, crossed-demographic effects or county-specific acute poverty interactions from the existing margins.

Five questions must remain separate: does mortality change with temperature; do two groups have different response curves; do fitted thresholds differ; do rules classify different dates; and do those classified dates imply different allocated population totals? A positive answer to one is not automatically a positive answer to the others.

# Direct demographic comparisons

For groups g and h, the relevant comparison is RRR=RRg/RRh at a common physical exposure, reference, lag window and fitting calendar. The standard error uses the difference of their shared-date influence contributions, retaining covariance from overlapping deaths and common weather. Comparing whether one group's confidence interval crosses one and another's does not is not a test of a group difference.

The clean implementation independently reproduces all twelve historical direct Tmax/DTR comparisons, including coefficients, standard errors and Holm12 arithmetic. None passes the historical Holm12 family at 0.05. Those historical decisions are model-conditional and retrospective; a numerical reproduction does not fix later calibration concerns. An uncertain comparison is not equivalence and does not establish absence of inequality. Age, race and recorded sex are separate margins, so these are not mutually adjusted demographic effects.

'''
    body += markdown_table(direct, [("Comparison", "Groups"), ("Exposure", "Metric"), ("Window", "Window"), ("RR", "RRR"), ("Lower", "Lower"), ("Upper", "Upper"), ("Historical_Holm12_recalculated", "Historical Holm12")])
    body += figure(THEMES[3], "demographic_contrast_forest", "Direct demographic ratios with nominal pointwise intervals.", parent)
    for exposure in ["Tmax", "DTR"]:
        body += figure(THEMES[3], "subgroup_curves_" + exposure, exposure + " subgroup response profiles; direct contrasts, rather than curve appearance, assess differences.", parent)
    body += r'''# Calendar, geography and model architecture

The later seven-region extension compared an aggregate series, a region-date panel with shared response, and harmonized region-specific responses. Its native seven-region fit calendar has 3,330 dates; the primary-eleven/seven-source intersection has 2,093; the common regional architecture comparison has 1,899. They are different estimands and cannot share denominators. Equal-region standardization averages log-responses, then exponentiates; it is not population-weighted statewide standardization.

The clean evidence library retains every regional estimate, leave-one-region-out result, covariance export and classification table. These specialized regional architectures are historical fitted evidence, not refitted by the new core workflow. Formal geographic p-values were withheld in that package because a valid coherent regional generator had not been established and prior covariance calibration was weak. The clean reports retain that decision.

The matched-source follow-up attenuated the nominal low-DTR Female/Male differences that appeared on native seven-region calendars. The common-breakpoint prediction example also showed that an arithmetic mean of fold scores differs from a per-subgroup-day weighted score when denominators differ. Both are methodological results: calendar selection and aggregation can change apparent evidence without a new biological finding.

# Pooled versus own-group classifications

Each rule evaluates the entire modeled qualifying set on the same observed assessment dates in P01–P99. Assessment requires observed current exposure and domain support, whereas model fitting additionally requires complete lag histories and outcomes. It is therefore legitimate for their counts to differ, provided their definitions are explicit.

The clean classification grid contains four metrics, six demographic groups, two lag windows and RR cutoffs 1.02/1.05/1.10. For every comparison it exports pooled, own, own-only, pooled-only, both, neither, signed net change and the assessed denominator. Positive net change means more dates under the own-group rule. It does not establish underestimated actual health risk under the pooled rule. Gross disagreement and signed net change answer different questions.

The clean rules in this chapter include both sides of the median whenever RR exceeds the cutoff. The older population-risk headlines separately restricted lower/upper branches, so their numbers should not be treated as a discrepancy without matching the rule and assessment mask. That older primary upper-DTR cumulative RR1.05 analysis changed the allocated Black-minus-White fraction gap by about +9.74 percentage points and Female-minus-Male by −13.05 points, with large conditional uncertainty. These are historical allocated fractions, not observed death inequalities.

'''
    headline = [row for row in classification if row["Rule"] == "RR1.05" and row["Exposure"] == "DTR" and row["Window"] == "Cumulative_0_21"]
    body += markdown_table(headline, [("Group", "Group"), ("Assessed", "Assessed"), ("Pooled_dates", "Pooled dates"), ("Own_dates", "Own dates"), ("Own_only", "Own only"), ("Pooled_only", "Pooled only"), ("Change_percentage_points", "Change pp")])
    body += figure(THEMES[3], "classification_change_heatmap", "Signed changes under RR1.05 using complete qualifying sets, not upper-branch-only rules.", parent)
    body += r'''# Population-day allocation scenarios

For a group, H is the sum of its interpolated Census population on flagged dates, D is the sum on all assessed dates, and F=H/D. H repeats residents over days; it is not unique people, observed deaths or attributable burden. The clean workflow rebuilds July-year populations before calculating these quantities. All H/F results are full-Virginia allocation scenarios under a shared station proxy and unverified workbook catchments.

Comparing own and pooled rules within the same group uses the same population weights and denominator. Comparing groups can retain an effect of different population trajectories. With common constant weights, the pooled rule cancels from the signed difference in group gap changes. With different time-varying weights, cancellation need not be exact. Sex, race and age allocations overlap and must not be summed across dimensions.

# Why uncertainty in classifications is difficult

Thresholding is discontinuous. A small shift in a smooth response can move many nearby dates across RR1.05 together; an entire qualifying component can disappear. Gaussian coefficient draws preserve specified model dependence, but the resulting percentile ranges are not automatically calibrated confidence intervals for a selected threshold, day count or signed gap.

The historical zero-endpoint audit found exact zero mass, not merely rounded small effects. Some saved draws were negative even when a displayed lower endpoint was zero. Reconstructed date flags did not exactly reproduce the earlier scalar draw series, so a target-specific explanation of zero mass could not be established from those reconstructions. The later reliability package saved a new, explicitly distinct coefficient-and-membership series.

That reliability study reports uncertainty in both the total number and identity of flagged dates. In the primary Black acute example, 92 dates were flagged at the fitted curve, the conditional count range was 0–249, and 352 dates changed membership at least once across the finite draw set. This is predominantly uncertainty in cutoff location and totals, rather than a precise total hiding extensive offsetting swaps.

The clean script independently recomputes the following coverage frequencies from the archived simulation replicate rows, rather than copying their printed summaries. These are known-truth synthetic designs, not additional Virginia estimates.

'''
    body += markdown_table(calibration, [("Scenario", "Scenario"), ("Successful", "Refits"), ("Coverage", "Coverage"), ("Excludes_zero", "Excludes zero"), ("Mean_width", "Mean width")])
    body += r'''
Under the equal-response scenario the nominal central 95% range covers the true gap in 88.5% of datasets. That undercoverage is a reason to keep the empirical ranges descriptive. The later four-scenario joint-region audit found very broad projected day-gap bounds with high gap coverage, while true coefficients lay in the nominal joint region only 65–77.5% of the time. Broad projection can cover a derived gap even when the underlying confidence-region claim fails. Neither outcome validates every regional or weighted classification target.

No new finite-family correction can repair invalid component tests. Holm controls family-wise error when its input tests are valid; it cannot cure covariance undercoverage, source error, retrospective selection or nonidentification. The old Holm12/Holm24/Holm21 families remain distinct, and the classification comparisons do not receive invented Gaussian p-values.

# Competing minima as a methodological result

The consolidated follow-up found two shallow female DTR troughs. The low/high minima were approximately 6.556/17.216°C on Workbook7 and 6.519/17.419°C on NOAA7. Tiny changes in their relative depths selected opposite branches. On 3,247 common assessment dates the above-minimum rule classified roughly 124 versus 2,977 dates when the fitted rule changed. This is a large classification reversal from an unstable optimum, not evidence that a physiological harm threshold shifted by ten degrees. Preserve response curves and all local minima when a scalar MMT drives decisions.

# Annual county poverty and racial-composition gradients

The clean regressions model log annual mortality rates, with year indicators and county-clustered covariance. Poverty is scaled per +5 percentage points and Black share per +10 points. Age-adjusted rates are standardized indices and are modeled with least squares on the log scale, not as Poisson counts. Geography specifications add longitude, latitude and independent-city status; these are sensitivity controls rather than full spatial identification.

'''
    body += markdown_table(county, [("Outcome", "Annual outcome"), ("Specification", "Specification"), ("Term", "Predictor"), ("Ratio", "Rate index ratio"), ("Lower", "Lower"), ("Upper", "Upper")])
    body += figure(THEMES[3], "county_regression_forest", "County-clustered annual ecological gradients.", parent)
    body += r'''
The poverty-only age-adjusted model reproduces an index ratio about 1.101 per +5 pp poverty (approximately 1.076–1.126). This is a baseline annual county association. It does not estimate whether poverty amplifies a daily temperature effect, and it does not measure the effect of an individual's poverty. Conditional racial-composition coefficients depend on the poverty/geography specification; they are not individual racial effects.

# Nonlinear temperature–DTR interactions

Two genuine penalized GAMs are newly fitted with REML: additive weather smooths, and the same terms plus a tensor interaction between midrange temperature and DTR. Both use the same eleven-station complete-case calendar, weekday, morning RH and a calendar smooth with basis dimension 112. A basis dimension is not the fitted effective degrees of freedom. The GAM is same-day; its surface cannot be compared numerically with a cumulative DLNM surface as if they were the same target.

The tensor term allows the DTR response to vary with temperature while retaining separate main smooths. Surfaces hold other covariates fixed and reference joint weather medians. Grid cells farther than 0.15 standardized units from an observed temperature–DTR pair are masked. This support rule is descriptive and not a confidence statement. Large concurvity warns that smooth terms can explain similar variation. Approximate mgcv smooth p-values are retained only as diagnostics because serial dependence is not corrected in those tests. These GAMs do not establish interaction causality or identify poverty modification.

'''
    body += markdown_table(gam, [("Model", "GAM"), ("Observations", "Days"), ("Converged", "Converged"), ("Deviance_explained", "Deviance explained"), ("Maximum_concurvity", "Maximum concurvity")])
    body += figure(THEMES[3], "GAM_temperature_DTR_surfaces", "Same-day additive and tensor-interaction surfaces; unsupported weather combinations are masked.", parent)
    body += figure(THEMES[3], "specification_sensitivity_forest", "Specification sensitivity across all primary weather metrics.", parent)
    body += figure(THEMES[3], "cause_specific_forest", "Respiratory and circulatory cumulative associations, without claiming direct cause differences from separate intervals.", parent)
    body += figure(THEMES[3], "lag_effect_heatmap", "Lag-specific effect magnitudes; fill colors are not significance decisions.", parent)
    body += figure(THEMES[3], "residual_autocorrelation", "Remaining temporal dependence on the real calendar.", parent)
    body += r'''# Composition adjustment and partial identification

The historical Stage 3 branch fits a composition-weighted marginal model (S1) and derives latent-allocation bounds (S3), with required simulation validation. Its observed-only annual-weight cumulative racial multiplier is approximately 0.975, with a model-conditional interval around 0.914–1.040. S1 assumes a common within-age response multiplier and compatible public demographic weights. It is assumption-bound adjustment, not identified age standardization.

The S3 observed-calendar cumulative identification region is approximately [0.745,1.274] and includes one. This describes feasible latent allocations conditional on fitted margins and compatible positive populations. It is not a causal interval, and its coefficient-draw envelope is not a validated confidence interval. The source decision note explicitly postpones transfer to the full nonlinear DLNM. The clean workflow preserves that branch and its optimizer failures, corrected Census mapping, LP checks and simulations; it does not quietly promote S1 latent cells into observed deaths.

# Other completed and unavailable branches

Historical evidence also includes pollution-support failures and local-anomaly sensitivities, between-day temperature variability, morbidity visits, regional cold timing, calibrated threshold searches, response-metric concordance, exposure/person-day accounting and recovery pilots. They are linked in the evidence review rather than mixed into a pooled significance narrative. Morbidity uses a different reporting/calendar process; it is not a direct replication of death models. An unavailable pollutant concentration, daily county outcome or crossed-demographic comparison remains unavailable, not an estimated zero effect.

A logistic model is not added solely to increase model count: these aggregated deaths do not provide an individual daily binary risk denominator. Historical conditional composition odds ratios are retained under their own exchangeability assumptions and are not relabeled as mortality risk ratios.

# What is learned and what remains unresolved

The strongest defensible synthesis is a methodological case study with a conditional delayed-cold association and a separate annual poverty gradient. Marginal group comparisons and classified-day changes are informative within their stated definitions, but do not establish hidden early racial/sex vulnerability or unequal causal mortality burden. The most consequential next data improvement is crossed daily deaths with verified residential exposure assignment; more elaborate models cannot substitute for those observations.
'''
    body += "\n\n" + historical("calendar_geography_extension_v1/reports/Technical_Appendix.md", "Latest regional-calendar technical appendix", parent) + " · " + historical("biometeorology_consolidated_followup_v1/reports/consolidated_followup.md", "Consolidated competing-minimum and calibration audit", parent) + " · " + historical("population_risk_comparison_v1/method_reliability/reports/REPORT.md", "Classification reliability report", parent)
    write_report(THEMES[3] + "/REPORT.md", TITLES[3], body)

    # A compact general report refers readers to the detailed thematic explanations.
    parent = REPORTS
    body = f'''# Project guide\n\nThis clean project separates immutable data, readable analysis scripts, outputs and reports into four folders. It contains a newly written core workflow plus {len(evidence):,} traceable historical evidence files. The historical packages are scientific snapshots, not additional independent studies. A report retained from an older package can have been qualified by a later audit; the evidence review records that progression.\n\nThere are 56 registered clean count-model requests ({dict(counts)}), two penalized GAM fits, two segmented profiles, complete demographic/date classifications, corrected Census population allocations and annual county regressions. {sum(row['Passed'] == 'TRUE' for row in checks)} of {len(checks)} registered numerical checks pass. The clean core independently reconstructs source cells and the corrected regional panel; it reproduces primary response checkpoints, all twelve direct demographic contrasts and the county poverty checkpoint. Reproduction of arithmetic and fitted estimates does not verify missing catchments or inferential calibration.\n\n'''
    body += r'''# Scientific synthesis

The data support conditional delayed-cold associations under the declared aggregate count models and a separate annual county poverty–mortality gradient. The primary direct demographic comparisons remain uncertain after their historical correction family. This is not equivalence and does not establish absence of social inequality.

Threshold and classification analyses expose methodological sensitivity. A curve minimum, RR1.05 crossing, statistical-support region and weather percentile describe different quantities. Small changes between competing troughs can substantially reverse above-minimum classifications. Pooled versus own-group rules can therefore allocate very different day totals even when a direct demographic response contrast is uncertain.

Later known-truth work found uncertainty undercoverage. Conditional coefficient-draw ranges remain descriptive, and some regional formal tests are withheld. A nominal interval, a simultaneous band, a conditional draw range and a partial-identification region must remain separately labeled.

The eleven-station DTR models and twelve-station temperature models differ in footprint, controls and available dates. Shared-weather classifications cannot identify differential personal thermal exposure. Corrected Census population allocations describe repeated population-days under a scenario, not unique people harmed or attributable deaths.

Original county memberships and possible overlap between workbook catchments remain unverified. Daily age×race×sex death cells are missing. Consequently marginal racial/sex responses are not age-standardized susceptibility, and annual county poverty gradients are not acute poverty modification. These design limits determine which claims are possible.

# Reading sequence

'''
    for theme, title in zip(THEMES, TITLES):
        body += "- " + link("4. Reports/" + theme + "/REPORT.html", title + " — HTML", parent) + " · " + link("4. Reports/" + theme + "/REPORT.pdf", "PDF", parent) + "\n"
    body += "\n" + link("4. Reports/MATHEMATICAL_APPENDIX.html", "Mathematical appendix", parent) + " · " + link("4. Reports/HISTORICAL_EVIDENCE_REVIEW.html", "Historical evidence review and result concordance", parent) + " · " + link("3. Outputs/index.html", "Visual output browser", parent)
    body += figure(THEMES[2], "cumulative_response_curves", "Clean cumulative all-cause response profiles; model-conditional nominal intervals.", parent)
    body += figure(THEMES[3], "demographic_contrast_forest", "Direct demographic comparisons, not comparisons of separate significance labels.", parent)
    body += figure(THEMES[3], "classification_change_heatmap", "How rule pooling changes date classification, without claiming identified health burden.", parent)
    body += r'''# Reproduction and scope

From the clean root run `python3 '2. Scripts/run.py'`. The importer is needed only for first assembly or an intentional refresh of historical evidence. Later runs use local frozen inputs. Successful count fits are cached by input/code hashes; failed optional negative-binomial optimizers remain explicitly stopped. Per-stage logs and software versions are retained. The pipeline installs nothing and sends no external messages.

The new core consolidates central data handling, count-response, demographic, threshold, GAM and county work. Specialized historical regional architectures, Stage 3 latent-allocation models, large calibration experiments and morbidity branches are retained with their original outputs and reports, not represented as freshly refitted or fully rewritten by this runner. The reproduction-status register distinguishes these scopes explicitly. This prevents an attractive folder reorganization from overstating scientific reproduction.

The scripts use numbered stages, named intermediate objects, section comments and English explanations inspired by the user's ENSAE script style. Git exclusions keep raw/derived datasets and historical evidence out of a normal code commit; generated reports and visual outputs can be reviewed separately for publication. No commit, push or publication is performed.
'''
    write_report("GENERAL_REPORT.md", "Virginia biometeorology: clean project and scientific synthesis", body)

    build_math()
    build_history(evidence)
    build_browser(evidence)


def build_math():
    body = r'''# From counts to relative responses

Let $Y_{gt}$ denote the observed supplied death count for group g on date t. The quasi-Poisson mean model has $E(Y_{gt}\mid X)=\mu_{gt}$ and working variance $\phi_g\mu_{gt}$. The log link guarantees positive fitted means. The count model omits an offset because the original catchment population is unknown. Its exposure ratios are conditional count-response ratios, not calibrated individual daily risks.

A population offset would add $\log P_{gt}$ and model a rate only if P describes the same catchment, demographic definition and dates as Y. Substituting statewide population for unverified regional catchments does not satisfy that requirement.

# Cross-basis construction

Let $b_j(x)$ be exposure-basis functions and $c_k(\ell)$ lag-basis functions. Define

$$Z_{t,jk}=\sum_{\ell=0}^{L}b_j(x_{t-\ell})c_k(\ell),\qquad L=21.$$

The linear predictor contains $\sum_{j,k}\beta_{jk}Z_{t,jk}$. This separates nonlinear exposure shape from smooth lag variation. All histories are constructed on the real daily calendar. If a required exposure is missing, the corresponding basis row is missing; dropping missing dates first would change the meaning of lag.

For reference exposure $x_0$ and lag window W, the contrast vector is

$$a_{jk}(x,W)=[b_j(x)-b_j(x_0)]\sum_{\ell\in W}c_k(\ell).$$

Thus $\log RR_g(x,W)=a(x,W)^\top\widehat\beta_g$. At $x=x_0$, every element of a is zero, so RR is exactly one. Summing single-lag contrast vectors must equal the cumulative contrast vector; these are useful implementation checks. The model profile is a standardized exposure-history contrast, not the actual predicted count on a day with mixed exposures.

# Nuisance adjustment and temporal covariance

The clean implementation partials target columns out of the weighted nuisance design. With design X, target residualized matrix R and working weights W, the target bread is $(R^\top R)^{-1}$. Target influence contributions combine the response score, residualized target and bread. Ineligible calendar dates retain zero rows in the influence matrix, not compressed positions.

For full-calendar influence rows $u_t$ and bandwidth q, the Bartlett HAC covariance is

$$\widehat V=\sum_tu_tu_t^\top+\sum_{h=1}^{q}\left(1-\frac{h}{q+1}\right)\sum_{t=h+1}^{T}(u_tu_{t-h}^\top+u_{t-h}u_t^\top).$$

The primary q is 28 days, with 14/56 sensitivities. A nominal pointwise standard error is $\sqrt{a^\top\widehat Va}$, and exponentiating $\log RR\pm1.96SE$ yields the conventional interval. The numerical check that a covariance is symmetric positive semidefinite is necessary for computation; it is not a coverage validation. The project's later simulations show why that distinction matters.

# Direct group contrasts and overlapping outcomes

The log ratio of relative responses is

$$d_{gh}=a_g^\top\widehat\beta_g-a_h^\top\widehat\beta_h.$$

Construct date-aligned scalar influence $v_t=u_{gt}^\top a_g-u_{ht}^\top a_h$, then estimate its HAC variance. This automatically includes the covariance term in $V(d)=V(g)+V(h)-2\operatorname{Cov}(g,h)$. Treating the group models as independent discards that term. The clean paired function also requires identical eligible masks and dates, avoiding accidental contrasts between different fitted populations or calendars.

Different percentile temperatures are allowed only when the target is explicitly a paired metric-specific percentile comparison. Such a contrast is not the effect of changing only day versus night temperature at a common physical dose.

# Multiplicity and its limits

For sorted p-values $p_{(1)}\le\cdots\le p_{(m)}$, Holm's adjusted values are the monotone capped sequence

$$\widetilde p_{(i)}=\min\{1,\max_{j\le i}(m-j+1)p_{(j)}\}.$$

Holm controls the chance of any false rejection in the declared family when its component p-values are valid. The family m=12 is retained when showing a subset of the old demographic contrasts. BH targets a different error criterion, the expected false-discovery proportion under its validity conditions. Neither correction repairs invalid component tests, unmeasured confounding, source mismatch or identification failure. The clean reruns preserve nominal p-values and old correction arithmetic for audit, without treating them as new confirmed discoveries.

See the [official R p.adjust documentation](https://stat.ethz.ch/R-manual/R-patched/library/stats/html/p.adjust.html) for the procedures and assumptions. The general project-wide history of extensive exploration is broader than any single local family.

# Minima, crossings and classification

A grid minimum is $\widehat x_g^*=\arg\min_{x\in\mathcal G}\widehat f_g(x)$. Multiple shallow minima and boundary solutions can make it unstable. A fixed-risk qualifying set is $\widehat S_g(c)=\{x:RR_g(x)\ge c\}$. Its first upper crossing is not sufficient to describe a disconnected set. A missing crossing means the set does not cross in the assessed domain, not a temperature of zero.

For assessed weather dates A, $I_{gt}=1\{x_t\in\widehat S_g\}$ and $N_g=\sum_{t\in A}I_{gt}$. Classified-date statistics are discontinuous functions of coefficients. Conditional Gaussian coefficient simulation evaluates that discontinuity under a fitted covariance, but does not automatically calibrate confidence limits. A support set based on a band is a statement about statistical support under a procedure, not a physiological damage threshold.

The clean segmented model is $\eta_t=\text{nuisance}+\beta_1x_t+\beta_2(x_t-k)_+$. Its below-threshold slope is $\beta_1$ and above-threshold slope $\beta_1+\beta_2$. Selecting k minimizes in-sample deviance over a finite grid. Standard errors conditional on that selected k ignore search and are not reported as breakpoint confidence intervals.

# Population-day algebra

With fixed scenario population weights $P_{gt}>0$,

$$H_g=\sum_{t\in A}P_{gt}I_{gt},\quad D_g=\sum_{t\in A}P_{gt},\quad F_g=H_g/D_g.$$

Define $\Delta_g=F_{g,own}-F_{g,pooled}$ and the change in a pair gap as $\Delta_g-\Delta_h$. For each group, net population-days equal own-only minus pooled-only population-days. Equal date weights make the shared pooled rule cancel from the pair gap change. Unequal population trajectories can prevent exact cancellation. Population-day uncertainty from response coefficients does not include unknown catchments, race-definition incompatibility or uncertainty in the population estimates themselves.

Date-membership draws are strongly dependent because all dates share coefficients. The variance of a total therefore cannot be calculated by pretending the date flags are independent Bernoulli observations. A zero percentile endpoint can arise from a block of exact-zero draws. Its sign frequency is neither a p-value nor a posterior health-risk probability.

# Penalized GAMs and interactions

A same-day GAM uses smooth functions $s_T(T_t)+s_D(D_t)+s_H(H_t)+s_{time}(t)$ and weekday. The interaction version adds a tensor term $s_{TD}(T_t,D_t)$. Penalized fitting balances fit and roughness, schematically minimizing

$$-2\ell(\beta)+\sum_j\lambda_j\beta^\top S_j\beta.$$

REML chooses smoothing penalties in the installed mgcv implementation. Basis dimension k limits available flexibility; the fitted effective degrees of freedom depend on the penalties. The tensor interaction allows departure from additive weather effects. It is not an independent DTR effect with Tmin/Tmax held fixed, because the variables share underlying extrema. Concurvity diagnoses overlap in what smooth terms can explain. The [mgcv documentation](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/gam.html) describes fitting and smoothness selection.

# Annual county ecological regression

The model is $\log r_{cy}=\alpha+\tau_y+\beta_P(Poverty_{cy}/5)+\beta_B(BlackShare_{cy}/10)+\cdots+\epsilon_{cy}$. The factor $\exp\beta_P$ describes the conditional multiplicative annual rate index associated with +5 poverty percentage points. County-clustered covariance sums score contributions within county, uses the CR1 finite-sample factor, and compares coefficients to t with G−1 degrees of freedom. It addresses repeated county observations, not every source of spatial dependence. No county fixed effects or individual interpretation is implied by this formula.

# Composition and partial identification

Marginal group rates are mixtures of latent cell rates. For age a and race r, a marginal expected count is $\mu_{rt}=\sum_aP_{art}\lambda_{art}$. Knowing the P cells does not identify each $\lambda_{art}$ from separate age and race death margins. Multiple latent tables can generate the same observed margins.

S1 restricts those tables with a common within-age racial multiplier; its point estimate is conditional on that restriction. S3 optimizes a chosen standardized contrast over feasible nonnegative latent allocations satisfying declared fitted margins and positive population weights. Such a feasible-allocation region is an identification statement conditional on the model inputs. Adding a coefficient-draw envelope introduces a separate estimation layer; it is not automatically a confidence interval or a causal bound.

# Validation vocabulary

Numerical reproduction means an independent calculation matches a saved value. Model-conditional uncertainty assumes the fitted model and covariance procedure. Known-truth calibration evaluates coverage in specified synthetic designs. Internal predictive validation checks held-out observed outcomes. Identification asks whether the target is uniquely determined by the observational design and assumptions. These are different achievements: none should be used as a synonym for the others.
'''
    write_report("MATHEMATICAL_APPENDIX.md", "Mathematical appendix: intuition, estimands and assumptions", body)


def build_history(evidence):
    parent = REPORTS
    packages = [
      ("Earlier corrected Tmin consolidation", "Project_main_results/REPORT.md", "Historical anchors show delayed cold and a calendar-sensitive age contrast. Frozen earlier corrections/families are preserved, not pooled with the new percentile models."),
      ("Reanalysis core", "human_biometeorology_reanalysis_v1/reports/FULL_REPORT.md", "Source-cell reconstruction, 64 historical fits, support failures, cause/pollution branches, broad inequality and threshold atlas. Clean primary checkpoints are independently reproduced."),
      ("Response-metric concordance", "3. Results/Course metric concordance/REPORT.md", "One original adjusted-supported observed-calendar age/cold cell; no broad conclusion of equivalence or metric disagreement."),
      ("Stage 3 composition", "3. Results/Stage 3 composition inference/REPORT.md", "Assumption-bound S1, conditional S3 feasible-allocation bounds, corrected Census years, required simulation and optimizer checks; not identified racial standardization."),
      ("Person-day exposure accounting", "3. Results/Person-day exposure assessment/README.md", "Fixed/percentile/MMT weather allocation and population partitions; population-days are distinct from deaths attributable to temperature."),
      ("Population classification", "population_risk_comparison_v1/reports/FULL_REPORT.md", "Pooled versus own-group classifications, paired draw ranges and 56-comparison registry. Earlier population-dependent outputs retain their historical inputs."),
      ("Exact-zero endpoints", "population_risk_comparison_v1/zero_bounds_audit/reports/ZERO_BOUNDS_REPORT.md", "Large exact-zero blocks explain highlighted percentile endpoints. Original per-draw flags were unavailable and reconstruction was not exact."),
      ("Classification reliability", "population_risk_comparison_v1/method_reliability/reports/REPORT.md", "A distinct saved membership-draw series and 800 synthetic refits. Equal-scenario central interval coverage is 88.5%; totals and date identities are uncertain."),
      ("DTR inequality extension", "dtr_inequality_exploration_v1/reports/REPORT.md", "Primary Holm24 points/Holm4 curves uncertain; seven-region native-calendar sex sensitivity and broad projection bounds. Later matched-calendar and calibration audits qualify inference."),
      ("DTR article benchmark", "article_replication_benchmark_v1/reports/CLEAR_RESULTS.md", "Fixed weather events, unique lag unions, matched approximation models and 21 conditional composition ORs. Catchments/source snapshot prevent exact paper reproduction."),
      ("Consolidated follow-up", "biometeorology_consolidated_followup_v1/reports/consolidated_followup.md", "Matched-calendar attenuation, competing female trough reversal, fold-score denominator reconciliation and confidence-region undercoverage."),
      ("Regional architectures", "calendar_geography_extension_v1/reports/Technical_Appendix.md", "M0 aggregate/M1 shared regional/M2 regional responses; 1,899 common fit dates; no new formal inference or regional population allocation."),
      ("Threshold search validation", "3. Results/Segmented threshold validation/README.md", "Calibrated search and held-year evidence are retained. Failure to detect a unique threshold does not imply identical smooth curves."),
      ("Morbidity visits", "3. Results/Morbidity exploration/README.md", "A distinct outcome/reporting process, weather calendar and support limitations; not additional independent evidence for mortality susceptibility."),
      ("County annual mortality", "3. Results/County baseline mortality/README.md", "Annual ecological gradients with clustered/spatial sensitivity; the age-adjusted poverty checkpoint is independently rerun in the clean project."),
    ]
    rows = []
    body = "# How to read inherited evidence\n\nReports are retained because they document completed scientific work. Their conclusions are not frozen against later corrections. The table distinguishes newly reproduced calculations, preserved historical fits and targets that remain unidentified. Old ZIP deliveries and duplicate review packages are historical snapshots, not independent samples. No historical file is edited.\n\n"
    for title, source, interpretation in packages:
        body += "## " + title + "\n\n" + interpretation + "\n\n" + historical(source, "Read the source report", parent) + "\n\n"
        rows.append(dict(Topic=title, Source=source, Clean_scope=interpretation))
    path = ROOT / "1. Data/metadata/reproduction_status.csv"
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=["Topic", "Source", "Clean_scope"])
        writer.writeheader()
        writer.writerows(rows)
    body += "# Important numerical distinctions\n\nThe earlier corrected Tmin result around RR1.088 and the new Tmax P10 result around RR1.083 use different exposures and anchors. Historical poverty ratios around 1.095 and the clean poverty-only value around 1.101 belong to different covariate specifications. Historic observed/short-gap calendars (3,657/5,244), reanalysis matched temperatures (3,604), primary adjusted DTR (2,155), seven-region native DTR (3,330), matched-source DTR (2,093) and common regional architectures (1,899) must not be interchanged. Every substantive comparison should retain metric, outcome, footprint, source, controls, reference, lag, calendar and uncertainty type.\n\n"
    body += "# Scope of rewriting\n\nThe clean runner replaces the central data, descriptive, distributed-lag, direct-comparison, threshold/classification, population-accounting, county and GAM stages with readable scripts. It does not execute the hundreds of exploratory/review scripts as though every historical branch were part of one current model. The specialized Stage 3, regional-architecture, morbidity and large simulation branches remain source-linked historical results. Their outputs are inspectable offline and their status is explicit; this consolidation does not claim independent regeneration of their original latent draws or nonlinear optimizers.\n\n"
    body += "# Complete source register\n\n" + link("1. Data/metadata/evidence_manifest.csv", "Full evidence manifest with source paths and SHA-256", parent) + "\n\n" + link("1. Data/metadata/historical_documents.csv", "Historical document register", parent) + "\n\n" + link("1. Data/metadata/reproduction_status.csv", "Major-result concordance", parent)
    write_report("HISTORICAL_EVIDENCE_REVIEW.md", "Historical evidence review and reproduction concordance", body)


def build_browser(evidence):
    def page(path, title, content):
        path.write_text('<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>' + html.escape(title) + '</title><style>body{max-width:1250px;margin:35px auto;padding:0 20px;font:16px/1.5 system-ui;color:#203343}a{color:#14637e}img{max-width:100%;height:auto}section{border:1px solid #d5e0e8;border-radius:8px;padding:18px;margin:20px 0}input{padding:12px;width:90%;font:inherit}.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(380px,1fr));gap:20px}summary{cursor:pointer}small{color:#526676}</style><h1>' + html.escape(title) + '</h1>' + content + '</html>', encoding="utf-8")
    root = ROOT / "3. Outputs"
    main = '<p>Newly generated plots and readable table images appear first. Historical snapshots retain their original package names and scientific limitations.</p><p><a href="' + url("4. Reports/GENERAL_REPORT.html", root) + '">Read the general report</a></p>'
    for theme, title in zip(THEMES, TITLES):
        folder = root / theme
        body = '<p><a href="../index.html">All themes</a> · <a href="' + url("4. Reports/" + theme + "/REPORT.html", folder) + '">Thematic report</a></p>'
        body += '<h2>New figures</h2><div class="grid">'
        for image in sorted((folder / "png").glob("*.png")):
            relative = quote(str(image.relative_to(folder)))
            body += '<section><h3>' + html.escape(image.stem.replace("_", " ")) + '</h3><a href="' + relative + '"><img loading="lazy" src="' + relative + '"></a></section>'
        body += '</div><details><summary>Paginated table images</summary><div class="grid">'
        for image in sorted((folder / "tables_png").glob("*.png")):
            relative = quote(str(image.relative_to(folder)))
            body += '<section><small>' + html.escape(image.name) + '</small><a href="' + relative + '"><img loading="lazy" src="' + relative + '"></a></section>'
        body += '</div></details><h2>CSV numerical exports</h2><ul>'
        for table in sorted((folder / "csv").glob("*.csv")):
            body += '<li><a href="' + quote(str(table.relative_to(folder))) + '">' + html.escape(table.name) + '</a></li>'
        body += '</ul><h2>Historical evidence</h2><p>Search by package, figure or report name. These files are snapshots, not newly refitted results.</p><input id="search" placeholder="Filter historical files…"><ul id="evidence">'
        for row in evidence:
            if row["Theme"] != theme:
                continue
            body += '<li><a href="' + url(row["Clean_path"], folder) + '">' + html.escape(row["Source"]) + '</a></li>'
        body += '</ul><script>document.getElementById("search").addEventListener("input",function(){const q=this.value.toLowerCase();document.querySelectorAll("#evidence li").forEach(x=>x.style.display=x.textContent.toLowerCase().includes(q)?"":"none")});</script>'
        page(folder / "index.html", title, body)
        main += '<section><h2><a href="' + theme + '/index.html">' + title + '</a></h2><p>' + str(len(list((folder / "png").glob("*.png")))) + ' new figures · ' + str(len(list((folder / "tables_png").glob("*.png")))) + ' table images · ' + str(sum(row["Theme"] == theme for row in evidence)) + ' historical files</p></section>'
    page(root / "index.html", "Virginia biometeorology — visual results browser", main)
    reports = '<p><a href="' + url("3. Outputs/index.html", REPORTS) + '">Open visual results browser</a></p><ul>'
    for document in sorted(REPORTS.rglob("*.html")):
        if document.name == "index.html": continue
        reports += '<li><a href="' + quote(str(document.relative_to(REPORTS))) + '">' + html.escape(str(document.relative_to(REPORTS))) + '</a> · <a href="' + quote(str(document.with_suffix(".pdf").relative_to(REPORTS))) + '">PDF</a></li>'
    page(REPORTS / "index.html", "Scientific reports", reports + '</ul>')


if __name__ == "__main__":
    main()
