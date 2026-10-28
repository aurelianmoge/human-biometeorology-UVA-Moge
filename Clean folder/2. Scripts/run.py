"""Reproduce the clean project from any working directory.

Examples:
    python3 'Clean folder/2. Scripts/run.py'
    python3 'Clean folder/2. Scripts/run.py' --reports-only
    python3 'Clean folder/2. Scripts/run.py' --verify-only

Only --import-exploration refreshes the frozen historical evidence after the
initial collection. Analytical scripts read the clean project's local inputs.
"""
from pathlib import Path
import argparse
import csv
import datetime as dt
import hashlib
import platform
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
LOGS = ROOT / "3. Outputs/01_data_description/logs/pipeline"
STAGES = [
    "01_data_description/02_extract_workbooks.py",
    "01_data_description/03_assemble_panels.R",
    "01_data_description/04_census_populations.py",
    "02_descriptive_statistics/01_describe_data.R",
    "03_modelization/01_fit_models.R",
    "03_modelization/02_minima_and_thresholds.R",
    "03_modelization/03_blocked_validation.R",
    "03_modelization/04_segmented_thresholds.R",
    "04_original_research/01_direct_comparisons.R",
    "04_original_research/02_classification_accounting.R",
    "04_original_research/03_county_inequalities.R",
    "04_original_research/04_gam_interactions.R",
    "04_original_research/05_population_allocations.R",
    "04_original_research/06_research_atlas.R",
    "00_shared/verify_results.R",
    "00_shared/render_tables.R",
    "00_shared/build_reports.py",
]


def run_stage(relative):
    script = ROOT / "2. Scripts" / relative
    log = LOGS / (relative.replace("/", "__") + ".log")
    if script.suffix == ".R":
        # Show every warning in the log immediately, rather than truncating at 50.
        expression = 'options(warn = 1); source(commandArgs(trailingOnly = TRUE)[1])'
        command = ["Rscript", "--vanilla", "-e", expression, str(script)]
    else:
        command = [sys.executable, str(script)]
    print("Running", relative, flush=True)
    start = time.monotonic()
    with log.open("w") as stream:
        completed = subprocess.run(command, cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT)
    if completed.returncode:
        print(log.read_text()[-7000:], file=sys.stderr)
        raise RuntimeError(f"Stage failed: {relative}. See {log}")
    return dict(Stage=relative, Exit_code=completed.returncode, Seconds=round(time.monotonic()-start, 3), Log=str(log.relative_to(ROOT)))


def write_csv(path, rows, fields):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--import-exploration", action="store_true", help="Refresh immutable data/evidence copies deliberately.")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--reports-only", action="store_true")
    mode.add_argument("--verify-only", action="store_true")
    arguments = parser.parse_args()
    LOGS.mkdir(parents=True, exist_ok=True)
    executions = []
    if arguments.import_exploration or not (ROOT / "1. Data/metadata/raw_manifest.csv").exists():
        executions.append(run_stage("01_data_description/01_collect_sources.py"))
    if arguments.reports_only:
        stages = ["00_shared/build_reports.py"]
    elif arguments.verify_only:
        stages = ["00_shared/verify_results.R"]
    else:
        stages = STAGES
    for stage in stages:
        executions.append(run_stage(stage))
    software = dict(Python=sys.version, Platform=platform.platform(), Run_UTC=dt.datetime.now(dt.timezone.utc).isoformat())
    (LOGS / "python_environment.txt").write_text("\n".join(f"{key}: {value}" for key,value in software.items()))
    write_csv(LOGS / "execution.csv", executions, ["Stage", "Exit_code", "Seconds", "Log"])
    scripts = []
    for path in sorted((ROOT / "2. Scripts").rglob("*")):
        if path.suffix in {".R", ".py"}:
            scripts.append(dict(Path=str(path.relative_to(ROOT)), SHA256=hashlib.sha256(path.read_bytes()).hexdigest()))
    write_csv(ROOT / "1. Data/metadata/clean_script_manifest.csv", scripts, ["Path", "SHA256"])
    print("All requested stages completed. Start with 4. Reports/GENERAL_REPORT.html or 3. Outputs/index.html.", flush=True)


if __name__ == "__main__":
    main()
