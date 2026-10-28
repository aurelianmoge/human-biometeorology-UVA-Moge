"""Copy immutable data and catalogue the exploration evidence without modifying it.

This is the only stage that needs the Exploration directory. Subsequent stages
read exclusively from the clean project's own data and evidence snapshots.
"""
from pathlib import Path
import csv
import hashlib
import shutil

ROOT = Path(__file__).resolve().parents[2]
EXPLORATION = ROOT.parent / "Exploration"
THEMES = ["01_data_description", "02_descriptive_statistics", "03_modelization", "04_original_research"]


def write_csv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def digest(path):
    hasher = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            hasher.update(block)
    return hasher.hexdigest()


def theme_for(path):
    # Route whole packages together so original relative report/image links survive.
    parts = path.parts
    container = Path(*parts[:2]) if parts[0] == "3. Results" and len(parts) > 1 else Path(parts[0])
    label = str(container).lower()
    if any(word in label for word in ["data audit", "crosswalk", "geography", "weather validation", "specification and geography"]):
        return THEMES[0]
    if any(word in label for word in ["descriptive", "composition context", "count rate diagnosis"]):
        return THEMES[1]
    if any(word in label for word in ["replication", "threshold", "dlnm", "population mmt"]):
        return THEMES[2]
    return THEMES[3]


def main():
    # A partial 845 MB download is not a usable input. Keep it in Exploration.
    raw_manifest = []
    for source in sorted((EXPLORATION / "1. Data").rglob("*")):
        if not source.is_file() or source.name.startswith(".") or ".download" in str(source):
            continue
        relative = source.relative_to(EXPLORATION / "1. Data")
        destination = ROOT / "1. Data/raw" / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        if not destination.exists() or digest(destination) != digest(source):
            shutil.copy2(source, destination)
        raw_manifest.append(dict(Source=str(source.relative_to(ROOT.parent)), Clean_path=str(destination.relative_to(ROOT)), Bytes=source.stat().st_size, SHA256=digest(source)))
    write_csv(ROOT / "1. Data/metadata/raw_manifest.csv", raw_manifest, ["Source", "Clean_path", "Bytes", "SHA256"])

    evidence = []
    report_rows = []
    # Retain reports with their original relative image/table structure: links
    # therefore continue to resolve inside each snapshot when assets are present.
    allowed = {".csv", ".png", ".pdf", ".md", ".html", ".txt", ".qmd", ".docx", ".json", ".geojson", ".css", ".gz"}
    for source in sorted(EXPLORATION.rglob("*")):
        if not source.is_file() or source.suffix.lower() not in allowed:
            continue
        relative = source.relative_to(EXPLORATION)
        if relative.parts[0] == "1. Data" or "2. Scripts" in relative.parts or "scripts" in relative.parts or source.name.startswith("."):
            continue
        theme = theme_for(relative)
        # Snapshots are historical evidence, never advertised as newly refitted.
        destination = ROOT / "3. Outputs" / theme / "historical_evidence" / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        if not destination.exists() or source.stat().st_size != destination.stat().st_size:
            shutil.copy2(source, destination)
        record = dict(Theme=theme, Source=str(relative), Clean_path=str(destination.relative_to(ROOT)), Extension=source.suffix.lower(), Bytes=source.stat().st_size, SHA256=digest(source), Status="historical snapshot; not independently refitted by collection")
        evidence.append(record)
        if source.suffix.lower() in {".md", ".pdf", ".docx", ".qmd"}:
            report_rows.append(record)
    fields = ["Theme", "Source", "Clean_path", "Extension", "Bytes", "SHA256", "Status"]
    write_csv(ROOT / "1. Data/metadata/evidence_manifest.csv", evidence, fields)
    write_csv(ROOT / "1. Data/metadata/historical_documents.csv", report_rows, fields)

    # Freeze the numerical/reference inputs used for independent reproduction.
    reference = ROOT / "1. Data/reference"
    reference.mkdir(parents=True, exist_ok=True)
    package = EXPLORATION / "human_biometeorology_reanalysis_v1"
    for folder in ["metadata", "config", "results/csv"]:
        for source in (package / folder).glob("*.csv"):
            destination = reference / folder / source.name
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, destination)
    for relative in ["3. Results/County baseline mortality/county_baseline_analysis_data.csv"]:
        shutil.copy2(EXPLORATION / relative, reference / Path(relative).name)
    print(f"Collected {len(raw_manifest)} data sources, {len(evidence)} evidence files, {len(report_rows)} historical documents.")


if __name__ == "__main__":
    main()
