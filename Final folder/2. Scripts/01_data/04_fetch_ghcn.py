"""04_fetch_ghcn.py - download public NOAA GHCN-Daily TMIN/TMAX, 2005-2020 (donor data for weather gaps).

For each of the 12 regions: the GHCN candidate listed in
1. Data/reference/metadata/stations.csv, plus up to 3 other GHCN
stations within ~0.35 deg that report TMIN over the period (backup donors for
gap filling, as Davis et al. 2023 did with GHCN/PRISM).
Output: 1. Data/raw/ghcn/<STATION>.csv and ghcn_stations.csv.
Only public station IDs and dates are sent to NOAA.  Standard library only.
Already-downloaded station files are kept (the step is skipped when all files exist).
Run from the "Final folder" root:  python3 -I "2. Scripts/01_data/04_fetch_ghcn.py"
"""
import csv, json, time, urllib.request, urllib.parse
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
META = ROOT / "1. Data/reference/metadata/stations.csv"
OUT = ROOT / "1. Data/raw/ghcn"
OUT.mkdir(parents=True, exist_ok=True)
API = "https://www.ncei.noaa.gov/access/services"

def get(url, tries=4):
    for k in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return r.read().decode("utf-8")
        except Exception as e:
            time.sleep(5 * (k + 1)); err = e
    raise err

if (OUT / "ghcn_stations.csv").exists():
    print("GHCN files already present; nothing to download."); raise SystemExit(0)
stations = list(csv.DictReader(open(META, encoding="utf-8")))
rows = []
for s in stations:
    lat, lon, cand = float(s["Latitude"]), float(s["Longitude"]), s["GHCN_candidate"]
    q = urllib.parse.urlencode({"dataset": "daily-summaries", "bbox": f"{lat+.35},{lon-.35},{lat-.35},{lon+.35}",
                                "startDate": "2005-01-01", "endDate": "2020-12-31", "dataTypes": "TMIN", "limit": 50})
    res = json.loads(get(f"{API}/search/v1/data?{q}")).get("results", [])
    near = []
    for r in res:
        sid = r["stations"][0]["id"]
        if sid == cand or not sid.startswith("US"): continue
        c = r.get("location", {}).get("coordinates", [None, None])
        near.append((sid, r["stations"][0]["name"], c[1], c[0], r.get("startDate", ""), r.get("endDate", "")))
    near = [n for n in near if n[4][:4] <= "2008" and n[5][:4] >= "2018"][:3]      # long records only
    ids = [(cand, "candidate", s["GHCN_name"], lat, lon)] + [(n[0], "donor", n[1], n[2], n[3]) for n in near]
    for sid, role, name, la, lo in ids:
        f = OUT / f"{sid}.csv"
        if not f.exists():
            q = urllib.parse.urlencode({"dataset": "daily-summaries", "stations": sid, "startDate": "2005-01-01",
                                        "endDate": "2020-12-31", "dataTypes": "TMIN,TMAX", "format": "csv", "units": "metric"})
            f.write_text(get(f"{API}/data/v1?{q}"), encoding="utf-8")
        rows.append({"Region": s["Region"], "Station": sid, "Role": role, "Name": name, "Latitude": la, "Longitude": lo})
    print(s["Region"], [i[0] for i in ids], flush=True)
with open(OUT / "ghcn_stations.csv", "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
