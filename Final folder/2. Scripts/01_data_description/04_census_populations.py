"""Rebuild public Census age × race × sex cells with audited July-year mapping.

These cells supply population descriptions and allocation scenarios. They do not
supply the missing crossed mortality counts or verified workbook denominators.
"""
from pathlib import Path
from collections import defaultdict
import csv

ROOT = Path(__file__).resolve().parents[2]


def write(path, rows, fields):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main():
    county = defaultdict(int)
    age_totals = defaultdict(int)
    published_totals = defaultdict(int)
    mappings = []
    sources = [
        ("co-est00int-alldata-51.csv", True, {code: code + 1998 for code in range(7, 12)}),
        ("CC-EST2020-ALLDATA6-51.csv", False, {**{code: code + 2007 for code in range(3, 13)}, 14: 2020}),
    ]
    for filename, old, years in sources:
        for code, year in years.items():
            mappings.append(dict(Source=filename, Input_YEAR=code, Calendar_year=year, Reference_date=f"{year}-07-01"))
        with (ROOT / "1. Data/raw" / filename).open() as stream:
            for row in csv.DictReader(stream):
                if int(row["SUMLEV"]) != 50 or int(row["STATE"]) != 51 or int(row["COUNTY"]) == 0:
                    continue
                code = int(row["YEAR"])
                if code not in years:
                    continue
                year = years[code]
                fips = f"51{int(row['COUNTY']):03d}"
                age = int(row["AGEGRP"])
                if age == (99 if old else 0):
                    published_totals[(year, fips)] += int(row["TOT_POP"])
                    continue
                valid_ages = range(0, 19) if old else range(1, 19)
                if age not in valid_ages:
                    continue
                age_totals[(year, fips)] += int(row["TOT_POP"])
                # The older upload has a separate infant bin (AGEGRP 0).
                band = "Under50" if age <= 10 else "Age50_74" if age <= 15 else "Age75plus"
                harmonized_fips = "51019" if fips == "51515" else fips
                for sex, suffix in [("Female", "FEMALE"), ("Male", "MALE")]:
                    total = int(row["TOT_" + suffix])
                    black = int(row["BA_" + suffix])
                    white = int(row["WA_" + suffix])
                    if total < black + white:
                        raise ValueError("Race-alone counts exceed total population.")
                    for race, count in [("Black", black), ("White", white), ("Residual", total - black - white)]:
                        county[(year, harmonized_fips, band, race, sex)] += count
    if age_totals != published_totals:
        raise ValueError("Age-bin sums do not reproduce source county-year totals.")
    county_rows = [dict(Year=key[0], FIPS=key[1], Age=key[2], Race=key[3], Sex=key[4], Population=count) for key, count in sorted(county.items())]
    statewide = defaultdict(int)
    for (year, fips, age, race, sex), count in county.items():
        for group in ["All", age, race, sex]:
            statewide[(year, group)] += count
    state_rows = [dict(Year=key[0], Group=key[1], Population=count) for key, count in sorted(statewide.items())]
    output = ROOT / "1. Data/derived"
    write(output / "county_joint_populations.csv", county_rows, ["Year", "FIPS", "Age", "Race", "Sex", "Population"])
    write(output / "statewide_population_margins.csv", state_rows, ["Year", "Group", "Population"])
    write(output / "census_year_mapping.csv", mappings, ["Source", "Input_YEAR", "Calendar_year", "Reference_date"])
    print(f"Verified {len(age_totals)} county-year age sums; rebuilt {len(county_rows)} harmonized population cells.")


if __name__ == "__main__":
    main()
