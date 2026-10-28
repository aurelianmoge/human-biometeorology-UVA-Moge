"""Extract original Excel cells and apply only corroborated date corrections.

The XML reader preserves Excel row numbers and error tokens. It does not use
an R session's automatic date rounding. Celsius weather values are never
reconstructed from spot temperatures or a conflicting Fahrenheit observation.
"""
from pathlib import Path
from collections import Counter
import csv
import datetime as dt
import math
import re
import xml.etree.ElementTree as ET
import zipfile

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "1. Data/derived"
NS = {"s": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
ORIGIN = dt.date(1899, 12, 30)


def number(value):
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (ValueError, TypeError):
        return None


def write(name, rows, fields):
    with (OUTPUT / (name + ".csv")).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def read_rows(path):
    with zipfile.ZipFile(path) as workbook:
        shared = []
        if "xl/sharedStrings.xml" in workbook.namelist():
            tree = ET.fromstring(workbook.read("xl/sharedStrings.xml"))
            shared = ["".join(node.itertext()) for node in tree.findall("s:si", NS)]
        rows = ET.fromstring(workbook.read("xl/worksheets/sheet1.xml")).findall("s:sheetData/s:row", NS)
        def value(cell):
            if cell.attrib.get("t") == "inlineStr":
                return "".join(cell.find("s:is", NS).itertext())
            node = cell.find("s:v", NS)
            token = node.text if node is not None else ""
            return shared[int(token)] if cell.attrib.get("t") == "s" else token
        headers = {re.sub(r"\d+$", "", cell.attrib["r"]): value(cell) for cell in rows[0]}
        # Duplicate source labels get explicit positional suffixes in dictionary.
        for row in rows[1:]:
            yield int(row.attrib["r"]), {headers.get(re.sub(r"\d+$", "", cell.attrib["r"]), cell.attrib["r"]): value(cell) for cell in row}


def main():
    records = []
    for path in sorted((ROOT / "1. Data/raw").glob("*.xlsx")):
        region = path.name.split(".Final")[0]
        for row_number, source in read_rows(path):
            serial = number(source.get("Date"))
            if serial is None:
                continue
            date = ORIGIN + dt.timedelta(days=math.floor(serial))
            if not dt.date(2005, 1, 1) <= date <= dt.date(2020, 12, 31):
                raise ValueError(f"Out-of-period date in {region}, row {row_number}")
            try:
                ymd = dt.date(*[int(float(source[field])) for field in ["Year", "Month", "Day"]])
            except (ValueError, KeyError):
                ymd = None
            record = dict(Region=region, Excel_row=row_number, Date=str(date), Original_date=str(date), Serial=serial, YMD=str(ymd) if ymd else "", Date_repaired=False)
            fields = {"All": "Mort", "Tmin": "MinT(C)", "Tmax": "MaxT(C)", "DTR_supplied": "DTR(C)", "T07": "T(C)7am", "Td07": "Td(C)7am", "RH07": "RH7am", "Tafternoon": "T(C)1pm", "Wind": "Speedkts1pm", "Pressure": "SLPhPa1pm", "PM": "PM2.5RAW", "Ozone": "OzoneRAW", "Respiratory": "MDC.04", "Circulatory": "MDC.05"}
            fields.update({group: group for group in ["Black", "White", "Female", "Male"]})
            for target, field in fields.items():
                record[target] = number(source.get(field))
            for target, fields in {"Under50": ["X0.11", "X12.17", "X18.29", "X30.39", "X40.49"], "Age50_74": ["X50.64", "X65.74"], "Age75plus": ["X75."]}.items():
                values = [number(source.get(field)) for field in fields]
                record[target] = sum(values) if all(value is not None for value in values) else None
            records.append(record)
    counts = Counter((row["Region"], row["Date"]) for row in records)
    lookup = {(row["Region"], row["Excel_row"]): row for row in records}
    repairs = []
    for row in records:
        if not row["YMD"] or row["YMD"] == row["Date"]:
            continue
        target = dt.date.fromisoformat(row["YMD"])
        previous = lookup.get((row["Region"], row["Excel_row"] - 1))
        following = lookup.get((row["Region"], row["Excel_row"] + 1))
        corroborated = (
            target == dt.date.fromisoformat(row["Date"]) + dt.timedelta(days=1)
            and row["Serial"] % 1 > .99
            and target == dt.date(2005, 1, 1) + dt.timedelta(days=row["Excel_row"] - 2)
            and counts[(row["Region"], row["YMD"])] == 0
            and counts[(row["Region"], row["Date"])] > 1
            and previous is not None and following is not None
            and previous["YMD"] == str(target - dt.timedelta(days=1))
            and following["YMD"] == str(target + dt.timedelta(days=1))
        )
        if corroborated:
            row["Date"] = row["YMD"]
            row["Date_repaired"] = True
            repairs.append({field: row[field] for field in ["Region", "Excel_row", "Original_date", "Date"]})
    OUTPUT.mkdir(parents=True, exist_ok=True)
    write("workbook_cells", records, list(records[0]))
    write("date_repairs", repairs, ["Region", "Excel_row", "Original_date", "Date"])
    print(f"Extracted {len(records)} source rows; {len(repairs)} internally corroborated date repairs.")


if __name__ == "__main__":
    main()
