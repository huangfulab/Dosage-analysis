#!/usr/bin/env python3
"""Validate every supported R baseline and emit all or one selected series."""
import argparse
import csv
import datetime
import json
from pathlib import Path
import re
import urllib.request

ROOT = Path(__file__).resolve().parent


def matrix(release_version, r_series="all"):
    with (ROOT / "r-versions.csv").open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    with (ROOT / "minimum-packages.csv").open(newline="") as stream:
        reader = csv.DictReader(stream)
        fields = reader.fieldnames
        packages = list(reader)
    series = [row["r"] for row in rows]
    if not series or len(set(series)) != len(series):
        raise ValueError("R matrix must be nonempty and unique")
    if series != [f"4.{minor}" for minor in range(len(series))]:
        raise ValueError("R 4.0 onward must be present without gaps")
    if fields != ["package", "source", *series]:
        raise ValueError("Dependency columns must match every R matrix entry")
    required = {"drc", "ggplot2", "ggrepel", "openxlsx", "pheatmap", "devtools", "matrixStats", "DESeq2", "GSVA"}
    if {row["package"] for row in packages} != required or len(packages) != len(required):
        raise ValueError("Dependency matrix must contain all nine direct dependencies exactly once")
    for row in rows:
        datetime.date.fromisoformat(row["cran_snapshot"])
        if not re.fullmatch(r"3\.\d+", row["bioconductor"]):
            raise ValueError("Invalid Bioconductor release")
    for row in packages:
        if row["source"] not in {"CRAN", "Bioconductor"}:
            raise ValueError("Unknown dependency source")
        if any(not re.fullmatch(r"\d+(?:[.-]\d+)+", row[r]) for r in series):
            raise ValueError("Every package must have an exact version for every R series")
    match = re.fullmatch(r"(\d+)\.(\d+)\.\d+", release_version)
    if not match:
        raise ValueError("Release service returned an invalid stable R version")
    current = f"{match[1]}.{match[2]}"
    if current not in series:
        raise ValueError(f"R {current} has no dependency baseline; extend BOTH CI CSV files before claiming support")
    if r_series != "all":
        if r_series not in series:
            raise ValueError(f"Unknown R series selection: {r_series!r}; use all or a declared minor series")
        rows = [row for row in rows if row["r"] == r_series]
    return {"include": rows}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release-version", help="Offline validation override; CI queries the R-hub stable-release service")
    parser.add_argument("--r-series", default="all", help="Emit all declared series (default) or one minor series, e.g. 4.6")
    args = parser.parse_args()
    release = args.release_version
    if release is None:
        with urllib.request.urlopen("https://api.r-hub.io/rversions/r-release", timeout=30) as response:
            release = json.load(response)["version"]
    print(json.dumps(matrix(release, r_series=args.r_series), separators=(",", ":")))


if __name__ == "__main__":
    main()
