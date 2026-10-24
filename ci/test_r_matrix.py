import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

import r_matrix


class MatrixTests(unittest.TestCase):
    def test_all_r_series_and_bioc_pairings(self):
        rows = r_matrix.matrix("4.6.1")["include"]
        self.assertEqual([row["r"] for row in rows], [f"4.{i}" for i in range(7)])
        self.assertEqual([row["bioconductor"] for row in rows], ["3.11", "3.13", "3.15", "3.17", "3.19", "3.21", "3.23"])

    def test_new_release_must_not_be_silently_omitted(self):
        for release in ("4.7.0", "5.0.0", "4.6.1-devel"):
            with self.assertRaises(ValueError):
                r_matrix.matrix(release)

    def test_cli_json_is_valid(self):
        result = subprocess.run([sys.executable, str(r_matrix.ROOT / "r_matrix.py"),
            "--release-version", "4.6.1"], capture_output=True, text=True, check=True)
        self.assertEqual(len(json.loads(result.stdout)["include"]), 7)

    def test_missing_or_duplicate_baseline_is_rejected(self):
        for change in ("missing", "duplicate", "empty-floor"):
            with tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                for name in ("r-versions.csv", "minimum-packages.csv"):
                    (root / name).write_bytes((r_matrix.ROOT / name).read_bytes())
                if change == "empty-floor":
                    path = root / "minimum-packages.csv"
                    with path.open(newline="") as stream:
                        reader = csv.DictReader(stream)
                        fields, rows = reader.fieldnames, list(reader)
                    rows[0]["4.0"] = ""
                else:
                    path = root / "r-versions.csv"
                    with path.open(newline="") as stream:
                        reader = csv.DictReader(stream)
                        fields, rows = reader.fieldnames, list(reader)
                    rows = rows[1:] if change == "missing" else rows + [rows[0]]
                with path.open("w", newline="") as stream:
                    writer = csv.DictWriter(stream, fieldnames=fields)
                    writer.writeheader()
                    writer.writerows(rows)
                with mock.patch.object(r_matrix, "ROOT", root), self.assertRaises(ValueError):
                    r_matrix.matrix("4.6.1")


if __name__ == "__main__":
    unittest.main()
