import csv
import json
from pathlib import Path
import re
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

    def test_job_environment_uses_runtime_runner_paths(self):
        workflow = (r_matrix.ROOT.parent / ".github/workflows/r-compatibility.yaml").read_text()
        # GitHub does not allow the runner expression context in job-level env.
        # Resolve temporary paths in a step through RUNNER_TEMP/GITHUB_ENV.
        env_block = re.search(r"^    env:\n(.*?)^    defaults:", workflow, re.M | re.S)
        self.assertIsNotNone(env_block)
        self.assertNotIn("${{ runner.", env_block[1])
        self.assertIn('echo "CI_R_LIBRARY=$RUNNER_TEMP/dosage-r-library"', workflow)
        self.assertIn('echo "CI_LOG_DIR=$RUNNER_TEMP/dosage-ci-logs"', workflow)
        self.assertLess(workflow.index("Initialize isolated runner paths"), workflow.index("r-lib/actions/setup-r"))

    def test_native_dependencies_are_checked_before_r_installation(self):
        workflow = (r_matrix.ROOT.parent / ".github/workflows/r-compatibility.yaml").read_text()
        step = re.search(
            r"^      - name: Install system libraries for the hard dependency closure\n"
            r"(.*?)(?=^      - )", workflow, re.M | re.S
        )
        self.assertIsNotNone(step)
        apt = re.search(
            r"sudo apt-get install -y --no-install-recommends[ \t]+([^\n]+)",
            step[1].replace("\\\n", " ")
        )
        self.assertIsNotNone(apt)
        # Recent fs/ragg transitive releases need system libuv and WebP headers.
        for package in ("libuv1-dev", "libwebp-dev"):
            self.assertIn(package, apt[1].split())
        check = "pkg-config --modversion libuv freetype2 libpng libtiff-4 libjpeg libwebp libwebpmux"
        self.assertIn(check, step[1])
        self.assertLess(step[1].index("sudo apt-get install"), step[1].index(check))
        self.assertLess(workflow.index(check), workflow.index("Rscript ci/dependencies.R --install"))

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
