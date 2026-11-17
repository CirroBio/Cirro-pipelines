"""Offline tests for the nf-core/mag 3.0 samplesheet handling.

preprocess.py imports the Cirro SDK only when it runs, so this is a plain
`python -m unittest` that needs pandas and nothing else.

The samplesheet is the structured input Cirro writes from the run form: a CSV with
CRLF line endings, a header from the form's columns, and blank cells left empty.
"""
import importlib.util
import io
import unittest
from pathlib import Path

import pandas as pd


def _load_preprocess():
    path = Path(__file__).with_name("preprocess.py")
    spec = importlib.util.spec_from_file_location("mag_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()

HEADER = "sample,run,group,short_reads_1,short_reads_2,long_reads\r\n"


class ReadSamplesheet(unittest.TestCase):

    def test_values_stay_text(self):
        samplesheet = preprocess.read_samplesheet(io.StringIO(
            HEADER + "01,,,s3://bucket/data/01_R1.fastq.gz,s3://bucket/data/01_R2.fastq.gz,\r\n"
        ))
        self.assertEqual(samplesheet["sample"].tolist(), ["01"])
        self.assertEqual(samplesheet["run"].tolist(), [""])
        self.assertEqual(samplesheet["group"].tolist(), [""])

    def test_no_param_fails_clearly(self):
        for path in (None, ""):
            with self.assertRaisesRegex(ValueError, "No samplesheet was provided"):
                preprocess.read_samplesheet(path)

    def test_empty_file_fails_clearly(self):
        with self.assertRaisesRegex(ValueError, "The samplesheet is empty"):
            preprocess.read_samplesheet(io.StringIO(""))

    def test_header_without_rows_fails_clearly(self):
        with self.assertRaisesRegex(ValueError, "The samplesheet is empty"):
            preprocess.read_samplesheet(io.StringIO(HEADER))


class FillGroup(unittest.TestCase):

    def test_blank_becomes_zero_and_set_values_are_kept(self):
        # A partly annotated table: only some samples were given a group
        samplesheet = pd.DataFrame({"sample": ["A", "B", "C"], "group": ["", "gut", ""]})
        filled = preprocess.fill_group(samplesheet)
        self.assertEqual(filled["group"].tolist(), ["0", "gut", "0"])
        self.assertEqual(filled.columns.tolist(), ["sample", "group"])


class CheckUniqueRuns(unittest.TestCase):

    def test_distinct_runs_and_single_rows_pass(self):
        preprocess.check_unique_runs(pd.DataFrame({
            "sample": ["A", "A", "B", "C"],
            "run": ["1", "2", "", "1"],
        }))

    def test_repeated_run_fails(self):
        with self.assertRaisesRegex(ValueError, "check sample\\(s\\) A$"):
            preprocess.check_unique_runs(pd.DataFrame({
                "sample": ["A", "A", "B"],
                "run": ["1", "1", "1"],
            }))

    def test_two_blank_runs_fail(self):
        # One row per lane with no run ID
        with self.assertRaisesRegex(ValueError, "check sample\\(s\\) A$"):
            preprocess.check_unique_runs(pd.DataFrame({"sample": ["A", "A"], "run": ["", ""]}))

    def test_blank_and_zero_collide(self):
        with self.assertRaisesRegex(ValueError, "check sample\\(s\\) A$"):
            preprocess.check_unique_runs(pd.DataFrame({"sample": ["A", "A"], "run": ["", "0"]}))


class IsSingleEnd(unittest.TestCase):

    def test_all_paired(self):
        self.assertFalse(preprocess.is_single_end(pd.DataFrame({"short_reads_2": ["a_R2.fq.gz", "b_R2.fq.gz"]})))

    def test_all_single(self):
        self.assertTrue(preprocess.is_single_end(pd.DataFrame({"short_reads_2": ["", ""]})))

    def test_mixed_fails(self):
        with self.assertRaisesRegex(ValueError, "all single-end or all paired-end"):
            preprocess.is_single_end(pd.DataFrame({"short_reads_2": ["a_R2.fq.gz", ""]}))


if __name__ == "__main__":
    unittest.main()
