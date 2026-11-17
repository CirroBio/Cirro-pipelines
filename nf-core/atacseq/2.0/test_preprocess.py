"""Offline tests for the nf-core/atacseq 2.0 samplesheet handling.

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
    spec = importlib.util.spec_from_file_location("atacseq_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()

HEADER = "sample,fastq_1,fastq_2,replicate\r\n"


def replicate_rows(*rows):
    return pd.DataFrame(rows, columns=["sample", "replicate"])


class ReadSamplesheet(unittest.TestCase):

    def test_values_stay_text(self):
        samplesheet = preprocess.read_samplesheet(io.StringIO(
            HEADER + "01,s3://bucket/data/01_R1.fastq.gz,,1\r\n"
        ))
        self.assertEqual(samplesheet["sample"].tolist(), ["01"])
        self.assertEqual(samplesheet["fastq_2"].tolist(), [""])
        self.assertEqual(samplesheet["replicate"].tolist(), ["1"])

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


class CheckReplicates(unittest.TestCase):

    def test_numbered_from_one_passes(self):
        # Two lanes of replicate 1 are merged by the pipeline, so a repeated number is fine
        preprocess.check_replicates(replicate_rows(
            ("A", "1"), ("A", "1"), ("A", "2"), ("B", "1"),
        ))

    def test_gap_fails(self):
        with self.assertRaisesRegex(ValueError, "sample A must be numbered 1 to 2, not 1, 3"):
            preprocess.check_replicates(replicate_rows(("A", "1"), ("A", "3"), ("B", "1")))

    def test_not_starting_at_one_fails(self):
        with self.assertRaisesRegex(ValueError, "sample B must be numbered 1 to 1, not 2"):
            preprocess.check_replicates(replicate_rows(("A", "1"), ("B", "2")))

    def test_blank_or_non_numeric_fails(self):
        for value in ("", "one", "1.5"):
            with self.assertRaisesRegex(ValueError, "whole-number replicate; check sample\\(s\\) B"):
                preprocess.check_replicates(replicate_rows(("A", "1"), ("B", value)))


if __name__ == "__main__":
    unittest.main()
