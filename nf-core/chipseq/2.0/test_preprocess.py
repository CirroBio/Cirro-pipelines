"""Offline tests for the nf-core/chipseq 2.0 samplesheet handling.

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
    spec = importlib.util.spec_from_file_location("chipseq_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()

HEADER = "sample,fastq_1,fastq_2,replicate,antibody,control,control_replicate\r\n"


def replicate_rows(*rows):
    return pd.DataFrame(rows, columns=["sample", "replicate"])


def control_rows(*rows):
    return pd.DataFrame(rows, columns=["sample", "replicate", "antibody", "control", "control_replicate"])


INPUT_ROWS = [
    ("WT_INPUT", "1", "", "", ""),
    ("WT_INPUT", "2", "", "", ""),
]


class ReadSamplesheet(unittest.TestCase):

    def test_values_stay_text(self):
        samplesheet = preprocess.read_samplesheet(io.StringIO(
            HEADER + "01,s3://bucket/data/01_R1.fastq.gz,,1,,,\r\n"
        ))
        self.assertEqual(samplesheet["sample"].tolist(), ["01"])
        self.assertEqual(samplesheet["fastq_2"].tolist(), [""])
        self.assertEqual(samplesheet["replicate"].tolist(), ["1"])
        self.assertEqual(samplesheet["control_replicate"].tolist(), [""])

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

    def test_blank_or_non_numeric_fails(self):
        for value in ("", "one", "1.5"):
            with self.assertRaisesRegex(ValueError, "whole-number replicate; check sample\\(s\\) B"):
                preprocess.check_replicates(replicate_rows(("A", "1"), ("B", value)))


class CheckControls(unittest.TestCase):

    def test_chip_rows_naming_input_replicates_pass(self):
        preprocess.check_controls(control_rows(
            *INPUT_ROWS,
            ("WT_IP", "1", "BCATENIN", "WT_INPUT", "1"),
            ("WT_IP", "2", "BCATENIN", "WT_INPUT", "2"),
        ))

    def test_inputs_only_pass(self):
        preprocess.check_controls(control_rows(*INPUT_ROWS))

    def test_antibody_without_control_fails(self):
        with self.assertRaisesRegex(ValueError, "both antibody and control, or neither; check sample\\(s\\) WT_IP"):
            preprocess.check_controls(control_rows(*INPUT_ROWS, ("WT_IP", "1", "BCATENIN", "", "")))

    def test_control_without_antibody_fails(self):
        with self.assertRaisesRegex(ValueError, "both antibody and control, or neither; check sample\\(s\\) WT_IP"):
            preprocess.check_controls(control_rows(*INPUT_ROWS, ("WT_IP", "1", "", "WT_INPUT", "1")))

    def test_blank_or_non_numeric_control_replicate_fails(self):
        for value in ("", "one"):
            with self.assertRaisesRegex(ValueError, "whole-number control replicate; check sample\\(s\\) WT_IP"):
                preprocess.check_controls(control_rows(
                    *INPUT_ROWS, ("WT_IP", "1", "BCATENIN", "WT_INPUT", value),
                ))

    def test_unknown_control_sample_fails(self):
        with self.assertRaisesRegex(ValueError, "WT_IP names replicate 1 of KO_INPUT as its control"):
            preprocess.check_controls(control_rows(*INPUT_ROWS, ("WT_IP", "1", "BCATENIN", "KO_INPUT", "1")))

    def test_unknown_control_replicate_fails(self):
        with self.assertRaisesRegex(ValueError, "WT_IP names replicate 3 of WT_INPUT as its control"):
            preprocess.check_controls(control_rows(*INPUT_ROWS, ("WT_IP", "1", "BCATENIN", "WT_INPUT", "3")))

    def test_own_sample_as_control_fails(self):
        with self.assertRaisesRegex(ValueError, "Sample WT_IP names itself as its control"):
            preprocess.check_controls(control_rows(
                *INPUT_ROWS,
                ("WT_IP", "1", "BCATENIN", "WT_INPUT", "1"),
                ("WT_IP", "2", "BCATENIN", "WT_IP", "1"),
            ))


if __name__ == "__main__":
    unittest.main()
