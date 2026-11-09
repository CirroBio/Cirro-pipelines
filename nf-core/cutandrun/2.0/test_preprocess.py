"""Offline tests for the nf-core/cutandrun samplesheet handling.

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
    spec = importlib.util.spec_from_file_location("cutandrun_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()

HEADER = "group,replicate,fastq_1,fastq_2,control\r\n"


def replicate_rows(*rows):
    return pd.DataFrame(rows, columns=["group", "replicate"])


def control_rows(*rows):
    return pd.DataFrame(rows, columns=["group", "control"])


class ReadSamplesheet(unittest.TestCase):

    def test_values_stay_text(self):
        samplesheet = preprocess.read_samplesheet(io.StringIO(
            HEADER + "01,1,s3://bucket/data/01_R1.fastq.gz,s3://bucket/data/01_R2.fastq.gz,\r\n"
        ))
        self.assertEqual(samplesheet["group"].tolist(), ["01"])
        self.assertEqual(samplesheet["replicate"].tolist(), ["1"])
        self.assertEqual(samplesheet["control"].tolist(), [""])

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
        # Technical replicates share a number and are merged by the pipeline
        preprocess.check_replicates(replicate_rows(
            ("WT", "1"), ("WT", "1"), ("WT", "2"), ("IgG", "1"),
        ))

    def test_gap_fails(self):
        with self.assertRaisesRegex(ValueError, "group WT must be numbered 1 to 2, not 1, 3"):
            preprocess.check_replicates(replicate_rows(("WT", "1"), ("WT", "3"), ("IgG", "1")))

    def test_blank_or_non_numeric_fails(self):
        for value in ("", "one", "1.5"):
            with self.assertRaisesRegex(ValueError, "whole-number replicate; check group\\(s\\) IgG"):
                preprocess.check_replicates(replicate_rows(("WT", "1"), ("IgG", value)))


class CheckControls(unittest.TestCase):

    def test_targets_naming_a_control_group_pass(self):
        for use_control in (True, False):
            preprocess.check_controls(
                control_rows(("WT", "IgG"), ("KO", "IgG"), ("IgG", "")), use_control,
            )

    def test_unknown_control_fails(self):
        with self.assertRaisesRegex(ValueError, "Control\\(s\\) IGG are not a group"):
            preprocess.check_controls(control_rows(("WT", "IGG"), ("IgG", "")), True)

    def test_own_group_as_control_fails(self):
        with self.assertRaisesRegex(ValueError, "Group\\(s\\) WT name themselves as their control"):
            preprocess.check_controls(control_rows(("WT", "WT"), ("IgG", "")), False)

    def test_target_without_control_fails_with_use_control(self):
        with self.assertRaisesRegex(ValueError, "Use Control is on.*group\\(s\\) KO$"):
            preprocess.check_controls(control_rows(("WT", "IgG"), ("KO", ""), ("IgG", "")), True)

    def test_no_controls_fails_with_use_control(self):
        with self.assertRaisesRegex(ValueError, "Use Control is on.*group\\(s\\) WT, KO$"):
            preprocess.check_controls(control_rows(("WT", ""), ("KO", "")), True)

    def test_no_controls_pass_without_use_control(self):
        preprocess.check_controls(control_rows(("WT", ""), ("KO", "")), False)


if __name__ == "__main__":
    unittest.main()
