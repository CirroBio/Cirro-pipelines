"""Offline tests for the nf-core/sarek custom genome samplesheet handling.

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
    spec = importlib.util.spec_from_file_location("sarek_custom_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()

HEADER = "patient,sex,sample,lane,fastq_1,fastq_2\r\n"


def samplesheet(rows):
    """A table of (patient, sex, sample, lane) rows, with placeholder FASTQ paths."""
    table = pd.DataFrame(rows, columns=["patient", "sex", "sample", "lane"])
    return table.assign(fastq_1="r1.fastq.gz", fastq_2="r2.fastq.gz")


class ReadSamplesheet(unittest.TestCase):

    def test_values_stay_text(self):
        table = preprocess.read_samplesheet(io.StringIO(
            HEADER + ",,01,001,s3://bucket/data/01_R1.fastq.gz,s3://bucket/data/01_R2.fastq.gz\r\n"
        ))
        self.assertEqual(table["sample"].tolist(), ["01"])
        self.assertEqual(table["lane"].tolist(), ["001"])
        self.assertEqual(table["patient"].tolist(), [""])

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


class NormalizeSex(unittest.TestCase):

    def test_spellings_map_to_sarek_values(self):
        cases = {
            "male": "XY", " M ": "XY", "XY": "XY",
            "Female": "XX", "f": "XX", "xx": "XX",
            "": "NA", "NA": "NA", "unknown": "NA", "other": "NA",
        }
        for value, expected in cases.items():
            self.assertEqual(preprocess.normalize_sex(value), expected, value)


class FillSamplesheet(unittest.TestCase):

    def test_blank_patient_takes_the_sample_name(self):
        filled = preprocess.fill_samplesheet(samplesheet([
            ("", "", "S1", "1"),
            ("P2", "", "S2", "1"),
        ]))
        self.assertEqual(filled["patient"].tolist(), ["S1", "P2"])

    def test_sex_is_normalized(self):
        filled = preprocess.fill_samplesheet(samplesheet([
            ("", "male", "S1", "1"),
            ("", "", "S2", "1"),
        ]))
        self.assertEqual(filled["sex"].tolist(), ["XY", "NA"])

    def test_blank_lanes_are_numbered_within_each_sample(self):
        filled = preprocess.fill_samplesheet(samplesheet([
            ("", "", "S1", ""),
            ("", "", "S1", ""),
            ("", "", "S2", ""),
        ]))
        self.assertEqual(filled["lane"].tolist(), ["1", "2", "1"])

    def test_blank_lanes_skip_lanes_the_sample_already_uses(self):
        filled = preprocess.fill_samplesheet(samplesheet([
            ("", "", "S1", "1"),
            ("", "", "S1", ""),
            ("", "", "S1", "L002"),
            ("", "", "S1", ""),
        ]))
        self.assertEqual(filled["lane"].tolist(), ["1", "2", "L002", "3"])

    def test_columns_keep_their_order(self):
        table = samplesheet([("P1", "XX", "S1", "1")])
        filled = preprocess.fill_samplesheet(table)
        self.assertEqual(filled.columns.tolist(), table.columns.tolist())
        self.assertTrue(filled.equals(table))


class CheckOnePatientPerSample(unittest.TestCase):

    def test_consistent_patients_pass(self):
        preprocess.check_one_patient_per_sample(samplesheet([
            ("P1", "NA", "S1", "1"),
            ("P1", "NA", "S1", "2"),
            ("P1", "NA", "S2", "1"),
        ]))

    def test_sample_with_two_patients_fails(self):
        table = samplesheet([
            ("P1", "NA", "S1", "1"),
            ("S1", "NA", "S1", "2"),
        ])
        with self.assertRaisesRegex(ValueError, r"these have several: S1 \(P1, S1\)"):
            preprocess.check_one_patient_per_sample(table)


class CheckUniqueLanes(unittest.TestCase):

    def test_distinct_lanes_pass(self):
        preprocess.check_unique_lanes(samplesheet([
            ("P1", "NA", "S1", "1"),
            ("P1", "NA", "S1", "2"),
            ("P1", "NA", "S2", "1"),
        ]))

    def test_repeated_lane_fails(self):
        table = samplesheet([
            ("P1", "NA", "S1", "1"),
            ("P1", "NA", "S1", "1"),
        ])
        with self.assertRaisesRegex(ValueError, "patient P1, sample S1, lane 1"):
            preprocess.check_unique_lanes(table)


if __name__ == "__main__":
    unittest.main()
