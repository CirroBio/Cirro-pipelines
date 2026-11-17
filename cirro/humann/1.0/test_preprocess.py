"""Offline tests for the HUMAnN 1.0 preprocess script.

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
    spec = importlib.util.spec_from_file_location("humann_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()

HEADER = "sample,fastq_1,fastq_2\r\n"
REFERENCES = "s3://refs-bucket"


class AbsoluteReference(unittest.TestCase):

    def test_relative_path_is_put_under_the_references_bucket(self):
        self.assertEqual(
            preprocess.absolute_reference("humann/uniref/uniref90_201901b_full.dmnd", REFERENCES),
            "s3://refs-bucket/humann/uniref/uniref90_201901b_full.dmnd",
        )

    def test_absolute_path_is_left_alone(self):
        value = "s3://other-bucket/humann/uniref/uniref90_201901b_full.dmnd"
        self.assertEqual(preprocess.absolute_reference(value, REFERENCES), value)

    def test_missing_value_is_left_alone(self):
        for value in (None, ""):
            self.assertEqual(preprocess.absolute_reference(value, REFERENCES), value)


class ReadSamplesheet(unittest.TestCase):

    def test_single_end_row_keeps_fastq_2_blank(self):
        samplesheet = preprocess.read_samplesheet(io.StringIO(
            HEADER + "01,s3://bucket/data/01_R1.fastq.gz,\r\n"
        ))
        self.assertEqual(samplesheet["sample"].tolist(), ["01"])
        # nf-humann runs a row single-end when fastq_2 is empty
        self.assertEqual(samplesheet["fastq_2"].tolist(), [""])

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


class CheckUniqueSamples(unittest.TestCase):

    def test_distinct_samples_pass(self):
        preprocess.check_unique_samples(pd.DataFrame({
            "sample": ["A", "B"],
            "fastq_1": ["a_R1.fastq.gz", "b_R1.fastq.gz"],
            "fastq_2": ["a_R2.fastq.gz", ""],
        }))

    def test_repeated_sample_is_named(self):
        samplesheet = pd.DataFrame({
            "sample": ["A", "A", "B"],
            "fastq_1": ["a_L1_R1.fastq.gz", "a_L2_R1.fastq.gz", "b_R1.fastq.gz"],
            "fastq_2": ["a_L1_R2.fastq.gz", "a_L2_R2.fastq.gz", "b_R2.fastq.gz"],
        })
        with self.assertRaisesRegex(ValueError, "repeated: A$"):
            preprocess.check_unique_samples(samplesheet)


if __name__ == "__main__":
    unittest.main()
