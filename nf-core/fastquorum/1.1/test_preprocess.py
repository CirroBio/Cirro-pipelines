"""Offline tests for the nf-core/fastquorum 1.1 samplesheet handling.

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
    spec = importlib.util.spec_from_file_location("fastquorum_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()

HEADER = "sample,fastq_1,fastq_2,fastq_3,fastq_4,read_structure\r\n"


def read_structure_rows(*rows):
    return pd.DataFrame(rows, columns=["sample", "read_structure"])


class ReadSamplesheet(unittest.TestCase):

    def test_values_stay_text(self):
        samplesheet = preprocess.read_samplesheet(io.StringIO(
            HEADER + "01,s3://bucket/data/01_R1.fastq.gz,s3://bucket/data/01_R2.fastq.gz,,,6M+T 6M+T\r\n"
        ))
        self.assertEqual(samplesheet["sample"].tolist(), ["01"])
        self.assertEqual(samplesheet["fastq_3"].tolist(), [""])
        self.assertEqual(samplesheet["read_structure"].tolist(), ["6M+T 6M+T"])

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


class FillReadStructure(unittest.TestCase):

    def test_blank_rows_take_the_form_value_and_set_rows_are_kept(self):
        samplesheet = read_structure_rows(("A", ""), ("B", "+T +T 8M"), ("C", " "))
        filled = preprocess.fill_read_structure(samplesheet, "6M+T 6M+T")
        self.assertEqual(filled["read_structure"].tolist(), ["6M+T 6M+T", "+T +T 8M", "6M+T 6M+T"])

    def test_every_row_set_needs_no_form_value(self):
        samplesheet = read_structure_rows(("A", "6M+T 6M+T"), ("B", "+T +T 8M"))
        for form_value in (None, ""):
            filled = preprocess.fill_read_structure(samplesheet, form_value)
            self.assertTrue(filled.equals(samplesheet))

    def test_blank_row_without_form_value_fails(self):
        samplesheet = read_structure_rows(("A", "6M+T 6M+T"), ("B", ""), ("B", ""))
        for form_value in (None, ""):
            with self.assertRaisesRegex(ValueError, "missing for sample\\(s\\) B$"):
                preprocess.fill_read_structure(samplesheet, form_value)


if __name__ == "__main__":
    unittest.main()
