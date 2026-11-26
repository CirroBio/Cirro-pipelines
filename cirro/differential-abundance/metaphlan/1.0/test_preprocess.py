"""Offline tests for the differential abundance (MetaPhlAn) preprocess script.

preprocess.py imports the Cirro SDK only when it runs, so this is a plain
`python -m unittest` that needs pandas and nothing else.

The samplesheet is the structured input Cirro writes from the run form: a CSV with
CRLF line endings, a header from the form's columns, and one row per sample. The
metadata is a CSV the user picks, with any columns.
"""
import importlib.util
import io
import unittest
from pathlib import Path

import pandas as pd


def _load_preprocess():
    path = Path(__file__).with_name("preprocess.py")
    spec = importlib.util.spec_from_file_location("differential_abundance_metaphlan_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()

HEADER = "sample,file\r\n"
FILE_1 = "s3://bucket/datasets/abc/data/mpl/S1.metaphlan"
FILE_2 = "s3://bucket/datasets/abc/data/mpl/S2.metaphlan"
VALID = {"formula": "group", "fdr_cutoff": 0.05}


def samplesheet():
    return preprocess.read_samplesheet(io.StringIO(HEADER + f"01,{FILE_1}\r\n02,{FILE_2}\r\n"))


def metadata(text: str) -> pd.DataFrame:
    return pd.read_csv(io.StringIO(text), dtype=str, keep_default_na=False)


class CheckParams(unittest.TestCase):

    def test_valid_params_pass(self):
        preprocess.check_params(VALID)

    def test_missing_formula_fails(self):
        for formula in (None, ""):
            with self.assertRaisesRegex(ValueError, "Must provide formula"):
                preprocess.check_params({**VALID, "formula": formula})

    def test_fdr_cutoff_out_of_range_fails(self):
        for cutoff in (-0.1, 1, 1.5):
            with self.assertRaisesRegex(ValueError, "fdr_cutoff"):
                preprocess.check_params({**VALID, "fdr_cutoff": cutoff})


class ReadSamplesheet(unittest.TestCase):

    def test_sample_names_are_read_as_text(self):
        self.assertEqual(samplesheet()["sample"].tolist(), ["01", "02"])

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


class AddMetadata(unittest.TestCase):

    def test_metadata_columns_are_added_to_each_sample(self):
        merged = preprocess.add_metadata(
            samplesheet(),
            metadata("sample,group,age\n02,b,40\n01,a,\n03,c,9\n")
        )
        self.assertEqual(merged.columns.tolist(), ["sample", "file", "group", "age"])
        self.assertEqual(merged.values.tolist(), [["01", FILE_1, "a", ""], ["02", FILE_2, "b", "40"]])

    def test_metadata_without_sample_column_fails(self):
        with self.assertRaisesRegex(ValueError, "must have a 'sample' column"):
            preprocess.add_metadata(samplesheet(), metadata("id,group\n01,a\n02,b\n"))

    def test_metadata_with_a_file_column_fails(self):
        with self.assertRaisesRegex(ValueError, "cannot have a column named file"):
            preprocess.add_metadata(samplesheet(), metadata("sample,file,group\n01,x,a\n02,y,b\n"))

    def test_sample_missing_from_metadata_fails(self):
        with self.assertRaisesRegex(ValueError, "no row in the metadata table: 02"):
            preprocess.add_metadata(samplesheet(), metadata("sample,group\n01,a\n"))

    def test_repeated_sample_fails(self):
        with self.assertRaisesRegex(ValueError, "metadata table, but these are repeated: 01"):
            preprocess.add_metadata(samplesheet(), metadata("sample,group\n01,a\n01,b\n02,c\n"))
        repeated = preprocess.read_samplesheet(io.StringIO(HEADER + f"01,{FILE_1}\r\n01,{FILE_2}\r\n"))
        with self.assertRaisesRegex(ValueError, "samplesheet, but these are repeated: 01"):
            preprocess.add_metadata(repeated, metadata("sample,group\n01,a\n"))


if __name__ == "__main__":
    unittest.main()
