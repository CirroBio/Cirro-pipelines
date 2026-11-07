"""Offline tests for the flow cytometry QC preprocess script.

preprocess.py imports the Cirro SDK only when it runs, so this is a plain
`python -m unittest` that needs pandas and nothing else.

The samplesheet is the structured input Cirro writes from the run form: a CSV with
CRLF line endings, a header from the form's columns, and one row per FCS file.
"""
import importlib.util
import io
import unittest
from pathlib import Path


def _load_preprocess():
    path = Path(__file__).with_name("preprocess.py")
    spec = importlib.util.spec_from_file_location("flow_cytometry_qc_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()

HEADER = "sample,fcs\r\n"
ROW_01 = "01,s3://bucket/data/plate1/01.fcs\r\n"
ROW_02 = "02,s3://bucket/data/plate1/02.FCS\r\n"


class ReadSamplesheet(unittest.TestCase):

    def test_values_stay_text(self):
        samplesheet = preprocess.read_samplesheet(io.StringIO(HEADER + ROW_01))
        self.assertEqual(samplesheet.columns.tolist(), ["sample", "fcs"])
        self.assertEqual(samplesheet["sample"].tolist(), ["01"])

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


class CheckUnique(unittest.TestCase):

    def test_distinct_rows_pass(self):
        preprocess.check_unique(preprocess.read_samplesheet(io.StringIO(HEADER + ROW_01 + ROW_02)))

    def test_repeated_sample_fails_naming_it(self):
        samplesheet = preprocess.read_samplesheet(io.StringIO(
            HEADER + ROW_01 + "01,s3://bucket/data/plate2/01.fcs\r\n" + ROW_02
        ))
        with self.assertRaisesRegex(ValueError, "Each sample name can be listed only once.*: 01$"):
            preprocess.check_unique(samplesheet)

    def test_repeated_fcs_file_fails_naming_it(self):
        samplesheet = preprocess.read_samplesheet(io.StringIO(
            HEADER + ROW_01 + "01b,s3://bucket/data/plate1/01.fcs\r\n"
        ))
        with self.assertRaisesRegex(
            ValueError, "Each FCS file can be listed only once.*: s3://bucket/data/plate1/01.fcs$"
        ):
            preprocess.check_unique(samplesheet)


if __name__ == "__main__":
    unittest.main()
