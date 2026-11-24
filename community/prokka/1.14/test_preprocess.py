"""Offline tests for the Prokka preprocess script.

preprocess.py imports the Cirro SDK only when it runs, so this is a plain
`python -m unittest` that needs pandas and nothing else.

The genome table is the structured input Cirro writes from the run form: a CSV with
CRLF line endings, a header from the form's columns, and one row per genome.
"""
import importlib.util
import io
import unittest
from pathlib import Path


def _load_preprocess():
    path = Path(__file__).with_name("preprocess.py")
    spec = importlib.util.spec_from_file_location("prokka_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()

HEADER = "name,fasta\r\n"
ROW_01 = "01,s3://bucket/data/01.fna.gz\r\n"
ROW_02 = "02,s3://bucket/data/02.fasta\r\n"


class ReadSamplesheet(unittest.TestCase):

    def test_values_stay_text(self):
        samplesheet = preprocess.read_samplesheet(io.StringIO(HEADER + ROW_01))
        self.assertEqual(samplesheet["name"].tolist(), ["01"])

    def test_no_param_fails_clearly(self):
        for path in (None, ""):
            with self.assertRaisesRegex(ValueError, "No genomes were provided"):
                preprocess.read_samplesheet(path)

    def test_empty_file_fails_clearly(self):
        with self.assertRaisesRegex(ValueError, "The genome table is empty"):
            preprocess.read_samplesheet(io.StringIO(""))

    def test_header_without_rows_fails_clearly(self):
        with self.assertRaisesRegex(ValueError, "The genome table is empty"):
            preprocess.read_samplesheet(io.StringIO(HEADER))


class CheckUniqueNames(unittest.TestCase):

    def test_distinct_names_pass(self):
        preprocess.check_unique_names(preprocess.read_samplesheet(io.StringIO(HEADER + ROW_01 + ROW_02)))

    def test_repeated_name_fails_naming_it(self):
        samplesheet = preprocess.read_samplesheet(io.StringIO(
            HEADER + ROW_01 + "01,s3://bucket/data/01b.fna.gz\r\n"
        ))
        with self.assertRaisesRegex(ValueError, "Only one FASTA is allowed per genome.*: 01$"):
            preprocess.check_unique_names(samplesheet)


if __name__ == "__main__":
    unittest.main()
