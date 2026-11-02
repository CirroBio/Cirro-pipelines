"""Offline tests for the nf-core/bacass 2.1.0 long-read preprocess script.

preprocess.py imports the Cirro SDK only when it runs, so this is a plain
`python -m unittest` that needs pandas and nothing else.

The samplesheet is the structured input Cirro writes from the run form: a TSV with
CRLF line endings, a header from the form's columns, and blank cells left empty.
"""
import importlib.util
import io
import unittest
from pathlib import Path

import pandas as pd


def _load_preprocess():
    path = Path(__file__).with_name("preprocess.py")
    spec = importlib.util.spec_from_file_location("bacass_long_read_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()

HEADER = "ID\tLongFastQ\r\n"


class ReadSamplesheet(unittest.TestCase):

    def test_values_stay_text(self):
        samplesheet = preprocess.read_samplesheet(io.StringIO(
            HEADER + "01\ts3://bucket/data/01.fastq.gz\r\n"
        ))
        self.assertEqual(samplesheet["ID"].tolist(), ["01"])
        self.assertEqual(samplesheet["LongFastQ"].tolist(), ["s3://bucket/data/01.fastq.gz"])

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


class CheckUniqueIds(unittest.TestCase):

    def test_distinct_ids_pass(self):
        preprocess.check_unique_ids(pd.DataFrame({"ID": ["A", "B"], "LongFastQ": ["a.fq.gz", "b.fq.gz"]}))

    def test_repeated_id_is_named(self):
        samplesheet = pd.DataFrame({
            "ID": ["A", "A", "B"],
            "LongFastQ": ["a1.fq.gz", "a2.fq.gz", "b.fq.gz"],
        })
        with self.assertRaisesRegex(ValueError, "repeated: A$"):
            preprocess.check_unique_ids(samplesheet)


class ToBacassSamplesheet(unittest.TestCase):

    def test_all_six_columns_in_order_with_absent_inputs_as_na(self):
        samplesheet = preprocess.to_bacass_samplesheet(pd.DataFrame({
            "ID": ["A"],
            "LongFastQ": ["s3://bucket/data/a.fastq.gz"],
        }))
        self.assertEqual(samplesheet.columns.tolist(), ["ID", "R1", "R2", "LongFastQ", "Fast5", "GenomeSize"])
        self.assertEqual(
            samplesheet.iloc[0].tolist(),
            ["A", "NA", "NA", "s3://bucket/data/a.fastq.gz", "NA", "NA"],
        )

    def test_blank_cell_becomes_na(self):
        samplesheet = preprocess.to_bacass_samplesheet(pd.DataFrame({"ID": ["A"], "LongFastQ": [""]}))
        self.assertEqual(samplesheet["LongFastQ"].tolist(), ["NA"])


if __name__ == "__main__":
    unittest.main()
