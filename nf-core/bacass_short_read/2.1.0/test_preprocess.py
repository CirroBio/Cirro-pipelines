"""Offline tests for the nf-core/bacass 2.4.0 short-read preprocess script.

preprocess.py imports the Cirro SDK only when it runs, so this is a plain
`python -m unittest` that needs pandas and nothing else.

The samplesheet is the structured input Cirro writes from the run form: a TSV with
CRLF line endings, a header from the form's columns, and one row per pair of FASTQ files.
"""
import importlib.util
import io
import unittest
from pathlib import Path

import pandas as pd


def _load_preprocess():
    path = Path(__file__).with_name("preprocess.py")
    spec = importlib.util.spec_from_file_location("bacass_short_read_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()

HEADER = "ID\tR1\tR2\r\n"


class ReadSamplesheet(unittest.TestCase):

    def test_values_stay_text(self):
        samplesheet = preprocess.read_samplesheet(io.StringIO(
            HEADER + "01\ts3://bucket/data/01_R1.fastq.gz\ts3://bucket/data/01_R2.fastq.gz\r\n"
        ))
        self.assertEqual(samplesheet["ID"].tolist(), ["01"])
        self.assertEqual(samplesheet["R2"].tolist(), ["s3://bucket/data/01_R2.fastq.gz"])

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


class ToBacassSamplesheet(unittest.TestCase):

    def test_all_six_columns_in_order_with_absent_inputs_as_na(self):
        samplesheet = preprocess.to_bacass_samplesheet(pd.DataFrame({
            "ID": ["A"],
            "R1": ["s3://bucket/data/a_R1.fastq.gz"],
            "R2": ["s3://bucket/data/a_R2.fastq.gz"],
        }))
        self.assertEqual(samplesheet.columns.tolist(), ["ID", "R1", "R2", "LongFastQ", "Fast5", "GenomeSize"])
        self.assertEqual(
            samplesheet.iloc[0].tolist(),
            ["A", "s3://bucket/data/a_R1.fastq.gz", "s3://bucket/data/a_R2.fastq.gz", "NA", "NA", "NA"],
        )

    def test_blank_cell_becomes_na(self):
        samplesheet = preprocess.to_bacass_samplesheet(pd.DataFrame({"ID": ["A"], "R1": ["a_R1.fq.gz"], "R2": [""]}))
        self.assertEqual(samplesheet["R2"].tolist(), ["NA"])


if __name__ == "__main__":
    unittest.main()
