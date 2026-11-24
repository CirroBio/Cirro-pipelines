"""Offline tests for the NGSCheckMate preprocess script.

preprocess.py imports the Cirro SDK only when it runs, so this is a plain
`python -m unittest` that needs pandas and nothing else.

The samplesheet is the structured input Cirro writes from the run form: a CSV with
CRLF line endings, a header from the form's columns, and one row per FASTQ pair.
"""
import importlib.util
import io
import unittest
from pathlib import Path


def _load_preprocess():
    path = Path(__file__).with_name("preprocess.py")
    spec = importlib.util.spec_from_file_location("ngscheckmate_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()

HEADER = "sample,fastq_1,fastq_2\r\n"
LANE_1 = "01,s3://bucket/data/01_L001_R1_001.fastq.gz,s3://bucket/data/01_L001_R2_001.fastq.gz\r\n"
LANE_2 = "01,s3://bucket/data/01_L002_R1_001.fastq.gz,s3://bucket/data/01_L002_R2_001.fastq.gz\r\n"
OTHER = "02,s3://bucket/data/02_R1_001.fastq.gz,s3://bucket/data/02_R2_001.fastq.gz\r\n"


class ReadSamplesheet(unittest.TestCase):

    def test_values_stay_text(self):
        samplesheet = preprocess.read_samplesheet(io.StringIO(HEADER + LANE_1))
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


class FirstPairPerSample(unittest.TestCase):

    def test_later_lanes_of_a_sample_are_dropped(self):
        samplesheet = preprocess.first_pair_per_sample(
            preprocess.read_samplesheet(io.StringIO(HEADER + LANE_1 + LANE_2 + OTHER))
        )
        self.assertEqual(samplesheet["sample"].tolist(), ["01", "02"])
        self.assertEqual(samplesheet["fastq_1"].tolist()[0], "s3://bucket/data/01_L001_R1_001.fastq.gz")


if __name__ == "__main__":
    unittest.main()
