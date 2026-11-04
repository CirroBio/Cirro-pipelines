"""Offline tests for the KneadData preprocess script.

preprocess.py imports the Cirro SDK only when it runs, so this is a plain
`python -m unittest` that needs pandas and nothing else.

The samplesheet is the structured input Cirro writes from the run form: a CSV with
CRLF line endings, a header from the form's columns, and one row per sample.
"""
import importlib.util
import io
import unittest
from pathlib import Path


def _load_preprocess():
    path = Path(__file__).with_name("preprocess.py")
    spec = importlib.util.spec_from_file_location("kneaddata_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()

HEADER = "sample,fastq_1,fastq_2\r\n"
ROW_01 = "01,s3://bucket/data/01_R1.fastq.gz,s3://bucket/data/01_R2.fastq.gz\r\n"
ROW_02 = "02,s3://bucket/data/02_R1.fastq.gz,s3://bucket/data/02_R2.fastq.gz\r\n"
REFERENCES = "s3://refs-bucket"


class AbsoluteReference(unittest.TestCase):

    def test_relative_path_is_put_under_the_references_bucket(self):
        self.assertEqual(
            preprocess.absolute_reference("kneaddata/human_genome/bowtie2/", REFERENCES),
            "s3://refs-bucket/kneaddata/human_genome/bowtie2/",
        )

    def test_absolute_path_is_left_alone(self):
        value = "s3://other-bucket/kneaddata/human_genome/bowtie2/"
        self.assertEqual(preprocess.absolute_reference(value, REFERENCES), value)

    def test_missing_value_is_left_alone(self):
        for value in (None, ""):
            self.assertEqual(preprocess.absolute_reference(value, REFERENCES), value)


class ReadSamplesheet(unittest.TestCase):

    def test_values_stay_text(self):
        samplesheet = preprocess.read_samplesheet(io.StringIO(HEADER + ROW_01))
        self.assertEqual(samplesheet.columns.tolist(), ["sample", "fastq_1", "fastq_2"])
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


class CheckUniqueSamples(unittest.TestCase):

    def test_distinct_samples_pass(self):
        preprocess.check_unique_samples(
            preprocess.read_samplesheet(io.StringIO(HEADER + ROW_01 + ROW_02))
        )

    def test_sample_listed_per_lane_fails_naming_it(self):
        samplesheet = preprocess.read_samplesheet(io.StringIO(
            HEADER
            + "01,s3://bucket/data/01_S1_L001_R1_001.fastq.gz,s3://bucket/data/01_S1_L001_R2_001.fastq.gz\r\n"
            + "01,s3://bucket/data/01_S1_L002_R1_001.fastq.gz,s3://bucket/data/01_S1_L002_R2_001.fastq.gz\r\n"
            + ROW_02
        ))
        with self.assertRaisesRegex(ValueError, "Listed more than once: 01$"):
            preprocess.check_unique_samples(samplesheet)


if __name__ == "__main__":
    unittest.main()
