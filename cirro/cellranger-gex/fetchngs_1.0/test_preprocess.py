"""Offline tests for the Cell Ranger GEX (Public Data) preprocess script.

preprocess.py imports the Cirro SDK only when it runs, so this is a plain
`python -m unittest` that needs pandas and nothing else.

The samplesheet is the structured input Cirro writes from the run form: a CSV with
CRLF line endings, a header from the form's columns, and one row per SRA run.
"""
import importlib.util
import io
import unittest
from pathlib import Path


def _load_preprocess():
    path = Path(__file__).with_name("preprocess.py")
    spec = importlib.util.spec_from_file_location("cellranger_gex_fetchngs_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()

HEADER = "sample,fastq_1,fastq_2\r\n"
REFERENCES = "s3://refs-bucket"


class AbsoluteReference(unittest.TestCase):

    def test_relative_path_is_put_under_the_references_bucket(self):
        self.assertEqual(
            preprocess.absolute_reference("cellranger/refdata-gex-GRCh38-2024-A/", REFERENCES),
            "s3://refs-bucket/cellranger/refdata-gex-GRCh38-2024-A/",
        )

    def test_absolute_path_is_left_alone(self):
        value = "s3://other-bucket/cellranger/refdata-gex-GRCh38-2024-A/"
        self.assertEqual(preprocess.absolute_reference(value, REFERENCES), value)

    def test_missing_value_is_left_alone(self):
        for value in (None, ""):
            self.assertEqual(preprocess.absolute_reference(value, REFERENCES), value)


class ReadSamplesheet(unittest.TestCase):

    def test_run_rows_are_read_as_text(self):
        samplesheet = preprocess.read_samplesheet(io.StringIO(
            HEADER
            + "SRX3815835,s3://bucket/data/fastq/SRX3815835_SRR6860803_1.fastq.gz,s3://bucket/data/fastq/SRX3815835_SRR6860803_2.fastq.gz\r\n"
            + "SRX3815835,s3://bucket/data/fastq/SRX3815835_SRR6860804_1.fastq.gz,s3://bucket/data/fastq/SRX3815835_SRR6860804_2.fastq.gz\r\n"
        ))
        self.assertEqual(samplesheet["sample"].tolist(), ["SRX3815835", "SRX3815835"])

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


class OrderReads(unittest.TestCase):

    SAMPLESHEET = preprocess.read_samplesheet(io.StringIO(
        HEADER
        + "SRX1,s3://b/SRX1_SRR1_1.fastq.gz,s3://b/SRX1_SRR1_2.fastq.gz\r\n"
        + "SRX1,s3://b/SRX1_SRR2_1.fastq.gz,s3://b/SRX1_SRR2_2.fastq.gz\r\n"
    ))

    def test_smaller_fastq_1_is_kept(self):
        sizes = {
            "s3://b/SRX1_SRR1_1.fastq.gz": 10, "s3://b/SRX1_SRR1_2.fastq.gz": 40,
            "s3://b/SRX1_SRR2_1.fastq.gz": 12, "s3://b/SRX1_SRR2_2.fastq.gz": 45,
        }
        samplesheet, swapped = preprocess.order_reads(self.SAMPLESHEET, sizes)
        self.assertFalse(swapped)
        self.assertEqual(samplesheet["fastq_1"].tolist(), ["s3://b/SRX1_SRR1_1.fastq.gz", "s3://b/SRX1_SRR2_1.fastq.gz"])

    def test_larger_fastq_1_is_swapped_with_columns_in_form_order(self):
        sizes = {
            "s3://b/SRX1_SRR1_1.fastq.gz": 40, "s3://b/SRX1_SRR1_2.fastq.gz": 10,
            "s3://b/SRX1_SRR2_1.fastq.gz": 45, "s3://b/SRX1_SRR2_2.fastq.gz": 12,
        }
        samplesheet, swapped = preprocess.order_reads(self.SAMPLESHEET, sizes)
        self.assertTrue(swapped)
        self.assertEqual(list(samplesheet.columns), ["sample", "fastq_1", "fastq_2"])
        self.assertEqual(samplesheet["fastq_1"].tolist(), ["s3://b/SRX1_SRR1_2.fastq.gz", "s3://b/SRX1_SRR2_2.fastq.gz"])
        self.assertEqual(samplesheet["fastq_2"].tolist(), ["s3://b/SRX1_SRR1_1.fastq.gz", "s3://b/SRX1_SRR2_1.fastq.gz"])


if __name__ == "__main__":
    unittest.main()
