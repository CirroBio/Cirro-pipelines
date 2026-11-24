"""Offline tests for the gig-map align_reads preprocess script.

preprocess.py imports the Cirro SDK only when it runs, so this is a plain
`python -m unittest` that needs pandas and nothing else.

The samplesheet is the structured input Cirro writes from the run form: a CSV with
CRLF line endings, a header from the form's columns, and blank cells left empty.
"""
import importlib.util
import io
import unittest
from pathlib import Path


def _load_preprocess():
    path = Path(__file__).with_name("preprocess.py")
    spec = importlib.util.spec_from_file_location("gig_map_align_reads_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()

HEADER = "sample,fastq_1,fastq_2\r\n"
PAIRED = "01,s3://bucket/data/01_R1.fastq.gz,s3://bucket/data/01_R2.fastq.gz\r\n"
SINGLE = "02,s3://bucket/data/02_R1.fastq.gz,\r\n"


class ReadSamplesheet(unittest.TestCase):

    def test_values_stay_text(self):
        samplesheet = preprocess.read_samplesheet(io.StringIO(HEADER + SINGLE))
        self.assertEqual(samplesheet["sample"].tolist(), ["02"])
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


class DropBlankFastq2(unittest.TestCase):

    def test_paired_end_is_kept(self):
        samplesheet = preprocess.read_samplesheet(io.StringIO(HEADER + PAIRED))
        self.assertEqual(
            preprocess.drop_blank_fastq_2(samplesheet).columns.tolist(),
            ["sample", "fastq_1", "fastq_2"],
        )

    def test_single_end_loses_the_r2_column(self):
        samplesheet = preprocess.read_samplesheet(io.StringIO(HEADER + SINGLE))
        self.assertEqual(
            preprocess.drop_blank_fastq_2(samplesheet).columns.tolist(),
            ["sample", "fastq_1"],
        )

    def test_mixed_rows_fail_clearly(self):
        samplesheet = preprocess.read_samplesheet(io.StringIO(HEADER + PAIRED + SINGLE))
        with self.assertRaisesRegex(ValueError, r"no FASTQ \(R2\) file \(02\)"):
            preprocess.drop_blank_fastq_2(samplesheet)


class CatalogGenes(unittest.TestCase):

    def test_folder_pickers_trailing_slash_is_not_doubled(self):
        self.assertEqual(
            preprocess.catalog_genes("s3://bucket/datasets/catalog/data/"),
            "s3://bucket/datasets/catalog/data/centroids.faa.gz",
        )


if __name__ == "__main__":
    unittest.main()
