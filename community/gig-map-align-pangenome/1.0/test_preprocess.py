"""Offline tests for the gig-map align_pangenome preprocess script.

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
    spec = importlib.util.spec_from_file_location("gig_map_align_pangenome_preprocess", path)
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


class PangenomePaths(unittest.TestCase):

    def test_folder_pickers_trailing_slash_is_not_doubled(self):
        paths = preprocess.pangenome_paths("s3://bucket/datasets/pg/data/")
        self.assertEqual(paths["genes"], "s3://bucket/datasets/pg/data/gene_catalog/centroids.faa.gz")
        self.assertEqual(paths["gene_bins"], "s3://bucket/datasets/pg/data/bin_pangenome/gene_bins.csv")
        self.assertEqual(paths["group_profile"], "s3://bucket/datasets/pg/data/bin_pangenome/group_profile.csv")
        self.assertEqual(paths["genome_groups"], "s3://bucket/datasets/pg/data/bin_pangenome/genome_groups.csv")


class CheckMetadata(unittest.TestCase):

    METADATA = pd.DataFrame({"sample": ["01"], "treatment": ["A"], "batch": ["1"]})

    def test_formula_columns_present(self):
        preprocess.check_metadata(self.METADATA, "treatment + batch")

    def test_no_formula_needs_only_a_sample_column(self):
        preprocess.check_metadata(self.METADATA[["sample"]], None)

    def test_missing_sample_column_fails_clearly(self):
        with self.assertRaisesRegex(ValueError, "no column sample"):
            preprocess.check_metadata(self.METADATA.drop(columns="sample"), None)

    def test_missing_formula_column_fails_clearly(self):
        with self.assertRaisesRegex(ValueError, "no column group"):
            preprocess.check_metadata(self.METADATA, "group")


if __name__ == "__main__":
    unittest.main()
