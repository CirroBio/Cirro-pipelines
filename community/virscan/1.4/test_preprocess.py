"""Offline tests for the VirScan 1.4 preprocess script.

preprocess.py imports the Cirro SDK only when it runs, so this is a plain
`python -m unittest` that needs pandas and nothing else.

The sample table is the structured input Cirro writes from the run form: a CSV with
CRLF line endings, a header from the form's columns, blank cells left empty, and an
absolute S3 URI in each FASTQ cell.
"""
import importlib.util
import io
import unittest
from pathlib import Path

import pandas as pd


def _load_preprocess():
    path = Path(__file__).with_name("preprocess.py")
    spec = importlib.util.spec_from_file_location("virscan_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()

HEADER = "sample,sample_name,control_status,fastq_filepath\r\n"
DATA = "s3://bucket/datasets/ds1/data"


def sample_table(*rows) -> pd.DataFrame:
    return preprocess.read_samplesheet(io.StringIO(HEADER + "".join(f"{row}\r\n" for row in rows)))


class ReadSamplesheet(unittest.TestCase):

    def test_values_stay_text(self):
        table = sample_table(f"01,,,{DATA}/01_rep1.R1.fastq.gz")
        self.assertEqual(table["sample"].tolist(), ["01"])
        self.assertEqual(table["control_status"].tolist(), [""])

    def test_no_param_fails_clearly(self):
        for path in (None, ""):
            with self.assertRaisesRegex(ValueError, "No sample table was provided"):
                preprocess.read_samplesheet(path)

    def test_empty_file_fails_clearly(self):
        with self.assertRaisesRegex(ValueError, "The sample table is empty"):
            preprocess.read_samplesheet(io.StringIO(""))

    def test_header_without_rows_fails_clearly(self):
        with self.assertRaisesRegex(ValueError, "The sample table is empty"):
            preprocess.read_samplesheet(io.StringIO(HEADER))


class PublicEpitopesCsv(unittest.TestCase):

    def test_sits_beside_the_library(self):
        self.assertEqual(
            preprocess.public_epitopes_csv("s3://refs/virscan_libraries/v2/library.csv"),
            "s3://refs/virscan_libraries/v2/public_epitopes.csv",
        )

    def test_missing_library_fails_clearly(self):
        for value in (None, ""):
            with self.assertRaisesRegex(ValueError, "No peptide library was selected"):
                preprocess.public_epitopes_csv(value)


class ControlStatus(unittest.TestCase):

    def test_blank_becomes_empirical_and_set_values_are_kept(self):
        table = preprocess.label_controls(sample_table(
            f"A,,,{DATA}/A.fastq.gz",
            f"B,,beads_only,{DATA}/B.fastq.gz",
            f"C,,library,{DATA}/C.fastq.gz",
        ))
        self.assertEqual(table["control_status"].tolist(), ["empirical", "beads_only", "library"])

    def test_two_beads_only_rows_pass(self):
        preprocess.check_beads_only(preprocess.label_controls(sample_table(
            f"A,,,{DATA}/A.fastq.gz",
            f"beads,,beads_only,{DATA}/beads_1.fastq.gz",
            f"beads,,beads_only,{DATA}/beads_2.fastq.gz",
        )))

    def test_a_single_beads_only_row_fails(self):
        table = preprocess.label_controls(sample_table(
            f"A,,,{DATA}/A.fastq.gz",
            f"beads,,beads_only,{DATA}/beads_1.fastq.gz",
        ))
        with self.assertRaisesRegex(ValueError, r"at least 2 rows as beads_only .*\(found 1\)"):
            preprocess.check_beads_only(table)


class NameReplicates(unittest.TestCase):

    def test_blank_name_comes_from_the_file_and_set_names_are_kept(self):
        table = preprocess.name_replicates(sample_table(
            f"A,,,{DATA}/run1/A_rep1.R1.fastq.gz",
            f"A,A-second,,{DATA}/run1/A_rep2.R1.fastq.gz",
            f"B,,,{DATA}/B_rep1.fq.gz",
        ))
        self.assertEqual(table["sample_name"].tolist(), ["A_rep1.R1", "A-second", "B_rep1"])


class SplitReadsPrefix(unittest.TestCase):

    def test_files_in_subfolders_keep_their_subfolder(self):
        reads_prefix, relative = preprocess.split_reads_prefix(pd.Series([
            f"{DATA}/run1/A.fastq.gz",
            f"{DATA}/run2/B.fastq.gz",
            f"{DATA}/C.fastq.gz",
        ]))
        self.assertEqual(reads_prefix, DATA)
        self.assertEqual(relative.tolist(), ["run1/A.fastq.gz", "run2/B.fastq.gz", "C.fastq.gz"])

    def test_files_from_several_datasets_share_the_bucket_prefix(self):
        reads_prefix, relative = preprocess.split_reads_prefix(pd.Series([
            "s3://bucket/datasets/ds1/data/A.fastq.gz",
            "s3://bucket/datasets/ds2/data/B.fastq.gz",
        ]))
        self.assertEqual(reads_prefix, "s3://bucket/datasets")
        self.assertEqual(relative.tolist(), ["ds1/data/A.fastq.gz", "ds2/data/B.fastq.gz"])

    def test_prefix_ends_at_a_folder_boundary(self):
        reads_prefix, relative = preprocess.split_reads_prefix(pd.Series([
            "s3://bucket/data/run1/A.fastq.gz",
            "s3://bucket/data/run10/B.fastq.gz",
        ]))
        self.assertEqual(reads_prefix, "s3://bucket/data")
        self.assertEqual(relative.tolist(), ["run1/A.fastq.gz", "run10/B.fastq.gz"])

    def test_files_at_the_bucket_root_use_the_bucket(self):
        reads_prefix, relative = preprocess.split_reads_prefix(pd.Series([
            "s3://bucket/A.fastq.gz",
            "s3://bucket/data/B.fastq.gz",
        ]))
        self.assertEqual(reads_prefix, "s3://bucket")
        self.assertEqual(relative.tolist(), ["A.fastq.gz", "data/B.fastq.gz"])

    def test_a_single_file_is_relative_to_its_folder(self):
        reads_prefix, relative = preprocess.split_reads_prefix(pd.Series([f"{DATA}/A.fastq.gz"]))
        self.assertEqual(reads_prefix, DATA)
        self.assertEqual(relative.tolist(), ["A.fastq.gz"])

    def test_files_in_different_buckets_fail_naming_them(self):
        with self.assertRaisesRegex(ValueError, "same S3 bucket.*: bucket, shared-s3alias$"):
            preprocess.split_reads_prefix(pd.Series([
                f"{DATA}/A.fastq.gz",
                "s3://shared-s3alias/datasets/ds9/data/B.fastq.gz",
            ]))

    def test_row_without_a_file_fails_clearly(self):
        with self.assertRaisesRegex(ValueError, "needs a FASTQ file"):
            preprocess.split_reads_prefix(pd.Series([f"{DATA}/A.fastq.gz", ""]))


if __name__ == "__main__":
    unittest.main()
