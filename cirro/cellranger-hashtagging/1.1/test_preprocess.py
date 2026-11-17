"""Offline tests for the CellRanger hashtagging 1.1 preprocess script.

preprocess.py imports the Cirro SDK only when it runs, so this is a plain
`python -m unittest` that needs pandas and nothing else.

The samplesheet and the sample hashtag table are structured inputs Cirro writes from
the run form: CSVs with CRLF line endings and a header from the form's columns.
"""
import importlib.util
import io
import unittest
from pathlib import Path


def _load_preprocess():
    path = Path(__file__).with_name("preprocess.py")
    spec = importlib.util.spec_from_file_location("cellranger_hashtagging_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()

HEADER = "sample,feature_types\r\n"
GEX = "01_GEX,Gene Expression\r\n"
HTO = "01_HTO,Multiplexing Capture\r\n"
REFERENCES = "s3://refs-bucket"


class ReferenceDir(unittest.TestCase):

    def test_relative_path_is_put_under_the_references_bucket_as_a_folder(self):
        self.assertEqual(
            preprocess.reference_dir("cellranger/refdata-gex-GRCh38-2024-A", REFERENCES),
            "s3://refs-bucket/cellranger/refdata-gex-GRCh38-2024-A/",
        )

    def test_absolute_path_only_gains_the_trailing_separator(self):
        self.assertEqual(
            preprocess.reference_dir("s3://other/cellranger/refdata-gex-GRCh38-2024-A", REFERENCES),
            "s3://other/cellranger/refdata-gex-GRCh38-2024-A/",
        )
        self.assertEqual(
            preprocess.reference_dir("s3://other/refdata/", REFERENCES),
            "s3://other/refdata/",
        )

    def test_missing_value_is_left_alone(self):
        for value in (None, ""):
            self.assertEqual(preprocess.reference_dir(value, REFERENCES), value)


class ReadTable(unittest.TestCase):

    def test_values_stay_text(self):
        samplesheet = preprocess.read_table(io.StringIO(HEADER + GEX + HTO), "samplesheet", "rows")
        self.assertEqual(samplesheet.columns.tolist(), ["sample", "feature_types"])
        self.assertEqual(samplesheet["sample"].tolist(), ["01_GEX", "01_HTO"])

    def test_no_param_fails_clearly(self):
        for path in (None, ""):
            with self.assertRaisesRegex(ValueError, "No sample hashtag table was provided: add rows"):
                preprocess.read_table(path, "sample hashtag table", "rows")

    def test_empty_file_fails_clearly(self):
        with self.assertRaisesRegex(ValueError, "The samplesheet is empty"):
            preprocess.read_table(io.StringIO(""), "samplesheet", "rows")

    def test_header_without_rows_fails_clearly(self):
        with self.assertRaisesRegex(ValueError, "The sample hashtag table is empty"):
            preprocess.read_table(io.StringIO("sample_id,cmo_ids\r\n"), "sample hashtag table", "rows")


class DropRepeatedRows(unittest.TestCase):

    def test_rows_repeated_per_lane_become_one_library(self):
        samplesheet = preprocess.read_table(
            io.StringIO(HEADER + GEX + GEX + HTO + HTO), "samplesheet", "rows"
        )
        libraries = preprocess.drop_repeated_rows(samplesheet)
        self.assertEqual(
            libraries.values.tolist(),
            [["01_GEX", "Gene Expression"], ["01_HTO", "Multiplexing Capture"]],
        )

    def test_distinct_rows_are_kept(self):
        samplesheet = preprocess.read_table(io.StringIO(HEADER + GEX + HTO), "samplesheet", "rows")
        self.assertEqual(len(preprocess.drop_repeated_rows(samplesheet)), 2)


class UnsetUmiThresholds(unittest.TestCase):

    def test_blank_thresholds_are_unset(self):
        params = {"cellranger_version": "7.0.1", "emptydrops_minimum_umis": None}
        self.assertEqual(
            preprocess.unset_umi_thresholds(params),
            ["emptydrops_minimum_umis", "global_minimum_umis"],
        )

    def test_thresholds_set_on_a_new_enough_version_are_kept(self):
        params = {"cellranger_version": "8.0.0", "emptydrops_minimum_umis": 100, "global_minimum_umis": 0}
        self.assertEqual(preprocess.unset_umi_thresholds(params), [])

    def test_threshold_newer_than_the_version_fails_clearly(self):
        params = {"cellranger_version": "7.2.0", "global_minimum_umis": 50}
        with self.assertRaisesRegex(ValueError, r"global_minimum_umis requires CellRanger 8\.0 or later \(selected: 7\.2\.0\)"):
            preprocess.unset_umi_thresholds(params)


if __name__ == "__main__":
    unittest.main()
