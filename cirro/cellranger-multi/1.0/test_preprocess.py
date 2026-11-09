"""Offline tests for the CellRanger multi preprocess script.

preprocess.py imports the Cirro SDK only when it runs, so this is a plain
`python -m unittest` that needs pandas and nothing else.

The samplesheet and the probe barcode table are structured inputs Cirro writes from
the run form: CSVs with CRLF line endings and a header from the form's columns.
"""
import importlib.util
import io
import unittest
from pathlib import Path


def _load_preprocess():
    path = Path(__file__).with_name("preprocess.py")
    spec = importlib.util.spec_from_file_location("cellranger_multi_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()

HEADER = "sample,grouping,feature_types\r\n"
GEX = "01_GEX,01,Gene Expression\r\n"
VDJ = "01_VDJ,01,VDJ-B\r\n"
REFERENCES = "s3://refs-bucket"


class ReferenceDir(unittest.TestCase):

    def test_relative_path_is_put_under_the_references_bucket_as_a_folder(self):
        self.assertEqual(
            preprocess.reference_dir("cellranger/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0", REFERENCES),
            "s3://refs-bucket/cellranger/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0/",
        )

    def test_absolute_path_only_gains_the_trailing_separator(self):
        self.assertEqual(
            preprocess.reference_dir("s3://other/refdata-gex-mm10-2020-A", REFERENCES),
            "s3://other/refdata-gex-mm10-2020-A/",
        )

    def test_missing_value_is_left_alone(self):
        for value in (None, ""):
            self.assertEqual(preprocess.reference_dir(value, REFERENCES), value)


class ReadSamplesheet(unittest.TestCase):

    def test_values_stay_text(self):
        samplesheet = preprocess.read_samplesheet(io.StringIO(HEADER + GEX + VDJ))
        self.assertEqual(samplesheet.columns.tolist(), ["sample", "grouping", "feature_types"])
        self.assertEqual(samplesheet["grouping"].tolist(), ["01", "01"])

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


class ReadOptionalTable(unittest.TestCase):

    def test_untouched_or_empty_table_reads_as_empty(self):
        for path in (None, "", io.StringIO(""), io.StringIO("sample_id,probe_barcode_ids,description\r\n")):
            self.assertTrue(preprocess.read_optional_table(path).empty)

    def test_rows_are_read_as_text_with_blank_descriptions(self):
        probe_barcodes = preprocess.read_optional_table(io.StringIO(
            "sample_id,probe_barcode_ids,description\r\n01,BC001+AB005,\r\n02,BC002,Control\r\n"
        ))
        self.assertEqual(
            probe_barcodes.values.tolist(),
            [["01", "BC001+AB005", ""], ["02", "BC002", "Control"]],
        )


class DropRepeatedRows(unittest.TestCase):

    def test_rows_repeated_per_lane_become_one_library(self):
        samplesheet = preprocess.read_samplesheet(io.StringIO(HEADER + GEX + VDJ + GEX + VDJ))
        self.assertEqual(
            preprocess.drop_repeated_rows(samplesheet).values.tolist(),
            [["01_GEX", "01", "Gene Expression"], ["01_VDJ", "01", "VDJ-B"]],
        )


class UnsetUmiThresholds(unittest.TestCase):

    def test_blank_thresholds_are_unset(self):
        self.assertEqual(
            preprocess.unset_umi_thresholds({"cellranger_version": "6.1.2"}),
            ["emptydrops_minimum_umis", "global_minimum_umis"],
        )

    def test_threshold_newer_than_the_version_fails_clearly(self):
        params = {"cellranger_version": "7.2.0", "global_minimum_umis": 50}
        with self.assertRaisesRegex(ValueError, r"global_minimum_umis requires CellRanger 8\.0 or later"):
            preprocess.unset_umi_thresholds(params)


if __name__ == "__main__":
    unittest.main()
