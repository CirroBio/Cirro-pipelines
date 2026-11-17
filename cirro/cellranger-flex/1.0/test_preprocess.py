"""Offline tests for the CellRanger Flex preprocess script.

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
    spec = importlib.util.spec_from_file_location("cellranger_flex_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()

BARCODES_HEADER = "sample_id,probe_barcode_ids\r\n"


def read(text: str):
    return preprocess.read_table(io.StringIO(text), "table", "rows")


class ReadTable(unittest.TestCase):

    def test_values_stay_text(self):
        samplesheet = read("sample\r\n01\r\n02\r\n")
        self.assertEqual(samplesheet["sample"].tolist(), ["01", "02"])

    def test_no_param_fails_clearly(self):
        for path in (None, ""):
            with self.assertRaisesRegex(ValueError, "No probe barcode table was provided: add one row per sample"):
                preprocess.read_table(path, "probe barcode table", "one row per sample")

    def test_empty_file_fails_clearly(self):
        with self.assertRaisesRegex(ValueError, "The table is empty"):
            read("")

    def test_header_without_rows_fails_clearly(self):
        with self.assertRaisesRegex(ValueError, "The table is empty"):
            read("sample\r\n")


class DropRepeatedRows(unittest.TestCase):

    def test_rows_repeated_per_lane_become_one_library(self):
        samplesheet = read("sample\r\npool1\r\npool1\r\npool2\r\npool1\r\n")
        self.assertEqual(preprocess.drop_repeated_rows(samplesheet)["sample"].tolist(), ["pool1", "pool2"])


class UnsetUmiThresholds(unittest.TestCase):

    def test_blank_thresholds_are_unset(self):
        self.assertEqual(
            preprocess.unset_umi_thresholds({"cellranger_version": "10.0.0", "global_minimum_umis": 10}),
            ["emptydrops_minimum_umis"],
        )

    def test_threshold_newer_than_the_version_fails_clearly(self):
        params = {"cellranger_version": "7.0.1", "emptydrops_minimum_umis": 100}
        with self.assertRaisesRegex(ValueError, r"emptydrops_minimum_umis requires CellRanger 7\.1 or later"):
            preprocess.unset_umi_thresholds(params)


class ManualProbeBarcodes(unittest.TestCase):

    def test_each_named_barcode_names_its_sample(self):
        params = {"bc002": "tumor", "bc001": "normal", "bc003": None, "bc004": "", "barcode_source": "manual"}
        self.assertEqual(
            preprocess.manual_probe_barcodes(params).values.tolist(),
            [["normal", "BC001"], ["tumor", "BC002"]],
        )

    def test_no_named_barcode_gives_an_empty_table(self):
        probe_barcodes = preprocess.manual_probe_barcodes({"bc001": None})
        self.assertEqual(probe_barcodes.columns.tolist(), ["sample_id", "probe_barcode_ids"])
        self.assertTrue(probe_barcodes.empty)


class SelectProbeBarcodes(unittest.TestCase):

    def test_form_table_is_kept_as_is(self):
        table = read(BARCODES_HEADER + "01,BC001\r\n02,BC002|BC003\r\n")
        self.assertEqual(
            preprocess.select_probe_barcodes(table).values.tolist(),
            [["01", "BC001"], ["02", "BC002|BC003"]],
        )

    def test_rows_leaving_a_value_blank_are_skipped(self):
        table = read(BARCODES_HEADER + "01,BC001\r\n02,\r\n,BC003\r\n")
        self.assertEqual(preprocess.select_probe_barcodes(table).values.tolist(), [["01", "BC001"]])

    def test_saved_table_may_name_the_column_sample_barcode_ids(self):
        table = read("sample_id,sample_barcode_ids,description\r\n01,A-A01,first\r\n")
        self.assertEqual(preprocess.select_probe_barcodes(table).values.tolist(), [["01", "A-A01"]])

    def test_missing_column_fails_clearly(self):
        with self.assertRaisesRegex(ValueError, "must have a column named probe_barcode_ids"):
            preprocess.select_probe_barcodes(read("sample_id,barcode\r\n01,BC001\r\n"))


class JoinProbeBarcodes(unittest.TestCase):

    def test_barcodes_of_a_sample_on_several_rows_are_piped(self):
        table = read(BARCODES_HEADER + "02,BC003\r\n01,BC001\r\n01,BC002\r\n")
        self.assertEqual(
            preprocess.join_probe_barcodes(table).to_csv(index=None),
            "sample_id,probe_barcode_ids\n01,BC001|BC002\n02,BC003\n",
        )

    def test_no_barcodes_fails_clearly(self):
        with self.assertRaisesRegex(ValueError, "No sample was assigned a probe barcode"):
            preprocess.join_probe_barcodes(preprocess.manual_probe_barcodes({}))


if __name__ == "__main__":
    unittest.main()
