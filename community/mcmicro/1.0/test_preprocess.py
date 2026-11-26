"""Offline tests for the MCMICRO preprocess script.

preprocess.py imports the Cirro SDK and boto3 only when it runs, so this is a plain
`python -m unittest` that needs pandas and nothing else.

The markers table is the structured input Cirro writes from the run form: a CSV with
CRLF line endings, a header from the form's columns, and one row per image.
"""
import importlib.util
import io
import unittest
from pathlib import Path

import pandas as pd


def _load_preprocess():
    path = Path(__file__).with_name("preprocess.py")
    spec = importlib.util.spec_from_file_location("mcmicro_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()

HEADER = "image,cycle,marker_name\r\n"


class ReadMarkers(unittest.TestCase):

    def test_rows_are_read_as_text(self):
        markers = preprocess.read_markers(io.StringIO(
            HEADER + "s3://bucket/data/cyc001_reg001.tif,1,DAPI\r\n"
        ))
        self.assertEqual(markers.to_dict("records"), [
            {"image": "s3://bucket/data/cyc001_reg001.tif", "cycle": "1", "marker_name": "DAPI"},
        ])

    def test_no_param_fails_clearly(self):
        for path in (None, ""):
            with self.assertRaisesRegex(ValueError, "No markers table was provided"):
                preprocess.read_markers(path)

    def test_empty_file_fails_clearly(self):
        with self.assertRaisesRegex(ValueError, "The markers table is empty"):
            preprocess.read_markers(io.StringIO(""))

    def test_header_without_rows_fails_clearly(self):
        with self.assertRaisesRegex(ValueError, "The markers table is empty"):
            preprocess.read_markers(io.StringIO(HEADER))


class MarkerTable(unittest.TestCase):

    def test_rows_follow_file_name_order_and_keep_every_marker(self):
        markers = pd.DataFrame({
            "image": [
                "s3://bucket/b/data/round2.tif",
                "s3://bucket/a/data/round3.tif",
                "s3://bucket/a/data/round1.tif",
            ],
            "cycle": ["2", "3", "1"],
            "marker_name": ["CD3", "Blank", "DAPI"],
        })
        self.assertEqual(preprocess.marker_table(markers).to_dict("records"), [
            {"image": "s3://bucket/a/data/round1.tif", "cycle": 1, "marker_name": "DAPI"},
            {"image": "s3://bucket/b/data/round2.tif", "cycle": 2, "marker_name": "CD3"},
            {"image": "s3://bucket/a/data/round3.tif", "cycle": 3, "marker_name": "Blank"},
        ])

    def test_repeated_cycles_are_allowed(self):
        markers = pd.DataFrame({
            "image": ["s3://bucket/data/c1_a.tif", "s3://bucket/data/c1_b.tif"],
            "cycle": ["1", "1"],
            "marker_name": ["DAPI", "CD3"],
        })
        self.assertEqual(preprocess.marker_table(markers)["cycle"].tolist(), [1, 1])

    def test_cycles_out_of_file_name_order_fail_clearly(self):
        markers = pd.DataFrame({
            "image": ["s3://bucket/data/a.tif", "s3://bucket/data/b.tif"],
            "cycle": ["2", "1"],
            "marker_name": ["DAPI", "CD3"],
        })
        with self.assertRaisesRegex(ValueError, r"must not go down.*a\.tif \(2\), b\.tif \(1\)"):
            preprocess.marker_table(markers)

    def test_cycle_that_is_not_a_whole_number_fails_clearly(self):
        markers = pd.DataFrame({
            "image": ["s3://bucket/data/a.tif", "s3://bucket/data/b.tif"],
            "cycle": ["1", "two"],
            "marker_name": ["DAPI", "CD3"],
        })
        with self.assertRaisesRegex(ValueError, "whole numbers, not: two"):
            preprocess.marker_table(markers)

    def test_repeated_file_names_fail_clearly(self):
        markers = pd.DataFrame({
            "image": ["s3://bucket/a/data/cyc001_reg001.tif", "s3://bucket/b/data/cyc001_reg001.tif"],
            "cycle": ["1", "2"],
            "marker_name": ["DAPI", "CD3"],
        })
        with self.assertRaisesRegex(ValueError, "must be unique.*cyc001_reg001.tif"):
            preprocess.marker_table(markers)


if __name__ == "__main__":
    unittest.main()
