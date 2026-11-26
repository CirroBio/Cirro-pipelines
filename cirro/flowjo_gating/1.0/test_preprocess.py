"""Offline tests for the Apply FlowJo Gates preprocess script.

preprocess.py imports the Cirro SDK only when it runs, so this is a plain
`python -m unittest` that needs pandas and nothing else.

The FCS file table is the structured input Cirro writes from the run form: a CSV with
CRLF line endings, a header from the form's columns, and one row per FCS file.
"""
import importlib.util
import io
import logging
import unittest
from pathlib import Path


def _load_preprocess():
    path = Path(__file__).with_name("preprocess.py")
    spec = importlib.util.spec_from_file_location("flowjo_gating_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()

HEADER = "fcs\r\n"
LOGGER = logging.getLogger(__name__)


class ReadFcsFiles(unittest.TestCase):

    def test_rows_are_read(self):
        fcs_files = preprocess.read_fcs_files(io.StringIO(
            HEADER + "s3://bucket/data/01.fcs\r\ns3://bucket/data/02.FCS\r\n"
        ))
        self.assertEqual(fcs_files["fcs"].tolist(), ["s3://bucket/data/01.fcs", "s3://bucket/data/02.FCS"])

    def test_no_param_fails_clearly(self):
        for path in (None, ""):
            with self.assertRaisesRegex(ValueError, "No FCS files were provided"):
                preprocess.read_fcs_files(path)

    def test_empty_file_fails_clearly(self):
        with self.assertRaisesRegex(ValueError, "The FCS file table is empty"):
            preprocess.read_fcs_files(io.StringIO(""))

    def test_header_without_rows_fails_clearly(self):
        with self.assertRaisesRegex(ValueError, "The FCS file table is empty"):
            preprocess.read_fcs_files(io.StringIO(HEADER))


class DeduplicateFcs(unittest.TestCase):

    def test_distinct_names_are_joined(self):
        self.assertEqual(
            preprocess.deduplicate_fcs(["s3://b/data/a/01.fcs", "s3://b/data/a/02.FCS"], LOGGER),
            "s3://b/data/a/01.fcs,s3://b/data/a/02.FCS",
        )

    def test_repeated_name_in_another_folder_keeps_the_first(self):
        with self.assertLogs(LOGGER, level="WARNING"):
            self.assertEqual(
                preprocess.deduplicate_fcs(["s3://b/data/a/01.fcs", "s3://b/data/b/01.fcs"], LOGGER),
                "s3://b/data/a/01.fcs",
            )


if __name__ == "__main__":
    unittest.main()
