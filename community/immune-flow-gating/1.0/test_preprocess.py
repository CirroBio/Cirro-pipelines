"""Offline tests for the immune flow gating preprocess script.

preprocess.py imports the Cirro SDK only when it runs, so this is a plain
`python -m unittest` that needs pandas and nothing else.

The samplesheet is the structured input Cirro writes from the run form: a CSV with
CRLF line endings, a header from the form's columns, and one row per FCS file.
"""
import importlib.util
import io
import unittest
from pathlib import Path


def _load_preprocess():
    path = Path(__file__).with_name("preprocess.py")
    spec = importlib.util.spec_from_file_location("immune_flow_gating_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()

HEADER = "sample,fcs,batch\r\n"


def read(rows: str):
    return preprocess.read_samplesheet(io.StringIO(HEADER + rows))


class BuildChannelMap(unittest.TestCase):

    def test_default_channels_need_no_alias(self):
        self.assertEqual(preprocess.build_channel_map(dict(preprocess.MARKER_CHANNELS)), [])

    def test_changed_channel_is_aliased(self):
        params = dict(preprocess.MARKER_CHANNELS, CD4="APC-eFluor780-A")
        self.assertEqual(
            preprocess.build_channel_map(params),
            [dict(channels="APC-eFluor780-A", alias="APC-eFluor 780-A")],
        )


class ReadSamplesheet(unittest.TestCase):

    def test_values_stay_text(self):
        samplesheet = read("01,s3://bucket/data/plate1/01.fcs,\r\n")
        self.assertEqual(samplesheet["sample"].tolist(), ["01"])
        self.assertEqual(samplesheet["batch"].tolist(), [""])

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


class BuildWorkflowSamplesheet(unittest.TestCase):

    def test_files_are_named_without_their_folder(self):
        samplesheet = preprocess.build_workflow_samplesheet(read(
            "01,s3://bucket/data/plate1/01.fcs,\r\n02,s3://bucket/data/plate2/02.FCS,\r\n"
        ))
        self.assertEqual(samplesheet.columns.tolist(), ["sample", "file_1"])
        self.assertEqual(samplesheet["file_1"].tolist(), ["01.fcs", "02.FCS"])

    def test_filled_batch_is_kept(self):
        samplesheet = preprocess.build_workflow_samplesheet(read(
            "01,s3://bucket/data/01.fcs,A\r\n02,s3://bucket/data/02.fcs,B\r\n"
        ))
        self.assertEqual(samplesheet["batch"].tolist(), ["A", "B"])

    def test_repeated_file_name_fails_naming_it(self):
        with self.assertRaisesRegex(ValueError, "listed only once.*: 01.fcs$"):
            preprocess.build_workflow_samplesheet(read(
                "01,s3://bucket/data/plate1/01.fcs,\r\n01b,s3://bucket/data/plate2/01.fcs,\r\n"
            ))


if __name__ == "__main__":
    unittest.main()
