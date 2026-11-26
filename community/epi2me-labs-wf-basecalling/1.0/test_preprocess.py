"""Offline tests for the wf-basecalling preprocess script.

preprocess.py imports the Cirro SDK only when it runs, so this is a plain
`python -m unittest` that needs nothing else.
"""
import importlib.util
import unittest
from pathlib import Path


def _load_preprocess():
    path = Path(__file__).with_name("preprocess.py")
    spec = importlib.util.spec_from_file_location("wf_basecalling_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()


class EmptyParams(unittest.TestCase):

    def test_blank_strings_are_dropped(self):
        params = {
            "input": "s3://bucket/datasets/pod5/data/",
            "barcode_kit": "",
            "remora_cfg": "  ",
            "ref": "",
        }
        self.assertEqual(preprocess.empty_params(params), ["barcode_kit", "remora_cfg", "ref"])

    def test_set_values_are_kept(self):
        params = {
            "input": "s3://bucket/datasets/pod5/data/",
            "qscore_filter": 0,
            "duplex": False,
            "barcode_kit": "SQK-NBD114-24",
        }
        self.assertEqual(preprocess.empty_params(params), [])


if __name__ == "__main__":
    unittest.main()
