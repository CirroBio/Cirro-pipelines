"""Offline tests for the QuPath preprocess script.

preprocess.py imports the Cirro SDK only when it runs, so this is a plain
`python -m unittest` that needs nothing else.
"""
import importlib.util
import unittest
from pathlib import Path


def _load_preprocess():
    path = Path(__file__).with_name("preprocess.py")
    spec = importlib.util.spec_from_file_location("qupath_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()


class EmptyParams(unittest.TestCase):

    def test_blank_args_are_dropped(self):
        params = {
            "input_folder": "s3://bucket/datasets/images/data/",
            "script": "s3://refs/qupath_script/count/script.groovy",
            "args": "",
        }
        self.assertEqual(preprocess.empty_params(params), ["args"])

    def test_whitespace_args_are_dropped(self):
        self.assertEqual(preprocess.empty_params({"args": "  "}), ["args"])

    def test_set_args_are_kept(self):
        params = {
            "input_folder": "s3://bucket/datasets/images/data/",
            "args": "threshold=0.5",
        }
        self.assertEqual(preprocess.empty_params(params), [])


if __name__ == "__main__":
    unittest.main()
