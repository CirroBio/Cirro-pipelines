"""Offline tests for the InterOp QC preprocess script.

preprocess.py imports the Cirro SDK only when it runs, so this is a plain
`python -m unittest` with no dependencies.
"""
import importlib.util
import unittest
from pathlib import Path


def _load_preprocess():
    path = Path(__file__).with_name("preprocess.py")
    spec = importlib.util.spec_from_file_location("interop_qc_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()


class RunFolder(unittest.TestCase):

    def test_whole_dataset_pick_loses_its_trailing_slash(self):
        self.assertEqual(
            preprocess.run_folder("s3://bucket/datasets/abc/data/"),
            "s3://bucket/datasets/abc/data",
        )

    def test_subfolder_pick_is_kept(self):
        self.assertEqual(
            preprocess.run_folder("s3://bucket/datasets/abc/data/run1"),
            "s3://bucket/datasets/abc/data/run1",
        )

    def test_missing_folder_fails_clearly(self):
        for value in (None, ""):
            with self.assertRaisesRegex(ValueError, "No sequencing run folder"):
                preprocess.run_folder(value)


if __name__ == "__main__":
    unittest.main()
