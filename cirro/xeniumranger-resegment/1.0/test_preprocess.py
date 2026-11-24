"""Offline tests for the Xenium Ranger resegment preprocess script.

preprocess.py imports the Cirro SDK only when it runs, so this is a plain
`python -m unittest` that needs nothing else.
"""
import importlib.util
import unittest
from pathlib import Path


def _load_preprocess():
    path = Path(__file__).with_name("preprocess.py")
    spec = importlib.util.spec_from_file_location("xeniumranger_resegment_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()


class BundleFolder(unittest.TestCase):

    def test_picked_file_gives_its_bundle(self):
        self.assertEqual(
            preprocess.bundle_folder("s3://bucket/datasets/xen-1/data/output-XETG/analysis.zarr.zip"),
            "s3://bucket/datasets/xen-1/data/output-XETG",
        )

    def test_bundle_at_the_top_of_a_dataset(self):
        self.assertEqual(
            preprocess.bundle_folder("s3://bucket/datasets/xen-1/data/analysis.zarr.zip"),
            "s3://bucket/datasets/xen-1/data",
        )

    def test_no_input_fails_clearly(self):
        for value in (None, ""):
            with self.assertRaisesRegex(ValueError, "No input was provided"):
                preprocess.bundle_folder(value)

    def test_other_file_fails_clearly(self):
        with self.assertRaisesRegex(ValueError, "Expected the input to end with"):
            preprocess.bundle_folder("s3://bucket/datasets/xen-1/data/output-XETG/cells.zarr.zip")


if __name__ == "__main__":
    unittest.main()
