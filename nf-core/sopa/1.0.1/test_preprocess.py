"""Offline tests for the sopa samplesheet built from the images picked in the form.

The Cirro SDK and boto3 are stubbed, so this is a plain `python -m unittest` that needs
pandas and nothing else.

Regression coverage for: the PhenoCycler picker allows one file and sends its path as a
string, which was iterated as if it were the OME-TIFF picker's list, so every character
of the path became a row of the samplesheet.
"""
import importlib.util
import sys
import types
import unittest
from pathlib import Path


def _load_preprocess():
    helpers = types.ModuleType("cirro.helpers.preprocess_dataset")
    helpers.PreprocessDataset = object
    s3_path_mod = types.ModuleType("cirro.models.s3_path")
    s3_path_mod.S3Path = object
    for name in ["cirro", "cirro.helpers", "cirro.models"]:
        sys.modules.setdefault(name, types.ModuleType(name))
    sys.modules["cirro.helpers.preprocess_dataset"] = helpers
    sys.modules["cirro.models.s3_path"] = s3_path_mod
    sys.modules.setdefault("boto3", types.ModuleType("boto3"))

    path = Path(__file__).with_name("preprocess.py")
    spec = importlib.util.spec_from_file_location("sopa_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()

DATA_PATH = "s3://bucket/datasets/ds1/data"
PHENOCYCLER_IMAGE = f"{DATA_PATH}/slide1/scan.qptiff"


class PickedImageFiles(unittest.TestCase):

    def test_single_file_picker_sends_a_string(self):
        self.assertEqual(preprocess.picked_image_files(PHENOCYCLER_IMAGE), [PHENOCYCLER_IMAGE])

    def test_multiple_file_picker_sends_a_list(self):
        picked = [f"{DATA_PATH}/a.ome.tif", f" {DATA_PATH}/b.ome.tif ", ""]
        self.assertEqual(
            preprocess.picked_image_files(picked),
            [f"{DATA_PATH}/a.ome.tif", f"{DATA_PATH}/b.ome.tif"],
        )

    def test_nothing_picked(self):
        for picked in (None, "", [], ["  "]):
            self.assertEqual(preprocess.picked_image_files(picked), [])


class BuildSamplesheet(unittest.TestCase):

    def test_phenocycler_image_is_one_row(self):
        samplesheet = preprocess.build_samplesheet(
            DATA_PATH, "phenocycler", preprocess.picked_image_files(PHENOCYCLER_IMAGE)
        )
        self.assertEqual(samplesheet.to_dict("records"), [
            {"sample": "scan", "data_path": PHENOCYCLER_IMAGE},
        ])

    def test_ome_tif_images_are_one_row_each(self):
        picked = [f"{DATA_PATH}/a.ome.tif", f"{DATA_PATH}/b.ome.tif"]
        samplesheet = preprocess.build_samplesheet(
            DATA_PATH, "ome_tif", preprocess.picked_image_files(picked)
        )
        self.assertEqual(samplesheet["sample"].tolist(), ["a", "b"])
        self.assertEqual(samplesheet["data_path"].tolist(), picked)

    def test_no_image_picked_fails_clearly(self):
        with self.assertRaisesRegex(ValueError, "needs an image file"):
            preprocess.build_samplesheet(DATA_PATH, "phenocycler", preprocess.picked_image_files(None))


if __name__ == "__main__":
    unittest.main()
