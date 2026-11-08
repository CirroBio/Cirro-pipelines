"""Offline tests for the Sage Bionetworks nf-artist preprocess script.

preprocess.py imports the Cirro SDK only when it runs, so this is a plain
`python -m unittest` that needs pandas and nothing else.

The two image pickers are multi-selects, so each param is a list of S3 URIs, or
absent when nothing was picked.
"""
import importlib.util
import io
import unittest
from pathlib import Path

import pandas as pd


def _load_preprocess():
    path = Path(__file__).with_name("preprocess.py")
    spec = importlib.util.spec_from_file_location("sage_artist_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()

CYCIF = "s3://bucket/datasets/ds/data/exemplar-001.ome.tif"
HE = "s3://bucket/datasets/ds/data/CMU-1-Small-Region.svs"


class MakeSamplesheet(unittest.TestCase):

    def test_both_kinds_of_image(self):
        samplesheet = preprocess.make_samplesheet([CYCIF], [HE])
        self.assertEqual(
            samplesheet.columns.tolist(),
            ["image", "convert", "he", "minerva", "miniature"],
        )
        self.assertEqual(samplesheet["image"].tolist(), [CYCIF, HE])
        self.assertEqual(samplesheet["he"].tolist(), [False, True])
        self.assertTrue(samplesheet[["convert", "minerva", "miniature"]].all(axis=None))

    def test_one_kind_of_image(self):
        for images, he_images, he in (([CYCIF], None, False), (None, [HE], True)):
            with self.subTest(he=he):
                samplesheet = preprocess.make_samplesheet(images, he_images)
                self.assertEqual(samplesheet["he"].tolist(), [he])

    def test_blank_entries_are_dropped(self):
        samplesheet = preprocess.make_samplesheet([CYCIF, ""], [None, HE])
        self.assertEqual(samplesheet["image"].tolist(), [CYCIF, HE])

    def test_no_images_fails_clearly(self):
        for images, he_images in ((None, None), ([], []), ([""], [None])):
            with self.subTest(images=images, he_images=he_images):
                with self.assertRaisesRegex(ValueError, "No images selected"):
                    preprocess.make_samplesheet(images, he_images)

    def test_booleans_read_back_as_nf_artist_expects(self):
        # nf-artist parses each flag with Groovy's String.toBoolean(), which accepts "True"
        written = preprocess.make_samplesheet([CYCIF], [HE]).to_csv(index=False)
        read_back = pd.read_csv(io.StringIO(written), dtype=str)
        self.assertEqual(read_back["he"].tolist(), ["False", "True"])
        self.assertEqual(read_back["convert"].tolist(), ["True", "True"])


if __name__ == "__main__":
    unittest.main()
