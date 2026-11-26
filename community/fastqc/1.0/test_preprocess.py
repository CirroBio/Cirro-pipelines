"""Offline tests for the FASTQC preprocess script.

preprocess.py imports the Cirro SDK only when it runs, so this is a plain
`python -m unittest` that needs nothing beyond the standard library.
"""
import importlib.util
import unittest
from pathlib import Path


def _load_preprocess():
    path = Path(__file__).with_name("preprocess.py")
    spec = importlib.util.spec_from_file_location("fastqc_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()


class InputFolders(unittest.TestCase):

    def test_every_folder_ends_in_one_slash(self):
        self.assertEqual(
            preprocess.input_folders("s3://bucket/a/data/,s3://bucket/b/data/reads"),
            "s3://bucket/a/data/,s3://bucket/b/data/reads/",
        )

    def test_no_folder_fails_clearly(self):
        for value in (None, ""):
            with self.assertRaisesRegex(ValueError, "No input folders were provided"):
                preprocess.input_folders(value)


if __name__ == "__main__":
    unittest.main()
