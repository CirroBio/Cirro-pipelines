"""Offline tests for the Cell Ranger V(D)J preprocess script.

preprocess.py imports the Cirro SDK only when it runs, so this is a plain
`python -m unittest` with no dependencies.
"""
import importlib.util
import unittest
from pathlib import Path


def _load_preprocess():
    path = Path(__file__).with_name("preprocess.py")
    spec = importlib.util.spec_from_file_location("cellranger_vdj_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()

REFERENCES = "s3://refs-bucket"


class ReferenceDir(unittest.TestCase):

    def test_relative_path_is_put_under_the_references_bucket_as_a_folder(self):
        self.assertEqual(
            preprocess.reference_dir("cellranger/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0", REFERENCES),
            "s3://refs-bucket/cellranger/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0/",
        )

    def test_absolute_path_is_kept_with_one_trailing_slash(self):
        for value in ("s3://other/vdj-ref", "s3://other/vdj-ref/"):
            self.assertEqual(preprocess.reference_dir(value, REFERENCES), "s3://other/vdj-ref/")


class FastqFolder(unittest.TestCase):

    def test_whole_dataset_pick_loses_its_trailing_slash(self):
        self.assertEqual(
            preprocess.fastq_folder("s3://bucket/datasets/fastq-1/data/"),
            "s3://bucket/datasets/fastq-1/data",
        )

    def test_subfolder_is_kept_as_is(self):
        self.assertEqual(
            preprocess.fastq_folder("s3://bucket/datasets/fastq-1/data/run1"),
            "s3://bucket/datasets/fastq-1/data/run1",
        )

    def test_no_folder_fails_clearly(self):
        for value in (None, ""):
            with self.assertRaisesRegex(ValueError, "No FASTQ folder was selected"):
                preprocess.fastq_folder(value)


if __name__ == "__main__":
    unittest.main()
