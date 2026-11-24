"""Offline tests for the gig-map align_genomes preprocess script.

preprocess.py imports the Cirro SDK only when it runs, so this is a plain
`python -m unittest` that needs nothing else.
"""
import importlib.util
import unittest
from pathlib import Path


def _load_preprocess():
    path = Path(__file__).with_name("preprocess.py")
    spec = importlib.util.spec_from_file_location("gig_map_align_genomes_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()


class CatalogGenes(unittest.TestCase):

    def test_folder_pickers_trailing_slash_is_not_doubled(self):
        self.assertEqual(
            preprocess.catalog_genes("s3://bucket/datasets/catalog/data/"),
            "s3://bucket/datasets/catalog/data/centroids.faa.gz",
        )

    def test_folder_without_trailing_slash(self):
        self.assertEqual(
            preprocess.catalog_genes("s3://bucket/datasets/catalog/data/sub"),
            "s3://bucket/datasets/catalog/data/sub/centroids.faa.gz",
        )


if __name__ == "__main__":
    unittest.main()
