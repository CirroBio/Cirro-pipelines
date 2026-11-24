"""Offline tests for the gig-map download-genes preprocess script.

preprocess.py imports the Cirro SDK only when it runs, so this is a plain
`python -m unittest` that needs nothing else.
"""
import importlib.util
import unittest
from pathlib import Path


def _load_preprocess():
    path = Path(__file__).with_name("preprocess.py")
    spec = importlib.util.spec_from_file_location("gig_map_download_genes_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()


class GenomeTableGlobs(unittest.TestCase):

    def test_one_dataset(self):
        self.assertEqual(
            preprocess.genome_table_globs("s3://bucket/datasets/a/data/"),
            {
                "genome_csv": "s3://bucket/datasets/a/data/*.csv",
                "genome_tsv": "s3://bucket/datasets/a/data/*.tsv",
            },
        )

    def test_several_folders_are_comma_separated(self):
        self.assertEqual(
            preprocess.genome_table_globs("s3://bucket/datasets/a/data/,s3://bucket/datasets/b/data/sub"),
            {
                "genome_csv": "s3://bucket/datasets/a/data/*.csv,s3://bucket/datasets/b/data/sub/*.csv",
                "genome_tsv": "s3://bucket/datasets/a/data/*.tsv,s3://bucket/datasets/b/data/sub/*.tsv",
            },
        )


if __name__ == "__main__":
    unittest.main()
