"""Offline tests for the Space Ranger mkref preprocess script.

preprocess.py imports the Cirro SDK only when it runs, so this is a plain
`python -m unittest` with no dependencies.
"""
import importlib.util
import unittest
from pathlib import Path


def _load_preprocess():
    path = Path(__file__).with_name("preprocess.py")
    spec = importlib.util.spec_from_file_location("spaceranger_mkref_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()

REFERENCES = "s3://refs-bucket"


class AbsoluteReference(unittest.TestCase):

    def test_relative_path_is_put_under_the_references_bucket(self):
        self.assertEqual(
            preprocess.absolute_reference("cellranger/refdata-gex-GRCh38-2024-A", REFERENCES),
            "s3://refs-bucket/cellranger/refdata-gex-GRCh38-2024-A",
        )

    def test_absolute_path_is_left_alone(self):
        value = "s3://other-bucket/cellranger/flex/probes.csv"
        self.assertEqual(preprocess.absolute_reference(value, REFERENCES), value)


class GenomeFiles(unittest.TestCase):

    def test_fasta_and_gtf_are_inside_the_reference(self):
        for genome_dir in ("s3://refs-bucket/ref", "s3://refs-bucket/ref/"):
            self.assertEqual(preprocess.genome_files(genome_dir), {
                "fasta": "s3://refs-bucket/ref/fasta/genome.fa",
                "genes": "s3://refs-bucket/ref/genes/genes.gtf",
            })


class CustomGlobs(unittest.TestCase):

    def test_whole_dataset_pick_gives_single_slashes(self):
        self.assertEqual(
            preprocess.custom_globs("s3://bucket/datasets/fasta-1/data/"),
            "s3://bucket/datasets/fasta-1/data/**",
        )

    def test_several_folders_stay_comma_separated(self):
        self.assertEqual(
            preprocess.custom_globs("s3://bucket/datasets/fasta-1/data/,s3://bucket/datasets/fasta-2/data/gfp"),
            "s3://bucket/datasets/fasta-1/data/**,s3://bucket/datasets/fasta-2/data/gfp/**",
        )

    def test_no_folder_fails_clearly(self):
        for value in (None, ""):
            with self.assertRaisesRegex(ValueError, "No folder of custom sequences was selected"):
                preprocess.custom_globs(value)


if __name__ == "__main__":
    unittest.main()
