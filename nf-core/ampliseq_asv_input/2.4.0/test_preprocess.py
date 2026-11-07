"""Offline tests for the nf-core/ampliseq ASV-input preprocess script.

preprocess.py imports the Cirro SDK only when it runs, so this is a plain
`python -m unittest` that needs nothing beyond the standard library.
"""
import importlib.util
import unittest
from pathlib import Path


def _load_preprocess():
    path = Path(__file__).with_name("preprocess.py")
    spec = importlib.util.spec_from_file_location("ampliseq_asv_input_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()

REFERENCES = "s3://refs-bucket"


class StagedReference(unittest.TestCase):

    def test_staged_taxonomy_is_swapped_in(self):
        self.assertEqual(
            preprocess.staged_reference({"dada_ref_taxonomy": "silva=138"}, REFERENCES),
            {
                "dada_ref_tax_custom":
                    "s3://refs-bucket/ampliseq/silva-138.1/silva_nr99_v138.1_wSpecies_train_set.fa.gz",
                "dada_ref_tax_custom_sp":
                    "s3://refs-bucket/ampliseq/silva-138.1/silva_species_assignment_v138.1.fa.gz",
            },
        )

    def test_taxonomy_without_a_staged_copy_still_downloads(self):
        for taxonomy in ("silva", "gtdb=R07-RS207", None):
            self.assertIsNone(
                preprocess.staged_reference({"dada_ref_taxonomy": taxonomy}, REFERENCES)
            )

    def test_custom_reference_set_by_the_user_wins(self):
        params = {
            "dada_ref_taxonomy": "silva=138",
            "dada_ref_tax_custom": "s3://bucket/my_taxonomy.fa.gz",
        }
        self.assertIsNone(preprocess.staged_reference(params, REFERENCES))


if __name__ == "__main__":
    unittest.main()
