"""Offline tests for the differential abundance (curatedMetagenomicData) preprocess script.

preprocess.py imports the Cirro SDK only when it runs, so this is a plain
`python -m unittest` with no dependencies.
"""
import importlib.util
import unittest
from pathlib import Path


def _load_preprocess():
    path = Path(__file__).with_name("preprocess.py")
    spec = importlib.util.spec_from_file_location("differential_abundance_cmd_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()

VALID = {"formula": "group", "fdr_cutoff": 0.05}


class CheckParams(unittest.TestCase):

    def test_valid_params_pass(self):
        preprocess.check_params(VALID)

    def test_missing_formula_fails(self):
        for formula in (None, ""):
            with self.assertRaisesRegex(ValueError, "Must provide formula"):
                preprocess.check_params({**VALID, "formula": formula})

    def test_fdr_cutoff_out_of_range_fails(self):
        for cutoff in (-0.1, 1, 1.5):
            with self.assertRaisesRegex(ValueError, "fdr_cutoff"):
                preprocess.check_params({**VALID, "fdr_cutoff": cutoff})


if __name__ == "__main__":
    unittest.main()
