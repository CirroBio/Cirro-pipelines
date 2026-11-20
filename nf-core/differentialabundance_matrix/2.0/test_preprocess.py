"""Offline tests for the samplesheet, contrast and matrix built from a combined dataset.

preprocess.py imports the Cirro SDK only when it runs, so this is a plain
`python -m unittest` that needs pandas and nothing else.
"""
import importlib.util
import json
import unittest
from pathlib import Path

import pandas as pd


def _load_preprocess():
    path = Path(__file__).with_name("preprocess.py")
    spec = importlib.util.spec_from_file_location("differentialabundance_matrix_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()

SAMPLES = pd.DataFrame({
    "sample_id": ["P1 W0", "P1 W4", "P2 W0", "P2 W4", "P3 W0"],
    "subject_id": ["P1", "P1", "P2", "P2", "P3"],
    "diagnosis": ["CD", "CD", "nonIBD", "nonIBD", "CD"],
    "age": ["30", "30", "52", "", "41"],
})


class ParseColumns(unittest.TestCase):

    def test_unset_and_blank_give_no_columns(self):
        for value in (None, "", " , "):
            self.assertEqual(preprocess.parse_columns(value), [])

    def test_names_are_trimmed(self):
        self.assertEqual(preprocess.parse_columns("age, batch"), ["age", "batch"])


class CheckModelColumns(unittest.TestCase):

    def test_missing_column_fails_clearly(self):
        with self.assertRaisesRegex(ValueError, "samples.csv has no column 'sex'"):
            preprocess.check_model_columns(SAMPLES, ["diagnosis", "sex"])

    def test_name_r_cannot_use_fails_clearly(self):
        samples = SAMPLES.assign(**{"week 2": "1"})
        with self.assertRaisesRegex(ValueError, "can't be used in a model"):
            preprocess.check_model_columns(samples, ["week 2"])


class MakeSamplesheet(unittest.TestCase):

    def test_keeps_samples_with_the_matrix_and_every_value(self):
        samplesheet = preprocess.make_samplesheet(
            SAMPLES, ["P1 W0", "P2 W0", "P2 W4", "P3 W0"], ["diagnosis", "age"], "diagnosis", "nonIBD", "CD"
        )
        # P1 W4 lacks the matrix, and P2 W4 has no age
        self.assertEqual(samplesheet["sample"].tolist(), ["P1 W0", "P2 W0", "P3 W0"])
        self.assertEqual(samplesheet.columns[0], "sample")

    def test_level_left_without_samples_fails_clearly(self):
        with self.assertRaisesRegex(ValueError, "has diagnosis = 'nonIBD'"):
            preprocess.make_samplesheet(
                SAMPLES, ["P1 W0", "P2 W4", "P3 W0"], ["diagnosis", "age"], "diagnosis", "nonIBD", "CD"
            )

    def test_variable_named_sample_fails_clearly(self):
        with self.assertRaisesRegex(ValueError, "reserves for sample IDs"):
            preprocess.make_samplesheet(
                SAMPLES.assign(sample="x"), list(SAMPLES["sample_id"]), ["diagnosis"], "diagnosis", "nonIBD", "CD"
            )


class RName(unittest.TestCase):

    def test_matches_r_make_names(self):
        self.assertEqual(preprocess.r_name("factor(diagnosis)non IBD"), "factor.diagnosis.non.IBD")
        self.assertEqual(preprocess.r_name("1st"), "X1st")
        self.assertEqual(preprocess.r_name(".5"), "X.5")
        self.assertEqual(preprocess.r_name(".x"), ".x")


class MakeContrast(unittest.TestCase):

    TYPES = {"diagnosis": "categorical", "age": "numeric", "batch": "categorical"}

    def test_categorical_covariates_block_a_comparison(self):
        contrast = preprocess.make_contrast("diagnosis", "nonIBD", "CD", ["batch"], None, self.TYPES)
        self.assertEqual(contrast, {
            "id": "diagnosis_CD_vs_nonIBD",
            "comparison": ["diagnosis", "nonIBD", "CD"],
            "blocking_factors": ["batch"],
        })

    def test_no_covariates_give_a_bare_comparison(self):
        contrast = preprocess.make_contrast("diagnosis", "nonIBD", "CD", [], None, self.TYPES)
        self.assertNotIn("blocking_factors", contrast)

    def test_numeric_covariate_needs_a_formula(self):
        contrast = preprocess.make_contrast("diagnosis", "non IBD", "CD", ["age", "batch"], None, self.TYPES)
        self.assertEqual(contrast, {
            "id": "diagnosis_CD_vs_non_IBD",
            "formula": "~ 0 + factor(diagnosis) + age + factor(batch)",
            "make_contrasts_str": "factor.diagnosis.CD - factor.diagnosis.non.IBD",
        })

    def test_random_effect_needs_a_formula(self):
        contrast = preprocess.make_contrast("week", "0", "4", [], "subject_id", {"week": "numeric"})
        self.assertEqual(contrast["formula"], "~ 0 + factor(week) + (1 | subject_id)")
        self.assertEqual(contrast["make_contrasts_str"], "factor.week.4 - factor.week.0")

    def test_written_as_yaml_the_pipeline_reads(self):
        contrast = preprocess.make_contrast("diagnosis", "nonIBD", "CD", [], "subject_id", self.TYPES)
        self.assertEqual(json.loads(preprocess.contrasts_yaml(contrast)), {"contrasts": [contrast]})


class PrepareMatrix(unittest.TestCase):

    def test_renames_the_id_column_and_samples(self):
        matrix = pd.DataFrame({"feature_id": ["g1"], "A 1": [1], "B": [2], "C": [3]})
        prepared = preprocess.prepare_matrix(matrix, {"C": "C", "A 1": "A_1"})
        self.assertEqual(prepared.columns.tolist(), ["gene_id", "C", "A_1"])
        self.assertEqual(prepared.iloc[0].tolist(), ["g1", 3, 1])

    def test_other_id_column_fails_clearly(self):
        with self.assertRaisesRegex(ValueError, "first column must be feature_id"):
            preprocess.prepare_matrix(pd.DataFrame({"id": ["g1"], "A": [1]}), {"A": "A"})


class PipelineSampleNames(unittest.TestCase):

    def test_whitespace_becomes_underscores(self):
        self.assertEqual(
            preprocess.pipeline_sample_names(["P1 W0", "P1\tW4", "01"]),
            {"P1 W0": "P1_W0", "P1\tW4": "P1_W4", "01": "01"},
        )

    def test_names_that_collide_fail_clearly(self):
        with self.assertRaisesRegex(ValueError, "same once their spaces are replaced"):
            preprocess.pipeline_sample_names(["P1 W0", "P1_W0"])


class MethodParams(unittest.TestCase):

    def test_every_combination_names_a_method(self):
        for (values, repeated), params in preprocess.METHOD_PARAMS.items():
            expected = "dream" if repeated else {"counts": "deseq2", "log_abundance": "limma"}[values]
            self.assertEqual(params["differential_method"], expected)

    def test_small_count_matrix_is_stabilised_with_rlog(self):
        self.assertEqual(preprocess.method_params("counts", False, 300)["deseq2_vs_method"], "rlog")
        self.assertNotIn("deseq2_vs_method", preprocess.method_params("counts", False, 20000))
        self.assertNotIn("deseq2_vs_method", preprocess.method_params("counts", True, 300))


if __name__ == "__main__":
    unittest.main()
