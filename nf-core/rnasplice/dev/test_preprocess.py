"""Offline tests for the nf-core/rnasplice samplesheet and contrasts handling.

preprocess.py imports the Cirro SDK only when it runs, so this is a plain
`python -m unittest` that needs pandas and nothing else.

The samplesheet is the structured input Cirro writes from the run form: a CSV with
CRLF line endings, a header from the form's columns, and blank cells left empty.
"""
import importlib.util
import io
import unittest
from pathlib import Path

import pandas as pd


def _load_preprocess():
    path = Path(__file__).with_name("preprocess.py")
    spec = importlib.util.spec_from_file_location("rnasplice_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()

HEADER = "sample,fastq_1,fastq_2,strandedness,condition\r\n"


def _samplesheet(samples, conditions):
    return pd.DataFrame({
        "sample": samples,
        "fastq_1": [f"{sample}_R1.fastq.gz" for sample in samples],
        "fastq_2": [f"{sample}_R2.fastq.gz" for sample in samples],
        "strandedness": ["auto"] * len(samples),
        "condition": conditions,
    })


class ReadSamplesheet(unittest.TestCase):

    def test_values_stay_text(self):
        samplesheet = preprocess.read_samplesheet(io.StringIO(
            HEADER + "S01,s3://bucket/data/S01_R1.fastq.gz,,,NA\r\n"
        ))
        self.assertEqual(samplesheet["fastq_2"].tolist(), [""])
        self.assertEqual(samplesheet["strandedness"].tolist(), [""])
        self.assertEqual(samplesheet["condition"].tolist(), ["NA"])

    def test_no_param_fails_clearly(self):
        for path in (None, ""):
            with self.assertRaisesRegex(ValueError, "No samplesheet was provided"):
                preprocess.read_samplesheet(path)

    def test_empty_file_fails_clearly(self):
        with self.assertRaisesRegex(ValueError, "The samplesheet is empty"):
            preprocess.read_samplesheet(io.StringIO(""))

    def test_header_without_rows_fails_clearly(self):
        with self.assertRaisesRegex(ValueError, "The samplesheet is empty"):
            preprocess.read_samplesheet(io.StringIO(HEADER))


class FillStrandedness(unittest.TestCase):

    def test_blank_becomes_auto_and_set_values_are_kept(self):
        samplesheet = _samplesheet(["A", "B"], ["ctrl", "trt"]).assign(strandedness=["", "reverse"])
        self.assertEqual(preprocess.fill_strandedness(samplesheet)["strandedness"].tolist(), ["auto", "reverse"])


class CheckConditions(unittest.TestCase):

    def test_two_conditions_with_repeated_samples_pass(self):
        preprocess.check_conditions(_samplesheet(["A", "A", "B"], ["ctrl", "ctrl", "trt"]))

    def test_missing_condition_names_the_samples(self):
        with self.assertRaisesRegex(ValueError, "missing for: B, C$"):
            preprocess.check_conditions(_samplesheet(["A", "B", "C", "C"], ["ctrl", "", "", ""]))

    def test_sample_with_two_conditions_is_named(self):
        with self.assertRaisesRegex(ValueError, "have several: A$"):
            preprocess.check_conditions(_samplesheet(["A", "A", "B"], ["ctrl", "trt", "trt"]))

    def test_one_condition_fails_clearly(self):
        with self.assertRaisesRegex(ValueError, "at least two conditions"):
            preprocess.check_conditions(_samplesheet(["A", "B"], ["ctrl", "ctrl"]))


class MakeContrasts(unittest.TestCase):

    def test_every_pair_of_conditions_once(self):
        contrasts = preprocess.make_contrasts(_samplesheet(["A", "B", "C", "D"], ["ctrl", "trt", "ctrl", "alt"]))
        self.assertEqual(contrasts.columns.tolist(), ["contrast", "treatment", "control"])
        self.assertEqual(
            contrasts.values.tolist(),
            [["ctrl_alt", "ctrl", "alt"], ["trt_ctrl", "trt", "ctrl"], ["trt_alt", "trt", "alt"]],
        )


if __name__ == "__main__":
    unittest.main()
