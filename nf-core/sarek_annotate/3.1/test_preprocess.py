"""Offline tests for the nf-core/sarek annotate samplesheet handling.

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
    spec = importlib.util.spec_from_file_location("sarek_annotate_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()

HEADER = "patient,sample,variantcaller,vcf\r\n"


class ReadSamplesheet(unittest.TestCase):

    def test_values_stay_text(self):
        samplesheet = preprocess.read_samplesheet(io.StringIO(
            HEADER + "01,001,,s3://bucket/data/001.vcf.gz\r\n"
        ))
        self.assertEqual(samplesheet["patient"].tolist(), ["01"])
        self.assertEqual(samplesheet["sample"].tolist(), ["001"])
        self.assertEqual(samplesheet["variantcaller"].tolist(), [""])

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


class CheckUniqueSamples(unittest.TestCase):

    def test_distinct_samples_pass(self):
        preprocess.check_unique_samples(pd.DataFrame({
            "patient": ["P1", "P1", "P2"],
            "sample": ["N", "T", "N"],
            "variantcaller": ["haplotypecaller", "haplotypecaller", "haplotypecaller"],
            "vcf": ["p1n.vcf.gz", "p1t.vcf.gz", "p2n.vcf.gz"],
        }))

    def test_same_sample_from_two_callers_fails(self):
        samplesheet = pd.DataFrame({
            "patient": ["P1", "P1", "P2"],
            "sample": ["N", "N", "N"],
            "variantcaller": ["strelka", "manta", "strelka"],
            "vcf": ["a.vcf.gz", "b.vcf.gz", "c.vcf.gz"],
        })
        with self.assertRaisesRegex(ValueError, r"on more than one row: P1/N\. Give"):
            preprocess.check_unique_samples(samplesheet)


if __name__ == "__main__":
    unittest.main()
