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


class DistinguishRepeatedSamples(unittest.TestCase):

    def test_distinct_samples_are_unchanged(self):
        samplesheet = pd.DataFrame({
            "patient": ["P1", "P1", "P2"],
            "sample": ["N", "T", "N"],
            "variantcaller": ["haplotypecaller", "haplotypecaller", "haplotypecaller"],
            "vcf": ["p1n.vcf.gz", "p1t.vcf.gz", "p2n.vcf.gz"],
        })
        self.assertTrue(preprocess.distinguish_repeated_samples(samplesheet).equals(samplesheet))

    def test_same_sample_from_two_callers_gets_the_caller_added(self):
        samplesheet = pd.DataFrame({
            "patient": ["P1", "P1", "P2"],
            "sample": ["N", "N", "N"],
            "variantcaller": ["strelka", "manta", "strelka"],
            "vcf": ["a.vcf.gz", "b.vcf.gz", "c.vcf.gz"],
        })
        distinct = preprocess.distinguish_repeated_samples(samplesheet)
        self.assertEqual(distinct["sample"].tolist(), ["N_strelka", "N_manta", "N"])
        self.assertEqual(distinct["vcf"].tolist(), samplesheet["vcf"].tolist())

    def test_same_caller_twice_is_named_after_the_vcf(self):
        samplesheet = pd.DataFrame({
            "patient": ["P1", "P1"],
            "sample": ["T_vs_N", "T_vs_N"],
            "variantcaller": ["strelka", "strelka"],
            "vcf": [
                "s3://bucket/strelka/T_vs_N/T_vs_N.strelka.somatic_snvs.vcf.gz",
                "s3://bucket/strelka/T_vs_N/T_vs_N.strelka.somatic_indels.vcf.gz",
            ],
        })
        distinct = preprocess.distinguish_repeated_samples(samplesheet)
        self.assertEqual(
            distinct["sample"].tolist(),
            ["T_vs_N.strelka.somatic_snvs", "T_vs_N.strelka.somatic_indels"],
        )

    def test_blank_caller_is_named_after_the_vcf(self):
        samplesheet = pd.DataFrame({
            "patient": ["P1", "P1"],
            "sample": ["N", "N"],
            "variantcaller": ["", ""],
            "vcf": ["s3://bucket/a/N.snv.vcf", "s3://bucket/b/N.sv.vcf.gz"],
        })
        distinct = preprocess.distinguish_repeated_samples(samplesheet)
        self.assertEqual(distinct["sample"].tolist(), ["N.snv", "N.sv"])

    def test_same_file_name_twice_fails(self):
        samplesheet = pd.DataFrame({
            "patient": ["P1", "P1"],
            "sample": ["N", "N"],
            "variantcaller": ["", ""],
            "vcf": ["s3://bucket/run1/N.vcf.gz", "s3://bucket/run2/N.vcf.gz"],
        })
        with self.assertRaisesRegex(ValueError, "even after adding the variant caller and VCF file name: P1/N"):
            preprocess.distinguish_repeated_samples(samplesheet)


if __name__ == "__main__":
    unittest.main()
