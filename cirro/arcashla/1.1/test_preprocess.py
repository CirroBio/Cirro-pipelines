"""Offline tests for the arcasHLA 1.1 preprocess script.

preprocess.py imports the Cirro SDK only when it runs, so this is a plain
`python -m unittest` that needs pandas and nothing else.

The samplesheet is the structured input Cirro writes from the run form: a CSV with
CRLF line endings, a header from the form's columns, and one row per sample.
"""
import importlib.util
import io
import unittest
from pathlib import Path


def _load_preprocess():
    path = Path(__file__).with_name("preprocess.py")
    spec = importlib.util.spec_from_file_location("arcashla_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()

HEADER = "sample,bam\r\n"
ROW_01 = "01,s3://bucket/data/star_salmon/01.markdup.sorted.bam\r\n"
ROW_02 = "02,s3://bucket/data/star_salmon/02.markdup.sorted.bam\r\n"


class ReadSamplesheet(unittest.TestCase):

    def test_values_stay_text(self):
        samplesheet = preprocess.read_samplesheet(io.StringIO(HEADER + ROW_01))
        self.assertEqual(samplesheet.columns.tolist(), ["sample", "bam"])
        self.assertEqual(samplesheet["sample"].tolist(), ["01"])

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
        preprocess.check_unique_samples(
            preprocess.read_samplesheet(io.StringIO(HEADER + ROW_01 + ROW_02))
        )

    def test_sample_from_two_datasets_fails_naming_it(self):
        samplesheet = preprocess.read_samplesheet(io.StringIO(
            HEADER + ROW_01 + "01,s3://bucket/other/star_salmon/01.markdup.sorted.bam\r\n" + ROW_02
        ))
        with self.assertRaisesRegex(ValueError, "Listed more than once: 01$"):
            preprocess.check_unique_samples(samplesheet)


class FormatGenes(unittest.TestCase):

    def test_all_wins_over_named_genes(self):
        self.assertEqual(preprocess.format_genes(["A", "all"]), "all")

    def test_named_genes_are_comma_joined(self):
        self.assertEqual(preprocess.format_genes(["A", "B", "DRB1"]), "A,B,DRB1")

    def test_no_genes_fails_clearly(self):
        with self.assertRaisesRegex(ValueError, "at least 1 gene"):
            preprocess.format_genes([])


if __name__ == "__main__":
    unittest.main()
