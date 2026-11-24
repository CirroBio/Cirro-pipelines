"""Offline tests for the Funcotator preprocess script, shared by Funcotator 1.0 and
its DRAGEN germline and somatic variants.

preprocess.py imports the Cirro SDK only when it runs, so this is a plain
`python -m unittest` that needs pandas and nothing else.

The samplesheet is the structured input Cirro writes from the run form: a CSV with
CRLF line endings and a header from the form's columns.
"""
import importlib.util
import io
import unittest
from pathlib import Path


def _load_preprocess():
    path = Path(__file__).with_name("preprocess.py")
    spec = importlib.util.spec_from_file_location("funcotator_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()

HEADER = "sample,vcf\r\n"
REFERENCES = "s3://refs-bucket"


class AbsoluteReference(unittest.TestCase):

    def test_relative_path_is_put_under_the_references_bucket(self):
        self.assertEqual(
            preprocess.absolute_reference("funcotator/funcotator_dataSources.v1.7.20200521g/", REFERENCES),
            "s3://refs-bucket/funcotator/funcotator_dataSources.v1.7.20200521g/",
        )

    def test_absolute_path_is_left_alone(self):
        value = "s3://other-bucket/funcotator/funcotator_dataSources.v1.7.20200521s/"
        self.assertEqual(preprocess.absolute_reference(value, REFERENCES), value)

    def test_missing_value_is_left_alone(self):
        for value in (None, ""):
            self.assertEqual(preprocess.absolute_reference(value, REFERENCES), value)


class ReadSamplesheet(unittest.TestCase):

    def test_values_stay_text(self):
        samplesheet = preprocess.read_samplesheet(io.StringIO(
            HEADER + "01,s3://bucket/data/01.hard-filtered.vcf.gz\r\n"
        ))
        self.assertEqual(samplesheet["sample"].tolist(), ["01"])
        self.assertEqual(samplesheet["vcf"].tolist(), ["s3://bucket/data/01.hard-filtered.vcf.gz"])

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


if __name__ == "__main__":
    unittest.main()
