"""Offline tests for the sourmash 1.0 preprocess script.

preprocess.py imports the Cirro SDK only when it runs, so this is a plain
`python -m unittest` that needs pandas and nothing else.

The samplesheet is the structured input Cirro writes from the run form: a CSV with
CRLF line endings, a header from the form's columns, and one row per file.
"""
import importlib.util
import io
import unittest
from pathlib import Path


def _load_preprocess():
    path = Path(__file__).with_name("preprocess.py")
    spec = importlib.util.spec_from_file_location("sourmash_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()

HEADER = "sample,file\r\n"
REFERENCES = "s3://refs-bucket"


class ReadSamplesheet(unittest.TestCase):

    def test_rows_of_one_sample_are_read_as_text(self):
        samplesheet = preprocess.read_samplesheet(io.StringIO(
            HEADER
            + "01,s3://bucket/data/01_R1.fastq.gz\r\n"
            + "01,s3://bucket/data/01_R2.fastq.gz\r\n"
        ))
        self.assertEqual(samplesheet["sample"].tolist(), ["01", "01"])

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


class FormatDb(unittest.TestCase):

    def test_selected_databases_are_globbed_at_the_ksize(self):
        db = preprocess.format_db(
            bact_db="85k: GTDB R08-RS214 bacterial representatives",
            viral_db="48k: NCBI GenBank (2022) - viruses",
            archaea_db="None",
            protozoa_db="None",
            fungi_db="None",
            ksize="21",
            prefix=f"{REFERENCES}/",
            threshold_bp=50000,
        )
        self.assertEqual(
            db,
            "s3://refs-bucket/{sourmash/gtdb-rs214-reps.k21,genbank/genbank-2022.03-viral-k21}.zip",
        )


class CheckThreshold(unittest.TestCase):

    def test_zero_and_positive_pass(self):
        for threshold_bp in (0, 50000):
            preprocess.check_threshold(threshold_bp)

    def test_negative_fails_with_the_value(self):
        with self.assertRaisesRegex(ValueError, r"cannot be negative \(-1\)"):
            preprocess.check_threshold(-1)


if __name__ == "__main__":
    unittest.main()
