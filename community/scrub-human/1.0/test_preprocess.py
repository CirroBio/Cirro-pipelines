"""Offline tests for the Scrub Human Sequences preprocess script.

preprocess.py imports the Cirro SDK only when it runs, so this is a plain
`python -m unittest` that needs pandas and nothing else.

The table is the structured input Cirro writes from the run form: a CSV with CRLF
line endings and a header from the form's columns.
"""
import importlib.util
import io
import unittest
from pathlib import Path


def _load_preprocess():
    path = Path(__file__).with_name("preprocess.py")
    spec = importlib.util.spec_from_file_location("scrub_human_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()

HEADER = "fastq_1,fastq_2\r\n"


class ReadSamplesheet(unittest.TestCase):

    def test_no_param_fails_clearly(self):
        for path in (None, ""):
            with self.assertRaisesRegex(ValueError, "No FASTQ files were provided"):
                preprocess.read_samplesheet(path)

    def test_empty_file_fails_clearly(self):
        with self.assertRaisesRegex(ValueError, "The table of FASTQ files is empty"):
            preprocess.read_samplesheet(io.StringIO(""))

    def test_header_without_rows_fails_clearly(self):
        with self.assertRaisesRegex(ValueError, "The table of FASTQ files is empty"):
            preprocess.read_samplesheet(io.StringIO(HEADER))


class FileList(unittest.TestCase):

    def test_one_path_per_line_ending_in_a_newline(self):
        samplesheet = preprocess.read_samplesheet(io.StringIO(
            HEADER
            + "s3://bucket/data/A_R1.fastq.gz,s3://bucket/data/A_R2.fastq.gz\r\n"
            + "s3://other/data/B_R1.fastq.gz,\r\n"
        ))
        written = io.StringIO()
        preprocess.file_list(samplesheet).to_csv(written, header=False, index=False)
        self.assertEqual(
            written.getvalue(),
            "s3://bucket/data/A_R1.fastq.gz\n"
            "s3://bucket/data/A_R2.fastq.gz\n"
            "s3://other/data/B_R1.fastq.gz\n",
        )


if __name__ == "__main__":
    unittest.main()
