"""Offline tests for the nf-core/ampliseq 2.9.0 preprocess script.

preprocess.py imports the Cirro SDK only when it runs, so this is a plain
`python -m unittest` that needs pandas and nothing else.

The samplesheet is the structured input Cirro writes from the run form: a TSV with
CRLF line endings, a header from the form's columns, and blank cells left empty.
"""
import importlib.util
import io
import unittest
from pathlib import Path


def _load_preprocess():
    path = Path(__file__).with_name("preprocess.py")
    spec = importlib.util.spec_from_file_location("ampliseq_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()

HEADER = "sampleID\tforwardReads\treverseReads\trun\r\n"
DATA = "s3://bucket/datasets/ds1/data"


def samplesheet(*rows):
    return preprocess.read_samplesheet(io.StringIO(HEADER + "".join(f"{row}\r\n" for row in rows)))


def paired(sample_id: str, run: str = "") -> str:
    return f"{sample_id}\t{DATA}/{sample_id}_R1_001.fastq.gz\t{DATA}/{sample_id}_R2_001.fastq.gz\t{run}"


def forward_only(sample_id: str) -> str:
    return f"{sample_id}\t{DATA}/{sample_id}_R1_001.fastq.gz\t\t"


class ReadSamplesheet(unittest.TestCase):

    def test_tab_separated_values_stay_text(self):
        table = samplesheet(paired("S01", run="01"))
        self.assertEqual(table.columns.tolist(), ["sampleID", "forwardReads", "reverseReads", "run"])
        self.assertEqual(table["run"].tolist(), ["01"])

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


class CheckUniqueSampleIds(unittest.TestCase):

    def test_distinct_ids_pass(self):
        preprocess.check_unique_sample_ids(samplesheet(paired("S1"), paired("S2")))

    def test_repeated_id_fails_naming_it(self):
        with self.assertRaisesRegex(ValueError, "repeated: S1$"):
            preprocess.check_unique_sample_ids(samplesheet(paired("S1"), paired("S2"), paired("S1", run="B")))


class IsSingleEnd(unittest.TestCase):

    def test_paired_rows_are_paired_end(self):
        self.assertFalse(preprocess.is_single_end(samplesheet(paired("S1"), paired("S2")), False))

    def test_forward_reads_only_option_makes_it_single_end(self):
        self.assertTrue(preprocess.is_single_end(samplesheet(paired("S1"), paired("S2")), True))

    def test_any_row_without_reverse_reads_makes_it_single_end(self):
        self.assertTrue(preprocess.is_single_end(samplesheet(paired("S1"), forward_only("S2")), False))


class CheckMetadataIds(unittest.TestCase):

    def test_every_sample_listed_passes(self):
        metadata = preprocess.read_metadata(io.StringIO(
            "ID\tcondition\r\nS1\tcontrol\r\nS2 \ttreatment\r\nS3\tcontrol\r\n"
        ))
        preprocess.check_metadata_ids(samplesheet(paired("S1"), paired("S2")), metadata)

    def test_missing_sample_fails_naming_it(self):
        metadata = preprocess.read_metadata(io.StringIO("ID\tcondition\r\nS1\tcontrol\r\n"))
        with self.assertRaisesRegex(ValueError, "Missing: S2$"):
            preprocess.check_metadata_ids(samplesheet(paired("S1"), paired("S2")), metadata)

    def test_ids_are_read_from_the_first_column(self):
        metadata = preprocess.read_metadata(io.StringIO("condition\tID\r\ncontrol\tS1\r\n"))
        with self.assertRaisesRegex(ValueError, "first column. Missing: S1$"):
            preprocess.check_metadata_ids(samplesheet(paired("S1")), metadata)


if __name__ == "__main__":
    unittest.main()
