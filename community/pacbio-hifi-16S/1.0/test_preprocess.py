"""Offline tests for the HiFi 16S preprocess script.

preprocess.py imports the Cirro SDK only when it runs, so this is a plain
`python -m unittest` that needs pandas and nothing else.

The tables are the structured inputs Cirro writes from the run form: TSVs with
CRLF line endings and a header from the form's columns.
"""
import importlib.util
import io
import unittest
from pathlib import Path


def _load_preprocess():
    path = Path(__file__).with_name("preprocess.py")
    spec = importlib.util.spec_from_file_location("pacbio_hifi_16s_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()

SAMPLES = (
    "sample-id\tabsolute-filepath\r\n"
    "01\ts3://bucket/data/01.hifi_reads.fastq.gz\r\n"
    "02\ts3://bucket/data/02.hifi_reads.fastq.gz\r\n"
    "03\ts3://bucket/data/03.hifi_reads.fastq.gz\r\n"
)
METADATA_HEADER = "sample_name\tcondition\tpool\r\n"


def samples():
    return preprocess.read_table(io.StringIO(SAMPLES))


class ReadTable(unittest.TestCase):

    def test_sample_ids_are_read_as_text(self):
        self.assertEqual(samples()["sample-id"].tolist(), ["01", "02", "03"])

    def test_missing_or_empty_table_reads_as_empty(self):
        for source in (None, "", io.StringIO(""), io.StringIO(METADATA_HEADER)):
            self.assertTrue(preprocess.read_table(source).empty)


class CheckSamples(unittest.TestCase):

    def test_empty_table_fails_clearly(self):
        with self.assertRaisesRegex(ValueError, "The samples table is empty"):
            preprocess.check_samples(preprocess.read_table(None))

    def test_repeated_sample_fails(self):
        table = preprocess.read_table(io.StringIO(SAMPLES + "01\ts3://bucket/data/01b.fastq.gz\r\n"))
        with self.assertRaisesRegex(ValueError, "repeated: 01"):
            preprocess.check_samples(table)

    def test_unique_samples_pass(self):
        preprocess.check_samples(samples())


class CompleteMetadata(unittest.TestCase):

    def test_empty_metadata_gets_placeholder_conditions(self):
        metadata = preprocess.complete_metadata(preprocess.read_table(None), samples())
        self.assertEqual(metadata["sample_name"].tolist(), ["01", "02", "03"])
        self.assertEqual(metadata["condition"].tolist(), ["cat_0", "cat_1", "cat_0"])

    def test_sample_without_metadata_fails(self):
        metadata = preprocess.read_table(io.StringIO(METADATA_HEADER + "01\tA\t\r\n02\tB\t\r\n"))
        with self.assertRaisesRegex(ValueError, "no row for these samples: 03"):
            preprocess.complete_metadata(metadata, samples())

    def test_blank_pool_column_is_dropped(self):
        metadata = preprocess.read_table(io.StringIO(METADATA_HEADER + "01\tA\t\r\n02\tB\t\r\n03\tB\t\r\n"))
        metadata = preprocess.complete_metadata(metadata, samples())
        self.assertEqual(list(metadata.columns), ["sample_name", "condition"])

    def test_pools_are_kept(self):
        metadata = preprocess.read_table(io.StringIO(METADATA_HEADER + "01\tA\tp1\r\n02\tB\tp1\r\n03\tB\tp2\r\n"))
        metadata = preprocess.complete_metadata(metadata, samples())
        self.assertEqual(metadata["pool"].tolist(), ["p1", "p1", "p2"])


if __name__ == "__main__":
    unittest.main()
