"""Offline tests for the OCTAPUS preprocess script.

preprocess.py imports the Cirro SDK only when it runs, so this is a plain
`python -m unittest` that needs pandas and nothing else.

The genome table is the structured input Cirro writes from the run form: a CSV with
CRLF line endings, a header from the form's columns, and one row per genome.
"""
import importlib.util
import io
import unittest
from pathlib import Path


def _load_preprocess():
    path = Path(__file__).with_name("preprocess.py")
    spec = importlib.util.spec_from_file_location("octapus_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()

HEADER = "organism,uri\r\n"
ROW = "E. coli K-12,s3://bucket/data/k12.fna.gz\r\n"


class ReadGenomeTable(unittest.TestCase):

    def test_values_stay_text(self):
        genomes = preprocess.read_genome_table(io.StringIO(HEADER + "001,s3://bucket/data/001.fasta\r\n"))
        self.assertEqual(genomes["organism"].tolist(), ["001"])

    def test_no_param_fails_clearly(self):
        for path in (None, ""):
            with self.assertRaisesRegex(ValueError, "No genomes were provided"):
                preprocess.read_genome_table(path)

    def test_empty_file_fails_clearly(self):
        with self.assertRaisesRegex(ValueError, "The genome table is empty"):
            preprocess.read_genome_table(io.StringIO(""))

    def test_header_without_rows_fails_clearly(self):
        with self.assertRaisesRegex(ValueError, "The genome table is empty"):
            preprocess.read_genome_table(io.StringIO(HEADER))


class FormatManifest(unittest.TestCase):

    def test_columns_match_the_workflow(self):
        manifest = preprocess.format_manifest(preprocess.read_genome_table(io.StringIO(HEADER + ROW)))
        self.assertEqual(manifest.columns.tolist(), ["#Organism Name", "uri", "GenBank FTP"])
        self.assertEqual(manifest.iloc[0].tolist(), ["E. coli K-12", "s3://bucket/data/k12.fna.gz", ""])


if __name__ == "__main__":
    unittest.main()
