"""Offline tests for the anvi'o pangenome preprocess script.

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
    spec = importlib.util.spec_from_file_location("anvio_pangenome_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()

HEADER = "genome,name,category\r\n"


def read(rows: str):
    return preprocess.read_samplesheet(io.StringIO(HEADER + rows))


class ReadSamplesheet(unittest.TestCase):

    def test_blank_cells_stay_empty_text(self):
        samplesheet = read("s3://bucket/data/a.fna.gz,Strain_A,\r\n")
        self.assertEqual(samplesheet["category"].tolist(), [""])

    def test_no_param_fails_clearly(self):
        for path in (None, ""):
            with self.assertRaisesRegex(ValueError, "No genomes were provided"):
                preprocess.read_samplesheet(path)

    def test_empty_file_fails_clearly(self):
        with self.assertRaisesRegex(ValueError, "The genome table is empty"):
            preprocess.read_samplesheet(io.StringIO(""))

    def test_header_without_rows_fails_clearly(self):
        with self.assertRaisesRegex(ValueError, "The genome table is empty"):
            read("")


class CheckUniqueNames(unittest.TestCase):

    def test_distinct_names_pass(self):
        preprocess.check_unique_names(read(
            "s3://bucket/data/a.fna,Strain_A,\r\ns3://bucket/data/b.fna,Strain_B,\r\n"
        ))

    def test_repeated_name_fails_naming_it(self):
        with self.assertRaisesRegex(ValueError, "listed only once.*: Strain_A$"):
            preprocess.check_unique_names(read(
                "s3://bucket/data/a.fna,Strain_A,\r\ns3://bucket/data/b.fna,Strain_A,\r\n"
            ))


class DropBlankCategory(unittest.TestCase):

    def test_blank_category_is_dropped(self):
        samplesheet = preprocess.drop_blank_category(read(
            "s3://bucket/data/a.fna,Strain_A,\r\ns3://bucket/data/b.fna,Strain_B,\r\n"
        ))
        self.assertEqual(samplesheet.columns.tolist(), ["genome", "name"])

    def test_filled_category_is_kept(self):
        samplesheet = preprocess.drop_blank_category(read(
            "s3://bucket/data/a.fna,Strain_A,gut\r\ns3://bucket/data/b.fna,Strain_B,skin\r\n"
        ))
        self.assertEqual(samplesheet["category"].tolist(), ["gut", "skin"])

    def test_single_category_value_fails(self):
        with self.assertRaisesRegex(ValueError, "at least two distinct values"):
            preprocess.drop_blank_category(read(
                "s3://bucket/data/a.fna,Strain_A,gut\r\ns3://bucket/data/b.fna,Strain_B,gut\r\n"
            ))


if __name__ == "__main__":
    unittest.main()
