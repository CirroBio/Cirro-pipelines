"""Offline tests for the PGAP 1.0 preprocess script.

preprocess.py imports the Cirro SDK and boto3 only when it writes to S3, so this is
a plain `python -m unittest` that needs nothing beyond the standard library.
"""
import importlib.util
import unittest
from pathlib import Path


def _load_preprocess():
    path = Path(__file__).with_name("preprocess.py")
    spec = importlib.util.spec_from_file_location("pgap_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()

PARAMS = {
    "genus_species": "Escherichia coli",
    "last_name": "Doe",
    "first_name": "Jane",
    "email": "jane_doe@example.org",
    "organization": "NIH",
    "department": "NCBI",
    "phone": "301-555-0245",
    "street": "9000 Rockville Pike",
    "city": "Bethesda",
    "state": "MD",
    "postal_code": "20850",
    "country": "USA",
}

CONTACT_INFO = """contact_info:
    last_name: 'Doe'
    first_name: 'Jane'
    email: 'jane_doe@example.org'
    organization: 'NIH'
    department: 'NCBI'
    phone: '301-555-0245'
    street: '9000 Rockville Pike'
    city: 'Bethesda'
    state: 'MD'
    postal_code: '20850'
    country: 'USA'"""


class YamlQuote(unittest.TestCase):

    def test_plain_value_is_single_quoted(self):
        self.assertEqual(preprocess.yaml_quote("Escherichia coli"), "'Escherichia coli'")

    def test_quote_is_doubled(self):
        self.assertEqual(preprocess.yaml_quote("St. Mary's"), "'St. Mary''s'")


class MakeSubmol(unittest.TestCase):

    def test_required_fields_only(self):
        self.assertEqual(preprocess.make_submol(PARAMS), "\n".join([
            "organism:",
            "    genus_species: 'Escherichia coli'",
            CONTACT_INFO,
            "authors:",
            "    - author:",
            "        last_name: 'Doe'",
            "        first_name: 'Jane'",
        ]))

    def test_strain_and_middle_initial(self):
        submol = preprocess.make_submol({**PARAMS, "strain": "K-12", "middle_initial": "A"})
        self.assertEqual(submol, "\n".join([
            "organism:",
            "    genus_species: 'Escherichia coli'",
            "    strain: 'K-12'",
            CONTACT_INFO,
            "authors:",
            "    - author:",
            "        last_name: 'Doe'",
            "        first_name: 'Jane'",
            "        middle_initial: 'A'",
        ]))

    def test_blank_middle_initial_is_left_out(self):
        self.assertNotIn("middle_initial", preprocess.make_submol({**PARAMS, "middle_initial": ""}))

    def test_quotes_in_values_are_escaped(self):
        submol = preprocess.make_submol({**PARAMS, "last_name": "O'Brien", "organization": "St. Mary's"})
        self.assertIn("    organization: 'St. Mary''s'", submol.split("\n"))
        self.assertEqual(submol.count("last_name: 'O''Brien'"), 2)


class MakeSampleSheet(unittest.TestCase):

    def test_pairs_the_genome_with_its_yaml(self):
        self.assertEqual(
            preprocess.make_sample_sheet(
                "s3://bucket/datasets/a/data/genome.fna",
                "s3://bucket/datasets/b/config/submol.yaml",
            ),
            "fasta,yaml\n"
            "s3://bucket/datasets/a/data/genome.fna,s3://bucket/datasets/b/config/submol.yaml\n",
        )


if __name__ == "__main__":
    unittest.main()
