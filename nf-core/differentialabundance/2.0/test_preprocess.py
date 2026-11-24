"""Offline tests for the nf-core/differentialabundance 2.0 preprocess script.

preprocess.py imports the Cirro SDK and boto3 only when it runs, so this is a plain
`python -m unittest` that needs pandas and nothing else.

The contrasts are the structured input Cirro writes from the run form: a CSV with CRLF
line endings, a header from the form's columns, and one row per contrast. The sample
metadata is a CSV or TSV the user picks, with any columns.
"""
import importlib.util
import io
import unittest
from pathlib import Path

import pandas as pd


def _load_preprocess():
    path = Path(__file__).with_name("preprocess.py")
    spec = importlib.util.spec_from_file_location("differentialabundance_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()

HEADER = "id,variable,reference,target,blocking\r\n"
RNASEQ_RESULTS = "s3://bucket/datasets/rnaseq/data/"


def contrasts(*rows: str) -> pd.DataFrame:
    return preprocess.read_contrasts(io.StringIO(HEADER + "".join(f"{row}\r\n" for row in rows)))


SAMPLE_METADATA = pd.read_csv(
    io.StringIO("sample,condition,batch\nA,control,1\nB,treated,1\nC,control,2\nD,treated,2\n"),
    dtype=str,
    keep_default_na=False,
)


class ReadContrasts(unittest.TestCase):

    def test_rows_are_read_as_text(self):
        table = contrasts("treated_vs_control,condition,control,treated,")
        self.assertEqual(table["blocking"].tolist(), [""])

    def test_no_param_fails_clearly(self):
        for path in (None, ""):
            with self.assertRaisesRegex(ValueError, "No contrasts were provided"):
                preprocess.read_contrasts(path)

    def test_empty_file_fails_clearly(self):
        with self.assertRaisesRegex(ValueError, "The contrasts table is empty"):
            preprocess.read_contrasts(io.StringIO(""))

    def test_header_without_rows_fails_clearly(self):
        with self.assertRaisesRegex(ValueError, "The contrasts table is empty"):
            preprocess.read_contrasts(io.StringIO(HEADER))


class CheckContrasts(unittest.TestCase):

    def test_matching_contrasts_pass(self):
        preprocess.check_contrasts(
            contrasts("treated_vs_control,condition,control,treated,", "blocked,condition,control,treated,batch"),
            SAMPLE_METADATA,
        )

    def test_repeated_id_fails(self):
        with self.assertRaisesRegex(ValueError, "repeated: c1"):
            preprocess.check_contrasts(
                contrasts("c1,condition,control,treated,", "c1,batch,1,2,"),
                SAMPLE_METADATA,
            )

    def test_every_mismatch_is_reported(self):
        with self.assertRaises(ValueError) as raised:
            preprocess.check_contrasts(
                contrasts("c1,genotype,wt,ko,", "c2,condition,control,mutant,", "c3,condition,control,treated,donor"),
                SAMPLE_METADATA,
            )
        message = str(raised.exception)
        self.assertIn("c1: no column genotype", message)
        self.assertIn("c2: no sample has condition = mutant", message)
        self.assertIn("c3: no blocking column donor", message)


class FindMatrices(unittest.TestCase):

    def test_first_published_matrix_is_used_with_its_lengths(self):
        published = {
            "s3://bucket/datasets/rnaseq/data/star_salmon/salmon.merged.gene_counts.tsv",
            "s3://bucket/datasets/rnaseq/data/star_salmon/salmon.merged.gene_lengths.tsv",
            "s3://bucket/datasets/rnaseq/data/star_rsem/rsem.merged.gene_counts.tsv",
        }
        self.assertEqual(
            preprocess.find_matrices(RNASEQ_RESULTS, published.__contains__),
            (
                "s3://bucket/datasets/rnaseq/data/star_salmon/salmon.merged.gene_counts.tsv",
                "s3://bucket/datasets/rnaseq/data/star_salmon/salmon.merged.gene_lengths.tsv",
            ),
        )

    def test_missing_lengths_are_left_out(self):
        published = {"s3://bucket/datasets/rnaseq/data/kallisto/kallisto.merged.gene_counts.tsv"}
        self.assertEqual(
            preprocess.find_matrices(RNASEQ_RESULTS, published.__contains__),
            ("s3://bucket/datasets/rnaseq/data/kallisto/kallisto.merged.gene_counts.tsv", None),
        )

    def test_no_matrix_fails_clearly(self):
        with self.assertRaisesRegex(ValueError, "No gene count matrix found in s3://bucket/datasets/rnaseq/data/ "):
            preprocess.find_matrices(RNASEQ_RESULTS, lambda uri: False)


class FindRnaseqDataset(unittest.TestCase):

    INPUTS = [
        {"id": "metadata", "s3": "s3://bucket/datasets/rnaseq-metadata"},
        {"id": "rnaseq", "s3": "s3://bucket/datasets/rnaseq"},
    ]

    def test_dataset_holding_the_picked_folder_is_found(self):
        self.assertEqual(preprocess.find_rnaseq_dataset(self.INPUTS, RNASEQ_RESULTS)["id"], "rnaseq")

    def test_no_matching_dataset_fails(self):
        with self.assertRaisesRegex(ValueError, "None of the run's input datasets"):
            preprocess.find_rnaseq_dataset(self.INPUTS[:1], RNASEQ_RESULTS)


class FindIgenomesGenome(unittest.TestCase):

    def test_igenomes_selection_of_the_chosen_aligner(self):
        params = {"aligner_and_reference_genome": {"reference_genome": {
            "aligner": "star_salmon",
            "star_salmon_genome_selection": {"genome_source": "igenomes", "igenomes": {"genome": "GRCm38"}},
        }}}
        self.assertEqual(preprocess.find_igenomes_genome(params), "GRCm38")

    def test_cirro_genome_index_has_no_igenomes_genome(self):
        params = {"aligner_and_reference_genome": {"reference_genome": {
            "aligner": "star_salmon",
            "star_salmon_genome_selection": {"genome_source": "dataset", "igenomes": {"genome": "GRCh38"}},
        }}}
        self.assertIsNone(preprocess.find_igenomes_genome(params))


if __name__ == "__main__":
    unittest.main()
