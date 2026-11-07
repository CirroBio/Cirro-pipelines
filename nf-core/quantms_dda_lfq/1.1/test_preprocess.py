"""Offline tests for the nf-core/quantms DDA-LFQ 1.1 SDRF handling.

preprocess.py imports the Cirro SDK only when it runs, so this is a plain
`python -m unittest` that needs pandas and nothing else.

The table is the structured input Cirro writes from the run form: a TSV with
CRLF line endings, a header from the form's columns, and blank cells left empty.
"""
import importlib.util
import io
import re
import unittest
from pathlib import Path

import pandas as pd


def _load_preprocess():
    path = Path(__file__).with_name("preprocess.py")
    spec = importlib.util.spec_from_file_location("quantms_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()

COLUMNS = [
    "source name",
    "comment[file uri]",
    "factor value[condition]",
    "comment[fraction identifier]",
    "comment[technical replicate]",
    "assay name",
    "characteristics[organism]",
    "characteristics[organism part]",
    "characteristics[disease]",
    "characteristics[cell type]",
    "characteristics[biological replicate]",
]
HEADER = "\t".join(COLUMNS) + "\r\n"
URI = "s3://bucket/datasets/ds1/data"


def table(*rows: dict) -> pd.DataFrame:
    """The table as read_samplesheet returns it: every column, with blanks as ''."""
    return pd.DataFrame([{column: row.get(column, "") for column in COLUMNS} for row in rows])


def make_sdrf(samplesheet: pd.DataFrame) -> pd.DataFrame:
    return preprocess.make_sdrf(
        samplesheet,
        organism="Homo sapiens",
        enzyme="Trypsin",
        fragment_tolerance="0.6 Da",
        precursor_tolerance="10 ppm",
    )


class ReadSamplesheet(unittest.TestCase):

    def test_values_stay_text(self):
        samplesheet = preprocess.read_samplesheet(io.StringIO(
            HEADER + "01\t" + URI + "/a.raw" + "\t" * (len(COLUMNS) - 2) + "\r\n"
        ))
        self.assertEqual(samplesheet.columns.tolist(), COLUMNS)
        self.assertEqual(samplesheet["source name"].tolist(), ["01"])
        self.assertEqual(samplesheet["characteristics[biological replicate]"].tolist(), [""])

    def test_no_param_fails_clearly(self):
        for path in (None, ""):
            with self.assertRaisesRegex(ValueError, "No SDRF was provided"):
                preprocess.read_samplesheet(path)

    def test_empty_file_fails_clearly(self):
        with self.assertRaisesRegex(ValueError, "The SDRF is empty"):
            preprocess.read_samplesheet(io.StringIO(""))

    def test_header_without_rows_fails_clearly(self):
        with self.assertRaisesRegex(ValueError, "The SDRF is empty"):
            preprocess.read_samplesheet(io.StringIO(HEADER))


class MakeSdrf(unittest.TestCase):

    def test_blanks_are_filled_and_run_columns_added(self):
        sdrf = make_sdrf(table(
            {"source name": "A", "comment[file uri]": f"{URI}/A_1.raw"},
            {"source name": "A", "comment[file uri]": f"{URI}/sub/A_2.raw"},
            {"source name": "B", "comment[file uri]": f"{URI}/B_1.mzML"},
        ))
        self.assertEqual(sdrf["comment[data file]"].tolist(), ["A_1.raw", "A_2.raw", "B_1.mzML"])
        self.assertEqual(sdrf["assay name"].tolist(), ["A", "A", "B"])
        self.assertEqual(sdrf["characteristics[organism]"].tolist(), ["Homo sapiens"] * 3)
        self.assertEqual(sdrf["characteristics[disease]"].tolist(), ["not applicable"] * 3)
        self.assertEqual(sdrf["comment[fraction identifier]"].tolist(), ["1", "1", "1"])
        self.assertEqual(sdrf["comment[technical replicate]"].tolist(), ["1", "2", "1"])
        self.assertEqual(sdrf["comment[label]"].unique().tolist(), ["NT=label free sample"])
        self.assertEqual(sdrf["comment[cleavage agent details]"].unique().tolist(), ["NT=Trypsin"])
        self.assertEqual(sdrf["comment[fragment mass tolerance]"].unique().tolist(), ["0.6 Da"])
        self.assertEqual(sdrf["comment[precursor mass tolerance]"].unique().tolist(), ["10 ppm"])
        self.assertNotIn("factor value[condition]", sdrf.columns)
        self.assertFalse((sdrf == "").any().any())

    def test_set_values_are_kept(self):
        sdrf = make_sdrf(table({
            "source name": "Sample 1",
            "comment[file uri]": f"{URI}/a.raw",
            "assay name": "run 1",
            "comment[technical replicate]": "3",
            "characteristics[organism]": "Mus musculus",
            "characteristics[disease]": "normal",
        }))
        row = sdrf.iloc[0]
        self.assertEqual(row["assay name"], "run 1")
        self.assertEqual(row["comment[technical replicate]"], "3")
        self.assertEqual(row["characteristics[organism]"], "Mus musculus")
        self.assertEqual(row["characteristics[disease]"], "normal")

    def test_fractions_of_a_sample_share_a_technical_replicate(self):
        sdrf = make_sdrf(table(*[
            {"source name": f"BSA{i}", "comment[file uri]": f"{URI}/BSA{i}_F{f}.mzML", "comment[fraction identifier]": f}
            for i in ("1", "2") for f in ("1", "2")
        ]))
        self.assertEqual(sdrf["comment[technical replicate]"].tolist(), ["1", "1", "1", "1"])
        self.assertEqual(sdrf["comment[fraction identifier]"].tolist(), ["1", "2", "1", "2"])

    def test_whitespace_is_stripped(self):
        sdrf = make_sdrf(table({"source name": " A ", "comment[file uri]": f"{URI}/a.raw "}))
        self.assertEqual(sdrf["source name"].tolist(), ["A"])
        self.assertEqual(sdrf["comment[data file]"].tolist(), ["a.raw"])

    def test_columns_are_in_the_order_the_sdrf_schema_requires(self):
        sdrf = make_sdrf(table(
            {"source name": "A", "comment[file uri]": f"{URI}/a.raw", "factor value[condition]": "treated"},
        ))
        columns = sdrf.columns.tolist()
        assay = columns.index("assay name")
        self.assertEqual(columns[0], "source name")
        self.assertEqual(columns[-1], "factor value[condition]")
        for i, column in enumerate(columns):
            if column.startswith("characteristics"):
                self.assertLess(i, assay, column)
            if column.startswith("comment"):
                self.assertGreater(i, assay, column)
        # sdrf-pipelines 0.0.24 (quantms 1.2.0) requires these, with or without --CHECK_MS
        for column in [
            "characteristics[organism]", "characteristics[organism part]", "characteristics[disease]",
            "characteristics[cell type]", "characteristics[biological replicate]",
            "comment[technical replicate]", "comment[fraction identifier]", "comment[data file]",
            "comment[label]", "comment[instrument]", "comment[cleavage agent details]",
        ]:
            self.assertIn(column, columns)
        # Every name must be one the SDRF schema accepts
        for column in columns:
            self.assertTrue(
                column in ("source name", "assay name")
                or re.fullmatch(r"(characteristics|comment|factor value)\[[^\]]+\]", column),
                column,
            )

    def test_enzyme_is_the_cleavage_agent(self):
        sdrf = preprocess.make_sdrf(
            table({"source name": "A", "comment[file uri]": f"{URI}/a.raw"}),
            organism="Homo sapiens", enzyme="Lys-C", fragment_tolerance="0.6 Da", precursor_tolerance="10 ppm",
        )
        self.assertEqual(sdrf["comment[cleavage agent details]"].tolist(), ["NT=Lys-C"])

    def test_condition_in_some_rows_only_fails(self):
        with self.assertRaisesRegex(ValueError, r"in every row or none.*row\(s\) 2"):
            make_sdrf(table(
                {"source name": "A", "comment[file uri]": f"{URI}/a.raw", "factor value[condition]": "treated"},
                {"source name": "B", "comment[file uri]": f"{URI}/b.raw"},
            ))

    def test_blank_source_name_or_file_fails(self):
        with self.assertRaisesRegex(ValueError, r"needs a source name.*row\(s\) 2"):
            make_sdrf(table(
                {"source name": "A", "comment[file uri]": f"{URI}/a.raw"},
                {"comment[file uri]": f"{URI}/b.raw"},
            ))
        with self.assertRaisesRegex(ValueError, r"needs a comment\[file uri\].*row\(s\) 1"):
            make_sdrf(table({"source name": "A"}))

    def test_repeated_file_name_fails(self):
        with self.assertRaisesRegex(ValueError, "repeated: a.raw"):
            make_sdrf(table(
                {"source name": "A", "comment[file uri]": f"{URI}/a.raw"},
                {"source name": "B", "comment[file uri]": "s3://bucket/datasets/ds2/data/a.raw"},
            ))

    def test_non_integer_replicate_fails(self):
        with self.assertRaisesRegex(ValueError, r"comment\[technical replicate\] must be a whole number.*not: one"):
            make_sdrf(table(
                {"source name": "A", "comment[file uri]": f"{URI}/a.raw", "comment[technical replicate]": "one"},
            ))


if __name__ == "__main__":
    unittest.main()
