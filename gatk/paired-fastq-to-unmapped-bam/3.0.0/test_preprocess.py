"""Offline tests for the GATK paired FASTQ to uBAM 3.0.0 preprocess script.

preprocess.py imports the Cirro SDK only when it runs, so this is a plain
`python -m unittest` that needs pandas and nothing else.

The samplesheet is the structured input Cirro writes from the run form: a CSV with
CRLF line endings, a header from the form's columns, and blank cells left empty.
"""
import importlib.util
import io
import unittest
from pathlib import Path


def _load_preprocess():
    path = Path(__file__).with_name("preprocess.py")
    spec = importlib.util.spec_from_file_location("paired_fastq_to_ubam_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()

HEADER = "sample_name,readgroup_name,library_name,platform_unit,fastq_1,fastq_2\r\n"
WF = "ConvertPairedFastQsToUnmappedBamWf"


def _samplesheet(rows):
    return preprocess.read_samplesheet(io.StringIO(HEADER + rows))


class ReadSamplesheet(unittest.TestCase):

    def test_values_stay_text(self):
        samplesheet = _samplesheet("01,,,,s3://b/data/01_R1.fastq.gz,s3://b/data/01_R2.fastq.gz\r\n")
        self.assertEqual(samplesheet["sample_name"].tolist(), ["01"])
        self.assertEqual(samplesheet["readgroup_name"].tolist(), [""])

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


class FillReadGroups(unittest.TestCase):

    def test_blank_fields_take_the_sample_name(self):
        samplesheet = preprocess.fill_read_groups(_samplesheet(
            "A,A.L1,libA,,s3://b/data/A_L1_R1.fq.gz,s3://b/data/A_L1_R2.fq.gz\r\n"
            + "B,,,FC1.2,s3://b/data/B_R1.fq.gz,s3://b/data/B_R2.fq.gz\r\n"
        ))
        self.assertEqual(samplesheet["readgroup_name"].tolist(), ["A.L1", "B"])
        self.assertEqual(samplesheet["library_name"].tolist(), ["libA", "B"])
        self.assertEqual(samplesheet["platform_unit"].tolist(), ["A", "FC1.2"])

    def test_repeated_read_group_fails_clearly(self):
        samplesheet = _samplesheet(
            "A,,,,s3://b/data/A_L1_R1.fq.gz,s3://b/data/A_L1_R2.fq.gz\r\n"
            + "A,,,,s3://b/data/A_L2_R1.fq.gz,s3://b/data/A_L2_R2.fq.gz\r\n"
        )
        with self.assertRaisesRegex(ValueError, "used more than once: A"):
            preprocess.fill_read_groups(samplesheet)


class SharedInputs(unittest.TestCase):

    def test_spaces_in_sequencing_center_are_replaced(self):
        shared = preprocess.shared_inputs({f"{WF}.sequencing_center": "Fred Hutch"}, "2026-01-01")
        self.assertEqual(shared[f"{WF}.sequencing_center"], "Fred_Hutch")

    def test_blank_run_date_becomes_today(self):
        for params in ({}, {f"{WF}.run_date": ""}, {f"{WF}.run_date": None}):
            shared = preprocess.shared_inputs(params, "2026-01-01")
            self.assertEqual(shared[f"{WF}.run_date"], "2026-01-01")

    def test_given_run_date_is_kept_and_options_are_left_out(self):
        shared = preprocess.shared_inputs(
            {f"{WF}.run_date": "2025-05-05", "samplesheet": "s3://b/config/x.csv"}, "2026-01-01"
        )
        self.assertEqual(shared[f"{WF}.run_date"], "2025-05-05")
        self.assertNotIn("samplesheet", shared)


class BuildInputs(unittest.TestCase):

    def test_one_input_set_per_row(self):
        samplesheet = preprocess.fill_read_groups(_samplesheet(
            "A,A.L1,,,s3://b/data/A_L1_R1.fq.gz,s3://b/data/A_L1_R2.fq.gz\r\n"
            + "A,A.L2,,,s3://b/data/A_L2_R1.fq.gz,s3://b/data/A_L2_R2.fq.gz\r\n"
        ))
        all_inputs = preprocess.build_inputs(samplesheet, {f"{WF}.platform_name": "Illumina"})
        self.assertEqual(len(all_inputs), 2)
        self.assertEqual(all_inputs[1], {
            f"{WF}.sample_name": "A",
            f"{WF}.readgroup_name": "A.L2",
            f"{WF}.library_name": "A",
            f"{WF}.platform_unit": "A",
            f"{WF}.fastq_1": "s3://b/data/A_L2_R1.fq.gz",
            f"{WF}.fastq_2": "s3://b/data/A_L2_R2.fq.gz",
            f"{WF}.platform_name": "Illumina",
        })


if __name__ == "__main__":
    unittest.main()
