"""Offline tests for the GATK HaplotypeCaller 2.3.1 preprocess script.

preprocess.py imports the Cirro SDK only when it runs, so this is a plain
`python -m unittest` that needs pandas and nothing else.

The samplesheet is the structured input Cirro writes from the run form: a CSV with
CRLF line endings, a header from the form's columns, and one row per sample.
"""
import importlib.util
import io
import unittest
from pathlib import Path


def _load_preprocess():
    path = Path(__file__).with_name("preprocess.py")
    spec = importlib.util.spec_from_file_location("haplotypecaller_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()

HEADER = "input_bam,input_bam_index\r\n"
ROWS = (
    "s3://bucket/data/A.bam,s3://bucket/data/A.bai\r\n"
    + "s3://bucket/data/B.cram,s3://bucket/data/B.cram.crai\r\n"
)


class ReadSamplesheet(unittest.TestCase):

    def test_rows_are_read(self):
        samplesheet = preprocess.read_samplesheet(io.StringIO(HEADER + ROWS))
        self.assertEqual(samplesheet["input_bam"].tolist(), ["s3://bucket/data/A.bam", "s3://bucket/data/B.cram"])

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


class BuildInputs(unittest.TestCase):

    def setUp(self):
        self.params = {
            "samplesheet": "s3://bucket/config/pipeline_samplesheet.csv",
            "HaplotypeCallerGvcf_GATK4.ref_fasta": "s3://refs/ref.fasta",
            "final_workflow_outputs_dir": "s3://bucket/data/",
        }
        self.inputs = preprocess.build_inputs(
            preprocess.read_samplesheet(io.StringIO(HEADER + ROWS)), self.params
        )

    def test_one_input_set_per_row(self):
        self.assertEqual(len(self.inputs), 2)
        self.assertEqual(self.inputs[1]["HaplotypeCallerGvcf_GATK4.input_bam"], "s3://bucket/data/B.cram")
        self.assertEqual(self.inputs[1]["HaplotypeCallerGvcf_GATK4.input_bam_index"], "s3://bucket/data/B.cram.crai")

    def test_workflow_params_are_shared_and_others_are_left_out(self):
        for single_input in self.inputs:
            self.assertEqual(single_input["HaplotypeCallerGvcf_GATK4.ref_fasta"], "s3://refs/ref.fasta")
            self.assertNotIn("samplesheet", single_input)
            self.assertNotIn("final_workflow_outputs_dir", single_input)

    def test_task_settings_are_set(self):
        self.assertEqual(self.inputs[0]["HaplotypeCallerGvcf_GATK4.HaplotypeCaller.disk_space_gb"], 100)
        self.assertEqual(self.inputs[0]["HaplotypeCallerGvcf_GATK4.MergeGVCFs.preemptible_attempts"], 3)


if __name__ == "__main__":
    unittest.main()
