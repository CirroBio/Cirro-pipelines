"""Offline tests for the GATK whole-genome germline single-sample 3.1.6 preprocess script.

preprocess.py imports the Cirro SDK only when it runs, so this is a plain
`python -m unittest` that needs pandas and nothing else.

The samplesheet is the structured input Cirro writes from the run form: a CSV with
CRLF line endings, a header from the form's columns, and one row per unmapped BAM.
"""
import importlib.util
import io
import unittest
from pathlib import Path


def _load_preprocess():
    path = Path(__file__).with_name("preprocess.py")
    spec = importlib.util.spec_from_file_location("wgs_germline_single_sample_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()

HEADER = "sample,flowcell_unmapped_bams\r\n"
WF = "WholeGenomeGermlineSingleSample"


class ReadSamplesheet(unittest.TestCase):

    def test_values_stay_text(self):
        samplesheet = preprocess.read_samplesheet(io.StringIO(HEADER + "01,s3://b/data/01.unmapped.bam\r\n"))
        self.assertEqual(samplesheet["sample"].tolist(), ["01"])

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

    def test_one_input_set_per_sample_in_table_order(self):
        samplesheet = preprocess.read_samplesheet(io.StringIO(
            HEADER
            + "B,s3://b/data/B.L1.unmapped.bam\r\n"
            + "A,s3://b/data/A.L1.unmapped.bam\r\n"
            + "B,s3://b/data/B.L2.unmapped.bam\r\n"
        ))
        all_inputs = preprocess.build_inputs(samplesheet, {
            f"{WF}.allow_empty_ref_alt": True,
            "samplesheet": "s3://b/config/pipeline_samplesheet.csv",
        })
        self.assertEqual(all_inputs, [
            {
                f"{WF}.sample_and_unmapped_bams": {
                    "sample_name": "B",
                    "base_file_name": "B",
                    "final_gvcf_base_name": "B",
                    "unmapped_bam_suffix": ".bam",
                    "flowcell_unmapped_bams": ["s3://b/data/B.L1.unmapped.bam", "s3://b/data/B.L2.unmapped.bam"],
                },
                f"{WF}.allow_empty_ref_alt": True,
            },
            {
                f"{WF}.sample_and_unmapped_bams": {
                    "sample_name": "A",
                    "base_file_name": "A",
                    "final_gvcf_base_name": "A",
                    "unmapped_bam_suffix": ".bam",
                    "flowcell_unmapped_bams": ["s3://b/data/A.L1.unmapped.bam"],
                },
                f"{WF}.allow_empty_ref_alt": True,
            },
        ])


if __name__ == "__main__":
    unittest.main()
