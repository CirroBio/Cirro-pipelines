"""Offline tests for the GATK processing-for-variant-discovery 2.1.1 preprocess script.

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
    spec = importlib.util.spec_from_file_location("processing_for_variant_discovery_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()

HEADER = "sample_name,flowcell_unmapped_bams\r\n"
WF = "PreProcessingForVariantDiscovery_GATK4"
ROWS = (
    "B,s3://b/data/B.L1.unmapped.bam\r\n"
    + "A,s3://b/data/A.L1.unmapped.bam\r\n"
    + "B,s3://b/data/B.L2.unmapped.bam\r\n"
)


class ReadSamplesheet(unittest.TestCase):

    def test_values_stay_text(self):
        samplesheet = preprocess.read_samplesheet(io.StringIO(HEADER + "01,s3://b/data/01.unmapped.bam\r\n"))
        self.assertEqual(samplesheet["sample_name"].tolist(), ["01"])

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


class GroupSamples(unittest.TestCase):

    def test_rows_are_grouped_by_sample_in_table_order(self):
        samples = preprocess.group_samples(preprocess.read_samplesheet(io.StringIO(HEADER + ROWS)))
        self.assertEqual(samples, [
            ("B", ["s3://b/data/B.L1.unmapped.bam", "s3://b/data/B.L2.unmapped.bam"]),
            ("A", ["s3://b/data/A.L1.unmapped.bam"]),
        ])


class KnownIndels(unittest.TestCase):

    def test_lists_are_under_the_references_bucket(self):
        indels = preprocess.known_indels("s3://refs")
        self.assertEqual(indels[f"{WF}.known_indels_sites_VCFs"], [
            "s3://refs/GATK/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
            "s3://refs/GATK/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz",
        ])
        self.assertEqual(indels[f"{WF}.known_indels_sites_indices"], [
            "s3://refs/GATK/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi",
            "s3://refs/GATK/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi",
        ])


class BuildInputs(unittest.TestCase):

    def test_one_input_set_per_sample(self):
        samples = preprocess.group_samples(preprocess.read_samplesheet(io.StringIO(HEADER + ROWS)))
        all_inputs = preprocess.build_inputs(
            samples,
            ["s3://b/config/flowcell_unmapped_bams.0.list", "s3://b/config/flowcell_unmapped_bams.1.list"],
            {f"{WF}.ref_name": "hg38", "samplesheet": "s3://b/config/pipeline_samplesheet.csv"},
        )
        self.assertEqual(all_inputs, [
            {
                f"{WF}.sample_name": "B",
                f"{WF}.flowcell_unmapped_bams_list": "s3://b/config/flowcell_unmapped_bams.0.list",
                f"{WF}.ref_name": "hg38",
            },
            {
                f"{WF}.sample_name": "A",
                f"{WF}.flowcell_unmapped_bams_list": "s3://b/config/flowcell_unmapped_bams.1.list",
                f"{WF}.ref_name": "hg38",
            },
        ])


if __name__ == "__main__":
    unittest.main()
