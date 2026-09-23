"""Offline tests for reference-genome/aligner resolution in sarek_align.

Runs without the Cirro SDK, so this is a plain `python -m unittest`.

Regression coverage for: the iGenomes branch of genome_source dropped the
`aligner` param entirely (it was only restored on the Custom Genome branch),
so selecting bwa-mem2 or the GPU-accelerated Parabricks aligner in the form
silently had no effect -- sarek always fell back to its own default (bwa-mem).
"""
import sys
import types
import unittest
from pathlib import Path
from unittest.mock import patch


class FakeS3Client:
    def __init__(self):
        self.copies = []

    def copy(self, copy_source, bucket, key):
        self.copies.append((copy_source["Bucket"], copy_source["Key"], bucket, key))


def _load_preprocess():
    helpers = types.ModuleType("cirro.helpers.preprocess_dataset")
    helpers.PreprocessDataset = object

    s3_path_mod = types.ModuleType("cirro.models.s3_path")

    class S3Path:
        def __init__(self, path):
            self.valid = path.startswith("s3://")
            rest = path.replace("s3://", "", 1)
            self.bucket, _, self.key = rest.partition("/")

    s3_path_mod.S3Path = S3Path

    for name in ["cirro", "cirro.helpers", "cirro.models"]:
        sys.modules.setdefault(name, types.ModuleType(name))
    sys.modules["cirro.helpers.preprocess_dataset"] = helpers
    sys.modules["cirro.models.s3_path"] = s3_path_mod
    sys.modules.setdefault("pandas", types.ModuleType("pandas"))

    boto3_stub = types.ModuleType("boto3")
    boto3_stub.client = lambda service: FakeS3Client()
    sys.modules.setdefault("boto3", boto3_stub)

    import importlib.util
    path = Path(__file__).with_name("preprocess.py")
    spec = importlib.util.spec_from_file_location("sarek_align_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()


class FakeLogger:
    def info(self, *args):
        pass

    def warning(self, *args):
        pass


class FakeDataset:
    """The slice of PreprocessDataset that resolve_reference_genome touches."""

    def __init__(self, params):
        self.params = dict(params)
        self.logger = FakeLogger()

    def add_param(self, name, value, overwrite=False):
        if not overwrite and name in self.params:
            return
        self.params[name] = value

    def remove_param(self, name, force=False):
        self.params.pop(name, None)


class ResolveReferenceGenomeTests(unittest.TestCase):

    def test_igenomes_keeps_selected_aligner(self):
        for aligner in ["bwa-mem", "bwa-mem2", "parabricks"]:
            ds = FakeDataset({
                "genome_source": "igenomes",
                "genome": "GATK.GRCh38",
                "aligner": aligner,
            })
            preprocess.resolve_reference_genome(ds)
            self.assertEqual(ds.params.get("aligner"), aligner)

    def test_igenomes_defaults_to_bwa_mem_when_aligner_unset(self):
        ds = FakeDataset({"genome_source": "igenomes", "genome": "GATK.GRCh38"})
        preprocess.resolve_reference_genome(ds)
        self.assertEqual(ds.params.get("aligner"), "bwa-mem")

    def test_igenomes_does_not_leak_custom_genome_params(self):
        ds = FakeDataset({
            "genome_source": "igenomes",
            "genome": "GATK.GRCh38",
            "aligner": "parabricks",
        })
        preprocess.resolve_reference_genome(ds)
        self.assertNotIn("genome_source", ds.params)
        self.assertNotIn("bwa_index", ds.params)
        self.assertNotIn("bwamem2_index", ds.params)
        self.assertEqual(ds.params.get("genome"), "GATK.GRCh38")

    def test_custom_genome_uses_matching_bwa_index(self):
        ds = FakeDataset({
            "genome_source": "dataset",
            "aligner": "bwa-mem",
            "bwa_index": "s3://bucket/bwa-index",
        })
        preprocess.resolve_reference_genome(ds)
        self.assertEqual(ds.params.get("aligner"), "bwa-mem")
        self.assertEqual(ds.params.get("bwa"), "s3://bucket/bwa-index")
        self.assertNotIn("bwamem2", ds.params)

    def test_custom_genome_uses_matching_bwamem2_index(self):
        ds = FakeDataset({
            "genome_source": "dataset",
            "aligner": "bwa-mem2",
            "bwamem2_index": "s3://bucket/bwamem2-index",
        })
        preprocess.resolve_reference_genome(ds)
        self.assertEqual(ds.params.get("bwamem2"), "s3://bucket/bwamem2-index")
        self.assertNotIn("bwa", ds.params)

    def test_custom_genome_missing_index_raises(self):
        ds = FakeDataset({"genome_source": "dataset", "aligner": "bwa-mem2"})
        with self.assertRaisesRegex(ValueError, "no matching"):
            preprocess.resolve_reference_genome(ds)


class SkipBaserecalibrationWithoutKnownSitesTests(unittest.TestCase):

    def test_custom_genome_without_known_sites_skips_baserecalibrator(self):
        ds = FakeDataset({})
        preprocess.skip_baserecalibration_without_known_sites(ds, is_custom_genome=True)
        self.assertEqual(ds.params.get("skip_tools"), "baserecalibrator")

    def test_custom_genome_with_dbsnp_does_not_skip(self):
        ds = FakeDataset({"dbsnp": "s3://bucket/dbsnp.vcf.gz"})
        preprocess.skip_baserecalibration_without_known_sites(ds, is_custom_genome=True)
        self.assertNotIn("skip_tools", ds.params)

    def test_custom_genome_with_known_indels_does_not_skip(self):
        ds = FakeDataset({"known_indels": "s3://bucket/known_indels.vcf.gz"})
        preprocess.skip_baserecalibration_without_known_sites(ds, is_custom_genome=True)
        self.assertNotIn("skip_tools", ds.params)

    def test_igenomes_is_a_no_op(self):
        ds = FakeDataset({})
        preprocess.skip_baserecalibration_without_known_sites(ds, is_custom_genome=False)
        self.assertNotIn("skip_tools", ds.params)

    def test_preserves_existing_skip_tools(self):
        ds = FakeDataset({"skip_tools": "fastqc"})
        preprocess.skip_baserecalibration_without_known_sites(ds, is_custom_genome=True)
        self.assertEqual(ds.params.get("skip_tools"), "fastqc,baserecalibrator")


class StageCollidingVcfParamsTests(unittest.TestCase):

    def _run(self, params):
        ds = FakeDataset(params)
        fake_client = FakeS3Client()
        with patch.object(preprocess.boto3, "client", return_value=fake_client):
            preprocess.stage_colliding_vcf_params(ds)
        return ds, fake_client

    def test_no_collision_when_basenames_differ(self):
        ds, client = self._run({
            "input": "s3://bucket/dataset/config/manifest.csv",
            "dbsnp": "s3://bucket/refs/dbsnp/dbsnp.vcf.gz",
            "known_indels": "s3://bucket/refs/indels/known_indels.vcf.gz",
        })
        self.assertEqual(ds.params["dbsnp"], "s3://bucket/refs/dbsnp/dbsnp.vcf.gz")
        self.assertEqual(ds.params["known_indels"], "s3://bucket/refs/indels/known_indels.vcf.gz")
        self.assertEqual(client.copies, [])

    def test_collision_stages_unique_copies(self):
        # Both resolve through the references library's shared germline_resource
        # reference type, so both land on a file literally named
        # germline_resource.vcf.gz -- exactly the scenario this exists to handle.
        ds, client = self._run({
            "input": "s3://bucket/dataset/config/manifest.csv",
            "dbsnp": "s3://bucket/refs/dog10k-af/germline_resource.vcf.gz",
            "known_indels": "s3://bucket/refs/dog10k-indels/germline_resource.vcf.gz",
        })
        self.assertEqual(
            ds.params["dbsnp"],
            "s3://bucket/dataset/config/dbsnp_germline_resource.vcf.gz",
        )
        self.assertEqual(
            ds.params["known_indels"],
            "s3://bucket/dataset/config/known_indels_germline_resource.vcf.gz",
        )
        self.assertEqual(len(client.copies), 2)

    def test_only_one_param_set_is_a_no_op(self):
        ds, client = self._run({
            "input": "s3://bucket/dataset/config/manifest.csv",
            "dbsnp": "s3://bucket/refs/dog10k-af/germline_resource.vcf.gz",
        })
        self.assertEqual(ds.params["dbsnp"], "s3://bucket/refs/dog10k-af/germline_resource.vcf.gz")
        self.assertEqual(client.copies, [])

    def test_neither_param_set_is_a_no_op(self):
        ds, client = self._run({"input": "s3://bucket/dataset/config/manifest.csv"})
        self.assertEqual(client.copies, [])
        self.assertNotIn("dbsnp", ds.params)
        self.assertNotIn("known_indels", ds.params)


if __name__ == "__main__":
    unittest.main()
