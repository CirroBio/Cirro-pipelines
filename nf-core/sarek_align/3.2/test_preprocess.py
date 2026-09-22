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


def _load_preprocess():
    helpers = types.ModuleType("cirro.helpers.preprocess_dataset")
    helpers.PreprocessDataset = object

    for name in ["cirro", "cirro.helpers"]:
        sys.modules.setdefault(name, types.ModuleType(name))
    sys.modules["cirro.helpers.preprocess_dataset"] = helpers
    sys.modules.setdefault("pandas", types.ModuleType("pandas"))

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


if __name__ == "__main__":
    unittest.main()
