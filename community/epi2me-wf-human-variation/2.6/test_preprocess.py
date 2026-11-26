"""Offline tests for the wf-human-variation preprocess script.

preprocess.py imports the Cirro SDK only when it runs, so this is a plain
`python -m unittest` that needs nothing else. FakeDataset stands in for the
slice of PreprocessDataset that the checks touch.
"""
import importlib.util
import unittest
from pathlib import Path


def _load_preprocess():
    path = Path(__file__).with_name("preprocess.py")
    spec = importlib.util.spec_from_file_location("wf_human_variation_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()

BAM = "s3://bucket/datasets/bam-1/data/sample.bam"
HG38 = "s3://bucket/datasets/ref-1/data/GRCh38.fa.gz"
HG19 = "s3://bucket/datasets/ref-2/data/hg19.fasta"


class FakeLogger:
    def info(self, *args):
        pass

    def warning(self, *args):
        pass


class FakeDataset:

    def __init__(self, params):
        self.params = dict(params)
        self.logger = FakeLogger()

    def remove_param(self, name, force=False):
        self.params.pop(name, None)


class DropParams(unittest.TestCase):

    def test_unset_pickers_are_removed(self):
        ds = FakeDataset({"bam": BAM, "bed": "", "coverage_bed": "", "snp": True})
        preprocess.drop_empty_params(ds)
        self.assertEqual(ds.params, {"bam": BAM, "snp": True})

    def test_dependents_of_disabled_toggles_are_removed(self):
        ds = FakeDataset({
            "downsample_coverage": False, "downsample_coverage_target": 60,
            "use_qdnaseq": True, "qdnaseq_bin_size": 500,
        })
        preprocess.drop_inactive_dependents(ds)
        self.assertNotIn("downsample_coverage_target", ds.params)
        self.assertEqual(ds.params["qdnaseq_bin_size"], 500)


class CheckSubworkflows(unittest.TestCase):

    def test_none_enabled_fails(self):
        with self.assertRaisesRegex(ValueError, "No analysis subworkflow selected"):
            preprocess.check_subworkflows(FakeDataset({"snp": False, "sv": False}))

    def test_enabled_are_returned_in_order(self):
        ds = FakeDataset({"snp": True, "sv": True, "cnv": False})
        self.assertEqual(preprocess.check_subworkflows(ds), ["sv", "snp"])


class CheckReference(unittest.TestCase):

    def test_build_is_read_from_the_picked_path(self):
        self.assertEqual(preprocess.infer_reference_build(HG38), "hg38")
        self.assertEqual(preprocess.infer_reference_build(HG19), "hg19")
        self.assertIsNone(preprocess.infer_reference_build("s3://bucket/data/genome.fa"))

    def test_str_on_hg19_fails(self):
        with self.assertRaisesRegex(ValueError, "requires an hg38 reference"):
            preprocess.check_reference(FakeDataset({"ref": HG19}), ["str"])

    def test_str_on_hg38_passes(self):
        preprocess.check_reference(FakeDataset({"ref": HG38}), ["str"])


class CheckOptionCombinations(unittest.TestCase):

    def test_phasing_without_snp_fails(self):
        with self.assertRaisesRegex(ValueError, "Phasing requires the SNP subworkflow"):
            preprocess.check_option_combinations(FakeDataset({"phased": True}), ["sv"])

    def test_custom_basecaller_model_fails(self):
        with self.assertRaisesRegex(ValueError, "custom Clair3 model"):
            preprocess.check_option_combinations(
                FakeDataset({"override_basecaller_cfg": "custom"}), ["snp"]
            )


if __name__ == "__main__":
    unittest.main()
