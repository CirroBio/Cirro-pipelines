"""Offline tests for the wf-somatic-variation preprocess script.

preprocess.py imports the Cirro SDK only when it runs, so this is a plain
`python -m unittest` that needs nothing else. FakeDataset stands in for the
slice of PreprocessDataset that the checks touch.
"""
import importlib.util
import unittest
from pathlib import Path


def _load_preprocess():
    path = Path(__file__).with_name("preprocess.py")
    spec = importlib.util.spec_from_file_location("wf_somatic_variation_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()

TUMOR = "s3://bucket/datasets/bam-1/data/tumor.bam"
NORMAL = "s3://bucket/datasets/bam-2/data/normal.bam"
HG38 = "s3://bucket/datasets/ref-1/data/GRCh38.fa.gz"


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


class DropEmptyParams(unittest.TestCase):

    def test_unset_pickers_are_removed(self):
        ds = FakeDataset({"bam_tumor": TUMOR, "bam_normal": "", "pon_file": "", "snv": True})
        preprocess.drop_empty_params(ds)
        self.assertEqual(ds.params, {"bam_tumor": TUMOR, "snv": True})


class CheckSubworkflows(unittest.TestCase):

    def test_none_enabled_fails(self):
        with self.assertRaisesRegex(ValueError, "No analysis subworkflow selected"):
            preprocess.check_subworkflows(FakeDataset({"snv": False, "sv": False, "mod": False}))

    def test_enabled_are_returned(self):
        self.assertEqual(preprocess.check_subworkflows(FakeDataset({"snv": True, "mod": True})), ["snv", "mod"])


class CheckTumorOnly(unittest.TestCase):

    def test_no_normal_is_tumor_only(self):
        self.assertTrue(preprocess.check_tumor_only(FakeDataset({"bam_tumor": TUMOR}), ["snv"]))

    def test_picked_normal_is_paired(self):
        ds = FakeDataset({"bam_tumor": TUMOR, "bam_normal": NORMAL})
        self.assertFalse(preprocess.check_tumor_only(ds, ["snv"]))


class CheckReference(unittest.TestCase):

    def test_build_is_read_from_the_picked_path(self):
        self.assertEqual(preprocess.infer_reference_build(HG38), "hg38")

    def test_any_reference_passes(self):
        preprocess.check_reference(FakeDataset({"ref": HG38, "annotation": True}), ["sv"], True)


class CheckOptionCombinations(unittest.TestCase):

    def test_hybrid_and_genotyping_vcfs_together_fail(self):
        ds = FakeDataset({
            "hybrid_mode_vcf": "s3://bucket/data/a.vcf.gz",
            "genotyping_mode_vcf": "s3://bucket/data/b.vcf.gz",
        })
        with self.assertRaisesRegex(ValueError, "not both"):
            preprocess.check_option_combinations(ds, ["snv"])

    def test_one_vcf_passes(self):
        ds = FakeDataset({"hybrid_mode_vcf": "s3://bucket/data/a.vcf.gz"})
        preprocess.check_option_combinations(ds, ["snv"])


if __name__ == "__main__":
    unittest.main()
