"""Offline tests for the alignment selection in sarek_call_variants.

Runs without the Cirro SDK or boto3: both are stubbed before preprocess is loaded,
so this is a plain `python -m unittest`.

What these cover is the part that caused real harm on the platform -- variant
callers were being run on reads from preprocessing/mapped/, taken before
MarkDuplicates and ApplyBQSR, because the format (BAM vs CRAM) was decided before
the processing stage and duplicate-marked files ranked below raw mapped ones.
"""
import sys
import types
import unittest
from pathlib import Path

import pandas as pd


def _load_preprocess():
    helpers = types.ModuleType("cirro.helpers.preprocess_dataset")
    helpers.PreprocessDataset = object
    s3_path_mod = types.ModuleType("cirro.models.s3_path")

    class S3Path:
        def __init__(self, path):
            rest = path.replace("s3://", "", 1)
            self.bucket, _, self.key = rest.partition("/")

    s3_path_mod.S3Path = S3Path

    for name in ["cirro", "cirro.helpers", "cirro.models"]:
        sys.modules.setdefault(name, types.ModuleType(name))
    sys.modules["cirro.helpers.preprocess_dataset"] = helpers
    sys.modules["cirro.models.s3_path"] = s3_path_mod
    sys.modules.setdefault("boto3", types.ModuleType("boto3"))

    import importlib.util
    path = Path(__file__).with_name("preprocess.py")
    spec = importlib.util.spec_from_file_location("sarek_call_variants_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()

PREFIX = "s3://bucket/dataset/data"


class FakeDataset:
    """The slice of PreprocessDataset that select_alignments touches."""

    def __init__(self, files):
        self.files = pd.DataFrame(
            [dict(sample=sample, file=f"{PREFIX}/{path}") for sample, path in files]
        )


def sarek_outputs(sample, stage, suffix, fmt):
    data = f"preprocessing/{stage}/{sample}/{sample}.{suffix}.{fmt}"
    index = f"{data}.{'bai' if fmt == 'bam' else 'crai'}"
    return [(sample, data), (sample, index)]


class AlignmentStageTests(unittest.TestCase):

    def test_stage_from_sarek_path(self):
        for stage in ["recalibrated", "markduplicates", "sentieon_dedup", "mapped"]:
            path = f"{PREFIX}/preprocessing/{stage}/S1/S1.recal.bam"
            self.assertEqual(preprocess.alignment_stage(path), stage)

    def test_stage_from_filename_when_no_preprocessing_path(self):
        cases = {
            "S1.recal.bam": "recalibrated",
            "S1.md.bam": "markduplicates",
            "S1.dedup.bam": "sentieon_dedup",
            "S1.sorted.bam": "mapped",
            "S1.bam": "unknown",
        }
        for name, stage in cases.items():
            self.assertEqual(preprocess.alignment_stage(f"{PREFIX}/{name}"), stage)


class SelectAlignmentsTests(unittest.TestCase):

    def select(self, files, requested="best"):
        return preprocess.select_alignments(FakeDataset(files), requested)

    def test_recalibrated_wins_over_earlier_stages(self):
        files = (
            sarek_outputs("S1", "mapped", "sorted", "bam")
            + sarek_outputs("S1", "markduplicates", "md", "bam")
            + sarek_outputs("S1", "recalibrated", "recal", "bam")
        )
        selected, stage, fmt = self.select(files)
        self.assertEqual((stage, fmt), ("recalibrated", "bam"))
        self.assertEqual(list(selected["data"]), [f"{PREFIX}/preprocessing/recalibrated/S1/S1.recal.bam"])

    def test_duplicate_marked_cram_beats_mapped_bam(self):
        # The reported failure: format was chosen before stage, so any BAM in the
        # dataset forced the raw mapped reads to win over duplicate-marked CRAMs.
        files = (
            sarek_outputs("S1", "mapped", "sorted", "bam")
            + sarek_outputs("S1", "markduplicates", "md", "cram")
        )
        selected, stage, fmt = self.select(files)
        self.assertEqual((stage, fmt), ("markduplicates", "cram"))
        self.assertEqual(list(selected["index"]), [f"{PREFIX}/preprocessing/markduplicates/S1/S1.md.cram.crai"])

    def test_markduplicates_used_when_bqsr_skipped(self):
        files = (
            sarek_outputs("S1", "mapped", "sorted", "bam")
            + sarek_outputs("S1", "markduplicates", "md", "bam")
        )
        _, stage, fmt = self.select(files)
        self.assertEqual((stage, fmt), ("markduplicates", "bam"))

    def test_flat_intake_dataset(self):
        files = [("S1", "S1.bam"), ("S1", "S1.bam.bai"), ("S2", "S2.bam"), ("S2", "S2.bam.bai")]
        selected, stage, fmt = self.select(files)
        self.assertEqual((stage, fmt), ("unknown", "bam"))
        self.assertEqual(sorted(selected["sample"]), ["S1", "S2"])

    def test_alignment_without_index_is_excluded(self):
        files = (
            sarek_outputs("S1", "markduplicates", "md", "bam")
            + [("S1", "preprocessing/recalibrated/S1/S1.recal.bam")]
        )
        _, stage, _ = self.select(files)
        self.assertEqual(stage, "markduplicates")

    def test_no_indexed_alignments(self):
        files = [("S1", "preprocessing/recalibrated/S1/S1.recal.bam")]
        with self.assertRaisesRegex(ValueError, "No indexed BAM or CRAM"):
            self.select(files)

    def test_stage_not_shared_by_all_samples(self):
        files = (
            sarek_outputs("S1", "recalibrated", "recal", "bam")
            + sarek_outputs("S2", "markduplicates", "md", "cram")
        )
        with self.assertRaisesRegex(ValueError, "No single alignment stage and format"):
            self.select(files)

    def test_common_stage_chosen_over_a_stage_one_sample_lacks(self):
        files = (
            sarek_outputs("S1", "recalibrated", "recal", "bam")
            + sarek_outputs("S1", "markduplicates", "md", "bam")
            + sarek_outputs("S2", "markduplicates", "md", "bam")
        )
        selected, stage, _ = self.select(files)
        self.assertEqual(stage, "markduplicates")
        self.assertEqual(sorted(selected["sample"]), ["S1", "S2"])

    def test_requested_stage_overrides_best_available(self):
        files = (
            sarek_outputs("S1", "mapped", "sorted", "bam")
            + sarek_outputs("S1", "recalibrated", "recal", "bam")
        )
        selected, stage, _ = self.select(files, requested="mapped")
        self.assertEqual(stage, "mapped")
        self.assertEqual(list(selected["data"]), [f"{PREFIX}/preprocessing/mapped/S1/S1.sorted.bam"])

    def test_requested_stage_accepts_sentieon_dedup_as_duplicate_marked(self):
        files = sarek_outputs("S1", "sentieon_dedup", "dedup", "cram")
        _, stage, fmt = self.select(files, requested="markduplicates")
        self.assertEqual((stage, fmt), ("sentieon_dedup", "cram"))

    def test_requested_stage_absent_from_dataset(self):
        files = sarek_outputs("S1", "markduplicates", "md", "bam")
        with self.assertRaisesRegex(ValueError, "no recalibrated alignments"):
            self.select(files, requested="recalibrated")


if __name__ == "__main__":
    unittest.main()
