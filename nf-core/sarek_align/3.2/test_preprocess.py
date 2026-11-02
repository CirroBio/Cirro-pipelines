"""Offline tests for the samplesheet and reference-genome/aligner handling in sarek_align.

preprocess.py imports the Cirro SDK only when it runs, so this is a plain
`python -m unittest` that needs pandas and nothing else.

The samplesheet is the structured input Cirro writes from the run form: a CSV with
CRLF line endings, a header from the form's columns, and blank cells left empty.

Regression coverage for: the iGenomes branch of genome_source dropped the
`aligner` param entirely (it was only restored on the Custom Genome branch),
so selecting bwa-mem2 or the GPU-accelerated Parabricks aligner in the form
silently had no effect -- sarek always fell back to its own default (bwa-mem).
"""
import importlib.util
import io
import unittest
from pathlib import Path

import pandas as pd


def _load_preprocess():
    path = Path(__file__).with_name("preprocess.py")
    spec = importlib.util.spec_from_file_location("sarek_align_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()

HEADER = "patient,sex,status,sample,lane,fastq_1,fastq_2\r\n"

CONFIG = "s3://bucket/dataset/config"
INPUT = f"{CONFIG}/pipeline_samplesheet.csv"


def samplesheet(rows):
    """A table of (patient, sex, status, sample, lane) rows, with placeholder FASTQ paths."""
    table = pd.DataFrame(rows, columns=["patient", "sex", "status", "sample", "lane"])
    return table.assign(fastq_1="r1.fastq.gz", fastq_2="r2.fastq.gz")


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


class CollidingVcfCopiesTests(unittest.TestCase):

    def _run(self, params):
        copies = preprocess.colliding_vcf_copies(params)
        return params, copies

    def test_no_collision_when_basenames_differ(self):
        params, copies = self._run({
            "input": INPUT,
            "dbsnp": "s3://bucket/refs/dbsnp/dbsnp.vcf.gz",
            "known_indels": "s3://bucket/refs/indels/known_indels.vcf.gz",
        })
        self.assertEqual(params["dbsnp"], "s3://bucket/refs/dbsnp/dbsnp.vcf.gz")
        self.assertEqual(params["known_indels"], "s3://bucket/refs/indels/known_indels.vcf.gz")
        self.assertEqual(copies, [])

    def test_collision_stages_unique_copies(self):
        params, copies = self._run({
            "input": INPUT,
            "dbsnp": "s3://bucket/refs/dog10k-af/germline_resource.vcf.gz",
            "known_indels": "s3://bucket/refs/dog10k-indels/germline_resource.vcf.gz",
        })
        self.assertEqual(copies, [
            (
                "dbsnp",
                "s3://bucket/refs/dog10k-af/germline_resource.vcf.gz",
                "s3://bucket/dataset/config/dbsnp_germline_resource.vcf.gz",
            ),
            (
                "known_indels",
                "s3://bucket/refs/dog10k-indels/germline_resource.vcf.gz",
                "s3://bucket/dataset/config/known_indels_germline_resource.vcf.gz",
            ),
        ])
        self.assertEqual(len(copies), 2)

    def test_collision_copies_each_index_beside_its_vcf(self):
        params, copies = self._run({
            "input": INPUT,
            "dbsnp": "s3://bucket/refs/dog10k-af/germline_resource.vcf.gz",
            "dbsnp_tbi": "s3://bucket/refs/dog10k-af/germline_resource.vcf.gz.tbi",
            "known_indels": "s3://bucket/refs/dog10k-indels/germline_resource.vcf.gz",
        })
        self.assertEqual(copies[:2], [
            (
                "dbsnp",
                "s3://bucket/refs/dog10k-af/germline_resource.vcf.gz",
                "s3://bucket/dataset/config/dbsnp_germline_resource.vcf.gz",
            ),
            (
                "dbsnp_tbi",
                "s3://bucket/refs/dog10k-af/germline_resource.vcf.gz.tbi",
                "s3://bucket/dataset/config/dbsnp_germline_resource.vcf.gz.tbi",
            ),
        ])
        self.assertEqual(len(copies), 3)

    def test_only_one_param_set_is_a_no_op(self):
        params, copies = self._run({
            "input": INPUT,
            "dbsnp": "s3://bucket/refs/dog10k-af/germline_resource.vcf.gz",
        })
        self.assertEqual(params["dbsnp"], "s3://bucket/refs/dog10k-af/germline_resource.vcf.gz")
        self.assertEqual(copies, [])

    def test_neither_param_set_is_a_no_op(self):
        params, copies = self._run({"input": INPUT})
        self.assertEqual(copies, [])
        self.assertNotIn("dbsnp", params)
        self.assertNotIn("known_indels", params)


class ReadSamplesheet(unittest.TestCase):

    def test_values_stay_text(self):
        table = preprocess.read_samplesheet(io.StringIO(
            HEADER + ",,,01,001,s3://bucket/data/01_R1.fastq.gz,s3://bucket/data/01_R2.fastq.gz\r\n"
        ))
        self.assertEqual(table["sample"].tolist(), ["01"])
        self.assertEqual(table["lane"].tolist(), ["001"])
        self.assertEqual(table["patient"].tolist(), [""])
        self.assertEqual(table["sex"].tolist(), [""])
        self.assertEqual(table["status"].tolist(), [""])

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


class NormalizeSex(unittest.TestCase):

    def test_spellings_map_to_sarek_values(self):
        cases = {
            "male": "XY", " M ": "XY", "XY": "XY",
            "Female": "XX", "f": "XX", "xx": "XX",
            "": "NA", "NA": "NA", "unknown": "NA", "other": "NA",
        }
        for value, expected in cases.items():
            self.assertEqual(preprocess.normalize_sex(value), expected, value)


class NormalizeStatus(unittest.TestCase):

    def test_spellings_map_to_sarek_values(self):
        cases = {
            "": "0", " ": "0", "0": "0", "0.0": "0", "Normal": "0", "normal": "0",
            "1": "1", "1.0": "1", "Tumor": "1", " TUMOR ": "1", "tumour": "1",
        }
        for value, expected in cases.items():
            self.assertEqual(preprocess.normalize_status(value), expected, value)

    def test_unknown_value_fails_clearly(self):
        for value in ("2", "relapse", "yes"):
            with self.assertRaisesRegex(ValueError, f"not '{value}'"):
                preprocess.normalize_status(value)


class FillSamplesheet(unittest.TestCase):

    def test_blank_patient_takes_the_sample_name(self):
        filled = preprocess.fill_samplesheet(samplesheet([
            ("", "", "", "S1", "1"),
            ("P2", "", "", "S2", "1"),
        ]))
        self.assertEqual(filled["patient"].tolist(), ["S1", "P2"])

    def test_blank_sex_stays_unknown(self):
        filled = preprocess.fill_samplesheet(samplesheet([
            ("", "", "", "S1", "1"),
            ("", "female", "", "S2", "1"),
        ]))
        self.assertEqual(filled["sex"].tolist(), ["NA", "XX"])

    def test_tumor_and_normal_spellings_become_status_codes(self):
        filled = preprocess.fill_samplesheet(samplesheet([
            ("P1", "", "Normal", "N1", "1"),
            ("P1", "", "Tumor", "T1", "1"),
            ("P2", "", "", "N2", "1"),
            ("P2", "", "1", "T2", "1"),
        ]))
        self.assertEqual(filled["status"].tolist(), ["0", "1", "0", "1"])

    def test_blank_lanes_are_numbered_within_each_sample(self):
        filled = preprocess.fill_samplesheet(samplesheet([
            ("", "", "", "S1", ""),
            ("", "", "", "S1", ""),
            ("", "", "", "S2", ""),
        ]))
        self.assertEqual(filled["lane"].tolist(), ["1", "2", "1"])

    def test_blank_lanes_skip_lanes_the_sample_already_uses(self):
        filled = preprocess.fill_samplesheet(samplesheet([
            ("", "", "", "S1", "1"),
            ("", "", "", "S1", ""),
            ("", "", "", "S1", "L002"),
            ("", "", "", "S1", ""),
        ]))
        self.assertEqual(filled["lane"].tolist(), ["1", "2", "L002", "3"])

    def test_complete_table_is_unchanged(self):
        table = samplesheet([("P1", "XX", "1", "S1", "1")])
        filled = preprocess.fill_samplesheet(table)
        self.assertEqual(filled.columns.tolist(), table.columns.tolist())
        self.assertTrue(filled.equals(table))


class SamplesheetChecks(unittest.TestCase):

    def test_patient_named_like_its_sample_passes(self):
        filled = preprocess.fill_samplesheet(samplesheet([
            ("P1", "", "0", "P1", "1"),
            ("P1", "", "1", "P1-T", "1"),
            ("", "", "", "P2", ""),
            ("", "", "", "P2", ""),
        ]))
        preprocess.check_one_patient_per_sample(filled)
        preprocess.check_one_status_per_sample(filled)
        preprocess.check_unique_lanes(filled)
        preprocess.check_one_normal_per_tumor_patient(filled)

    def test_sample_with_two_patients_fails(self):
        table = samplesheet([
            ("P1", "NA", "0", "S1", "1"),
            ("S1", "NA", "0", "S1", "2"),
        ])
        with self.assertRaisesRegex(ValueError, r"these have several: S1 \(P1, S1\)"):
            preprocess.check_one_patient_per_sample(table)

    def test_sample_with_both_statuses_fails(self):
        table = samplesheet([
            ("P1", "NA", "1", "T1", "1"),
            ("P1", "NA", "0", "T1", "2"),
        ])
        with self.assertRaisesRegex(ValueError, "rows of both: T1"):
            preprocess.check_one_status_per_sample(table)

    def test_repeated_lane_fails(self):
        table = samplesheet([
            ("P1", "NA", "0", "S1", "1"),
            ("P1", "NA", "0", "S1", "1"),
        ])
        with self.assertRaisesRegex(ValueError, "patient P1, sample S1, lane 1"):
            preprocess.check_unique_lanes(table)

    def test_tumor_patient_with_two_normals_fails(self):
        table = samplesheet([
            ("P1", "NA", "0", "N1", "1"),
            ("P1", "NA", "0", "N2", "1"),
            ("P1", "NA", "1", "T1", "1"),
        ])
        with self.assertRaisesRegex(ValueError, r"these patients have several: P1 \(N1, N2\)"):
            preprocess.check_one_normal_per_tumor_patient(table)


if __name__ == "__main__":
    unittest.main()
