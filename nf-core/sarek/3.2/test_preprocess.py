"""Offline tests for the nf-core/sarek samplesheet handling.

preprocess.py imports the Cirro SDK only when it runs, so this is a plain
`python -m unittest` that needs pandas and nothing else.

The samplesheet is the structured input Cirro writes from the run form: a CSV with
CRLF line endings, a header from the form's columns, and blank cells left empty.
"""
import importlib.util
import io
import unittest
from pathlib import Path

import pandas as pd


def _load_preprocess():
    path = Path(__file__).with_name("preprocess.py")
    spec = importlib.util.spec_from_file_location("sarek_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()

HEADER = "patient,sex,status,sample,lane,fastq_1,fastq_2\r\n"

CONFIG = "s3://bucket/datasets/ds1/config"
GERMLINE_RESOURCE = "s3://refs/germline_resource/gnomad/germline_resource.vcf.gz"


def samplesheet(rows):
    """A table of (patient, sex, status, sample, lane) rows, with placeholder FASTQ paths."""
    table = pd.DataFrame(rows, columns=["patient", "sex", "status", "sample", "lane"])
    return table.assign(fastq_1="r1.fastq.gz", fastq_2="r2.fastq.gz")


class ReadSamplesheet(unittest.TestCase):

    def test_values_stay_text(self):
        table = preprocess.read_samplesheet(io.StringIO(
            HEADER + ",,,01,001,s3://bucket/data/01_R1.fastq.gz,s3://bucket/data/01_R2.fastq.gz\r\n"
        ))
        self.assertEqual(table["sample"].tolist(), ["01"])
        self.assertEqual(table["lane"].tolist(), ["001"])
        self.assertEqual(table["patient"].tolist(), [""])
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

    def test_sex_is_normalized(self):
        filled = preprocess.fill_samplesheet(samplesheet([
            ("", "male", "", "S1", "1"),
            ("", "", "", "S2", "1"),
        ]))
        self.assertEqual(filled["sex"].tolist(), ["XY", "NA"])

    def test_status_is_normalized_and_blank_is_normal(self):
        filled = preprocess.fill_samplesheet(samplesheet([
            ("P1", "", "Normal", "N1", "1"),
            ("P1", "", "Tumor", "T1", "1"),
            ("P2", "", "", "N2", "1"),
            ("P2", "", "1.0", "T2", "1"),
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

    def test_columns_keep_their_order(self):
        table = samplesheet([("P1", "XX", "1", "S1", "1")])
        filled = preprocess.fill_samplesheet(table)
        self.assertEqual(filled.columns.tolist(), table.columns.tolist())
        self.assertTrue(filled.equals(table))


class PatientNamedLikeItsSample(unittest.TestCase):

    def test_passes_every_check(self):
        filled = preprocess.fill_samplesheet(samplesheet([
            ("P1", "XX", "0", "P1", "1"),
            ("P1", "XX", "1", "P1-T", "1"),
            ("", "XY", "0", "P2", "1"),
        ]))
        self.assertEqual(filled["patient"].tolist(), ["P1", "P1", "P2"])
        preprocess.check_one_patient_per_sample(filled)
        preprocess.check_one_status_per_sample(filled)
        preprocess.check_unique_lanes(filled)
        preprocess.check_one_normal_per_tumor_patient(filled)
        preprocess.check_status_requirements(filled, ["mutect2", "haplotypecaller"])
        preprocess.check_sex_requirements(filled, ["ascat"])


class CheckOnePatientPerSample(unittest.TestCase):

    def test_consistent_patients_pass(self):
        preprocess.check_one_patient_per_sample(samplesheet([
            ("P1", "NA", "0", "S1", "1"),
            ("P1", "NA", "0", "S1", "2"),
            ("P1", "NA", "1", "S2", "1"),
        ]))

    def test_sample_with_two_patients_fails(self):
        table = samplesheet([
            ("P1", "NA", "0", "S1", "1"),
            ("S1", "NA", "0", "S1", "2"),
        ])
        with self.assertRaisesRegex(ValueError, r"these have several: S1 \(P1, S1\)"):
            preprocess.check_one_patient_per_sample(table)


class CheckOneStatusPerSample(unittest.TestCase):

    def test_consistent_statuses_pass(self):
        preprocess.check_one_status_per_sample(samplesheet([
            ("P1", "NA", "1", "T1", "1"),
            ("P1", "NA", "1", "T1", "2"),
            ("P1", "NA", "0", "N1", "1"),
        ]))

    def test_sample_with_both_statuses_fails(self):
        table = samplesheet([
            ("P1", "NA", "1", "T1", "1"),
            ("P1", "NA", "0", "T1", "2"),
        ])
        with self.assertRaisesRegex(ValueError, "rows of both: T1"):
            preprocess.check_one_status_per_sample(table)


class CheckUniqueLanes(unittest.TestCase):

    def test_distinct_lanes_pass(self):
        preprocess.check_unique_lanes(samplesheet([
            ("P1", "NA", "0", "S1", "1"),
            ("P1", "NA", "0", "S1", "2"),
            ("P1", "NA", "1", "S2", "1"),
        ]))

    def test_repeated_lane_fails(self):
        table = samplesheet([
            ("P1", "NA", "0", "S1", "1"),
            ("P1", "NA", "0", "S1", "1"),
        ])
        with self.assertRaisesRegex(ValueError, "patient P1, sample S1, lane 1"):
            preprocess.check_unique_lanes(table)


class CheckOneNormalPerTumorPatient(unittest.TestCase):

    def test_one_normal_with_several_tumors_passes(self):
        preprocess.check_one_normal_per_tumor_patient(samplesheet([
            ("P1", "NA", "0", "N1", "1"),
            ("P1", "NA", "0", "N1", "2"),
            ("P1", "NA", "1", "T1", "1"),
            ("P1", "NA", "1", "T1-relapse", "1"),
        ]))

    def test_several_normals_without_a_tumor_pass(self):
        preprocess.check_one_normal_per_tumor_patient(samplesheet([
            ("P1", "NA", "0", "N1", "1"),
            ("P1", "NA", "0", "N2", "1"),
        ]))

    def test_tumor_patient_with_two_normals_fails(self):
        table = samplesheet([
            ("P1", "NA", "0", "N1", "1"),
            ("P1", "NA", "0", "N2", "1"),
            ("P1", "NA", "1", "T1", "1"),
        ])
        with self.assertRaisesRegex(ValueError, r"these patients have several: P1 \(N1, N2\)"):
            preprocess.check_one_normal_per_tumor_patient(table)


class CheckStatusRequirements(unittest.TestCase):

    def test_somatic_callers_without_a_tumor_fail(self):
        table = samplesheet([("P1", "XX", "0", "N1", "1")])
        with self.assertRaisesRegex(
            ValueError, r"ascat, controlfreec, msisensorpro, mutect2 needs at least one tumor sample \(status 1\)"
        ):
            preprocess.check_status_requirements(table, ["strelka", "mutect2", "msisensorpro", "controlfreec", "ascat"])

    def test_callers_needing_a_normal_fail_without_one(self):
        table = samplesheet([("P1", "XX", "1", "T1", "1")])
        with self.assertRaisesRegex(
            ValueError, r"ascat, deepvariant, haplotypecaller, msisensorpro needs at least one normal sample \(status 0\)"
        ):
            preprocess.check_status_requirements(
                table, ["haplotypecaller", "mutect2", "deepvariant", "ascat", "msisensorpro"]
            )

    def test_tools_without_requirements_pass(self):
        for status in ("0", "1"):
            preprocess.check_status_requirements(
                samplesheet([("P1", "NA", status, "S1", "1")]), ["freebayes", "manta", "strelka", "tiddit"]
            )

    def test_paired_samplesheet_passes_every_status_requirement(self):
        table = samplesheet([
            ("P1", "XX", "0", "N1", "1"),
            ("P1", "XX", "1", "T1", "1"),
        ])
        preprocess.check_status_requirements(
            table, ["ascat", "controlfreec", "deepvariant", "haplotypecaller", "msisensorpro", "mutect2"]
        )


class CheckSexRequirements(unittest.TestCase):

    def test_cnv_callers_need_every_sex(self):
        table = samplesheet([
            ("P1", "XX", "0", "N1", "1"),
            ("P1", "NA", "1", "T1", "1"),
            ("P1", "NA", "1", "T1", "2"),
        ])
        with self.assertRaisesRegex(ValueError, r"ascat, controlfreec needs the sex .* these have none: T1\."):
            preprocess.check_sex_requirements(table, ["controlfreec", "ascat", "mutect2"])

    def test_known_sexes_pass(self):
        table = samplesheet([
            ("P1", "XX", "0", "N1", "1"),
            ("P1", "XX", "1", "T1", "1"),
        ])
        preprocess.check_sex_requirements(table, ["ascat", "controlfreec"])

    def test_other_tools_ignore_sex(self):
        preprocess.check_sex_requirements(
            samplesheet([("P1", "NA", "1", "T1", "1")]), ["mutect2", "strelka", "manta"]
        )


class CollidingVcfCopies(unittest.TestCase):

    INPUT = f"{CONFIG}/pipeline_samplesheet.csv"

    def test_distinct_names_need_no_copy(self):
        params = {
            "input": self.INPUT,
            "germline_resource": GERMLINE_RESOURCE,
            "pon": "s3://refs/pon/pon.vcf.gz",
        }
        self.assertEqual(preprocess.colliding_vcf_copies(params), [])

    def test_colliding_names_are_copied_beside_the_samplesheet(self):
        pon = "s3://refs/germline_resource/my-pon/germline_resource.vcf.gz"
        params = {"input": self.INPUT, "germline_resource": GERMLINE_RESOURCE, "pon": pon}
        self.assertEqual(preprocess.colliding_vcf_copies(params), [
            ("germline_resource", GERMLINE_RESOURCE, f"{CONFIG}/germline_resource_germline_resource.vcf.gz"),
            ("pon", pon, f"{CONFIG}/pon_germline_resource.vcf.gz"),
        ])

    def test_an_index_is_copied_beside_its_vcf(self):
        pon = "s3://refs/germline_resource/my-pon/germline_resource.vcf.gz"
        params = {
            "input": self.INPUT,
            "germline_resource": GERMLINE_RESOURCE,
            "germline_resource_tbi": f"{GERMLINE_RESOURCE}.tbi",
            "pon": pon,
        }
        self.assertEqual(preprocess.colliding_vcf_copies(params), [
            ("germline_resource", GERMLINE_RESOURCE, f"{CONFIG}/germline_resource_germline_resource.vcf.gz"),
            (
                "germline_resource_tbi",
                f"{GERMLINE_RESOURCE}.tbi",
                f"{CONFIG}/germline_resource_germline_resource.vcf.gz.tbi",
            ),
            ("pon", pon, f"{CONFIG}/pon_germline_resource.vcf.gz"),
        ])

    def test_a_single_resource_needs_no_copy(self):
        self.assertEqual(preprocess.colliding_vcf_copies({"input": self.INPUT, "pon": GERMLINE_RESOURCE}), [])
        self.assertEqual(preprocess.colliding_vcf_copies({"input": self.INPUT}), [])


if __name__ == "__main__":
    unittest.main()
