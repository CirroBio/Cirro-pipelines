"""Offline tests for the nf-core/sarek variant-calling samplesheet handling.

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
    spec = importlib.util.spec_from_file_location("sarek_call_variants_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()

COLUMNS = ["patient", "sex", "status", "sample", "bam", "bai", "cram", "crai"]
HEADER = ",".join(COLUMNS) + "\r\n"

DATA = "s3://bucket/datasets/ds1/data"
CONFIG = "s3://bucket/datasets/ds2/config"
GERMLINE_RESOURCE = "s3://refs/germline_resource/gnomad/germline_resource.vcf.gz"


def samplesheet(rows):
    """A table of (patient, sex, status, sample) rows, each with a recalibrated BAM and its index."""
    records = []
    for patient, sex, status, sample in rows:
        bam = f"{DATA}/preprocessing/recalibrated/{sample}/{sample}.recal.bam"
        records.append(dict(
            patient=patient, sex=sex, status=status, sample=sample,
            bam=bam, bai=f"{bam}.bai", cram="", crai="",
        ))
    return pd.DataFrame(records, columns=COLUMNS)


def alignment(**cells):
    """The alignment cells of a row, blank except for those given."""
    return {"sample": "S1", "bam": "", "bai": "", "cram": "", "crai": "", **cells}


class ReadSamplesheet(unittest.TestCase):

    def test_values_stay_text(self):
        table = preprocess.read_samplesheet(io.StringIO(
            HEADER + ",,,01,s3://bucket/data/01.bam,s3://bucket/data/01.bam.bai,,\r\n"
        ))
        self.assertEqual(table["sample"].tolist(), ["01"])
        self.assertEqual(table["patient"].tolist(), [""])
        self.assertEqual(table["status"].tolist(), [""])
        self.assertEqual(table["cram"].tolist(), [""])

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
            ("", "", "", "S1"),
            ("P2", "", "", "S2"),
        ]))
        self.assertEqual(filled["patient"].tolist(), ["S1", "P2"])

    def test_sex_is_normalized(self):
        filled = preprocess.fill_samplesheet(samplesheet([
            ("", "male", "", "S1"),
            ("", "", "", "S2"),
        ]))
        self.assertEqual(filled["sex"].tolist(), ["XY", "NA"])

    def test_status_is_normalized_and_blank_is_normal(self):
        filled = preprocess.fill_samplesheet(samplesheet([
            ("P1", "", "Normal", "N1"),
            ("P1", "", "Tumor", "T1"),
            ("P2", "", "", "N2"),
            ("P2", "", "1.0", "T2"),
        ]))
        self.assertEqual(filled["status"].tolist(), ["0", "1", "0", "1"])

    def test_a_complete_table_is_unchanged(self):
        # sarek's own csv/recalibrated.csv, which the form fills from, needs no changes
        table = samplesheet([("P1", "XY", "0", "N1"), ("P1", "XY", "1", "T1")])
        filled = preprocess.fill_samplesheet(table)
        self.assertEqual(filled.columns.tolist(), COLUMNS)
        self.assertTrue(filled.equals(table))


class AlignmentProblem(unittest.TestCase):

    def test_bam_with_its_index_passes(self):
        for bai in ("S1.recal.bam.bai", "S1.recal.bai"):
            row = alignment(bam=f"{DATA}/S1.recal.bam", bai=f"{DATA}/{bai}")
            self.assertIsNone(preprocess.alignment_problem(row), bai)

    def test_cram_with_its_index_passes(self):
        for crai in ("S1.recal.cram.crai", "S1.recal.crai"):
            row = alignment(cram=f"{DATA}/S1.recal.cram", crai=f"{DATA}/{crai}")
            self.assertIsNone(preprocess.alignment_problem(row), crai)

    def test_index_in_another_folder_passes(self):
        row = alignment(bam=f"{DATA}/bams/S1.bam", bai=f"{DATA}/indexes/S1.bam.bai")
        self.assertIsNone(preprocess.alignment_problem(row))

    def test_both_formats_fail(self):
        row = alignment(bam="S1.bam", bai="S1.bam.bai", cram="S1.cram", crai="S1.cram.crai")
        self.assertRegex(preprocess.alignment_problem(row), "not both or neither")

    def test_no_alignment_fails(self):
        self.assertRegex(preprocess.alignment_problem(alignment()), "not both or neither")

    def test_alignment_without_its_index_fails(self):
        self.assertEqual(preprocess.alignment_problem(alignment(bam="S1.bam")), "the BAI is missing")
        self.assertEqual(preprocess.alignment_problem(alignment(crai="S1.cram.crai")), "the CRAM is missing")

    def test_index_of_another_alignment_fails(self):
        row = alignment(bam=f"{DATA}/S1.bam", bai=f"{DATA}/S2.bam.bai")
        self.assertEqual(
            preprocess.alignment_problem(row),
            "S2.bam.bai is not the index of S1.bam; expected S1.bam.bai or S1.bai",
        )

    def test_index_name_from_markduplicates_no_table_csv_fails(self):
        # sarek writes BAM indexes into csv/markduplicates_no_table.csv with .bam doubled
        row = alignment(bam=f"{DATA}/S1.md.bam", bai=f"{DATA}/S1.md.bam.bam.bai")
        self.assertRegex(preprocess.alignment_problem(row), "expected S1.md.bam.bai or S1.md.bai")


class CheckOneAlignmentPerRow(unittest.TestCase):

    def test_bam_and_cram_rows_together_pass(self):
        table = samplesheet([("P1", "XX", "0", "N1"), ("P1", "XX", "1", "T1")])
        table.loc[1, ["bam", "bai"]] = ""
        table.loc[1, ["cram", "crai"]] = [f"{DATA}/T1.cram", f"{DATA}/T1.cram.crai"]
        preprocess.check_one_alignment_per_row(table)

    def test_every_bad_row_is_listed(self):
        table = samplesheet([("P1", "XX", "0", "N1"), ("P1", "XX", "1", "T1"), ("P2", "XX", "0", "N2")])
        table.loc[0, "bai"] = ""
        table.loc[2, "cram"] = f"{DATA}/N2.cram"
        with self.assertRaisesRegex(
            ValueError, "N1: the BAI is missing; N2: give either a BAM and its BAI or a CRAM and its CRAI"
        ):
            preprocess.check_one_alignment_per_row(table)


class CheckOnePatientPerSample(unittest.TestCase):

    def test_distinct_samples_pass(self):
        preprocess.check_one_patient_per_sample(samplesheet([
            ("P1", "NA", "0", "N1"),
            ("P1", "NA", "1", "T1"),
        ]))

    def test_sample_with_two_patients_fails(self):
        table = samplesheet([
            ("P1", "NA", "0", "S1"),
            ("P2", "NA", "0", "S1"),
        ])
        with self.assertRaisesRegex(ValueError, r"these have several: S1 \(P1, P2\)"):
            preprocess.check_one_patient_per_sample(table)


class CheckOneRowPerSample(unittest.TestCase):

    def test_distinct_samples_pass(self):
        preprocess.check_one_row_per_sample(samplesheet([
            ("P1", "NA", "0", "N1"),
            ("P1", "NA", "1", "T1"),
        ]))

    def test_sample_on_two_rows_fails(self):
        table = samplesheet([
            ("P1", "NA", "0", "N1"),
            ("P1", "NA", "0", "N1"),
            ("P1", "NA", "1", "T1"),
        ])
        with self.assertRaisesRegex(ValueError, "on more than one: N1\\."):
            preprocess.check_one_row_per_sample(table)


class CheckOneNormalPerTumorPatient(unittest.TestCase):

    def test_one_normal_with_several_tumors_passes(self):
        preprocess.check_one_normal_per_tumor_patient(samplesheet([
            ("P1", "NA", "0", "N1"),
            ("P1", "NA", "1", "T1"),
            ("P1", "NA", "1", "T1-relapse"),
        ]))

    def test_several_normals_without_a_tumor_pass(self):
        preprocess.check_one_normal_per_tumor_patient(samplesheet([
            ("P1", "NA", "0", "N1"),
            ("P1", "NA", "0", "N2"),
        ]))

    def test_tumor_patient_with_two_normals_fails(self):
        table = samplesheet([
            ("P1", "NA", "0", "N1"),
            ("P1", "NA", "0", "N2"),
            ("P1", "NA", "1", "T1"),
        ])
        with self.assertRaisesRegex(ValueError, r"these patients have several: P1 \(N1, N2\)"):
            preprocess.check_one_normal_per_tumor_patient(table)


class CheckStatusRequirements(unittest.TestCase):

    def test_somatic_callers_without_a_tumor_fail(self):
        table = samplesheet([("P1", "XX", "0", "N1")])
        with self.assertRaisesRegex(
            ValueError, r"ascat, controlfreec, msisensorpro, mutect2 needs at least one tumor sample \(status 1\)"
        ):
            preprocess.check_status_requirements(table, ["strelka", "mutect2", "msisensorpro", "controlfreec", "ascat"])

    def test_callers_needing_a_normal_fail_without_one(self):
        table = samplesheet([("P1", "XX", "1", "T1")])
        with self.assertRaisesRegex(
            ValueError, r"ascat, deepvariant, haplotypecaller, msisensorpro needs at least one normal sample \(status 0\)"
        ):
            preprocess.check_status_requirements(
                table, ["haplotypecaller", "mutect2", "deepvariant", "ascat", "msisensorpro"]
            )

    def test_tools_without_requirements_pass(self):
        for status in ("0", "1"):
            preprocess.check_status_requirements(
                samplesheet([("P1", "NA", status, "S1")]), ["freebayes", "manta", "strelka", "tiddit"]
            )

    def test_paired_samplesheet_passes_every_status_requirement(self):
        table = samplesheet([("P1", "XX", "0", "N1"), ("P1", "XX", "1", "T1")])
        preprocess.check_status_requirements(
            table, ["ascat", "controlfreec", "deepvariant", "haplotypecaller", "msisensorpro", "mutect2"]
        )


class CheckSexRequirements(unittest.TestCase):

    def test_cnv_callers_need_every_sex(self):
        table = samplesheet([("P1", "XX", "0", "N1"), ("P1", "NA", "1", "T1")])
        with self.assertRaisesRegex(ValueError, r"ascat, controlfreec needs the sex .* these have none: T1\."):
            preprocess.check_sex_requirements(table, ["controlfreec", "ascat", "mutect2"])

    def test_known_sexes_pass(self):
        table = samplesheet([("P1", "XX", "0", "N1"), ("P1", "XX", "1", "T1")])
        preprocess.check_sex_requirements(table, ["ascat", "controlfreec"])

    def test_other_tools_ignore_sex(self):
        preprocess.check_sex_requirements(
            samplesheet([("P1", "NA", "1", "T1")]), ["mutect2", "strelka", "manta"]
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

    def test_known_sites_for_a_custom_genome_are_copied(self):
        # The references library offers VCFs only as germline_resource.vcf.gz, so a custom
        # genome's dbSNP and known indels always share a name
        dbsnp = "s3://refs/germline_resource/dbsnp/germline_resource.vcf.gz"
        known_indels = "s3://refs/germline_resource/indels/germline_resource.vcf.gz"
        params = {"input": self.INPUT, "dbsnp": dbsnp, "known_indels": known_indels}
        self.assertEqual(preprocess.colliding_vcf_copies(params), [
            ("dbsnp", dbsnp, f"{CONFIG}/dbsnp_germline_resource.vcf.gz"),
            ("known_indels", known_indels, f"{CONFIG}/known_indels_germline_resource.vcf.gz"),
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
