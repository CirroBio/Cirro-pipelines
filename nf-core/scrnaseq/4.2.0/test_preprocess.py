"""Offline tests for the nf-core/scrnaseq 4.2.0 preprocess script.

preprocess.py imports the Cirro SDK and boto3 only when it runs, and takes the S3
client as an argument, so this is a plain `python -m unittest` that needs pandas
and nothing else.

The samplesheet is the structured input Cirro writes from the run form: a CSV with
CRLF line endings, a header from the form's columns, and blank cells left empty. A
reference param is the value of a folder picker: the S3 URI of the picked dataset's
data/ folder, or of a folder inside it, with no trailing slash.

Beyond the samplesheet, the point is the part that is easy to get wrong and
expensive to discover on the platform -- that each aligner promotes the right
index params, that no other branch's values leak through, and that a Cell Ranger
reference is found by content rather than by an assumed directory name.
"""
import importlib.util
import io
import unittest
from pathlib import Path

import pandas as pd


def _load_preprocess():
    path = Path(__file__).with_name("preprocess.py")
    spec = importlib.util.spec_from_file_location("scrnaseq_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()

HEADER = "sample,fastq_1,fastq_2,fastq_barcode,expected_cells,seq_center,sample_type,feature_type\r\n"
COLUMNS = HEADER.strip().split(",")

REF = "s3://bucket/datasets/abc/data"
VDJ_REF = "s3://bucket/datasets/vdj/data"


def table(*rows):
    """Rows as read_samplesheet returns them: every column present, blank cells as ''."""
    return pd.DataFrame([{column: row.get(column, "") for column in COLUMNS} for row in rows])


class FakeLogger:
    def __init__(self):
        self.messages = []

    def info(self, message):
        self.messages.append(str(message))

    def warning(self, message):
        self.messages.append(str(message))


class FakeDataset:
    def __init__(self, params):
        self.params = dict(params)
        self.logger = FakeLogger()

    def add_param(self, key, value, overwrite=False):
        self.params[key] = value

    def remove_param(self, key, force=False):
        self.params.pop(key, None)


class FakeS3:
    """Holds the given objects and lists them the way boto3's list_objects_v2 paginator does."""

    def __init__(self, *uris):
        self.objects = [uri.removeprefix("s3://").split("/", 1) for uri in uris]

    def get_paginator(self, operation):
        assert operation == "list_objects_v2"
        return self

    def paginate(self, Bucket, Prefix):
        contents = [
            {"Key": key} for bucket, key in self.objects
            if bucket == Bucket and key.startswith(Prefix)
        ]
        # S3 leaves Contents out of a page that has no objects
        yield {"Contents": contents} if contents else {}


def reference_package(folder, name="cellranger_reference"):
    """The objects of a Cell Ranger reference package named `name` inside `folder`."""
    return [
        f"{folder}/{name}/reference.json",
        f"{folder}/{name}/fasta/genome.fa",
        f"{folder}/{name}/star/SA",
    ]


def resolve(aligner, params, s3=FakeS3()):
    """Run resolve_reference plus the cleanup the entrypoint does after it."""
    ds = FakeDataset(params)
    kept = preprocess.resolve_reference(ds, aligner, s3)
    for param in preprocess.INDEX_PARAMS:
        if param not in kept:
            ds.remove_param(param, force=True)
    for param in preprocess.SELECTION_PARAMS:
        ds.remove_param(param, force=True)
    return ds, kept


class ReadSamplesheet(unittest.TestCase):

    def test_values_stay_text(self):
        samplesheet = preprocess.read_samplesheet(io.StringIO(
            HEADER + "01,s3://bucket/data/01_R1.fastq.gz,s3://bucket/data/01_R2.fastq.gz,,3000,,,\r\n"
        ))
        self.assertEqual(samplesheet["sample"].tolist(), ["01"])
        self.assertEqual(samplesheet["expected_cells"].tolist(), ["3000"])
        self.assertEqual(samplesheet["fastq_barcode"].tolist(), [""])
        self.assertEqual(samplesheet["feature_type"].tolist(), [""])

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


class MultiplexedSamples(unittest.TestCase):

    BARCODE_HEADER = "sample,multiplexed_sample_id,probe_barcode_ids,cmo_ids,ocm_ids,description\r\n"

    def test_untouched_or_empty_table_is_empty(self):
        for path in (None, ""):
            self.assertTrue(preprocess.is_empty_table(path))
        self.assertTrue(preprocess.is_empty_table(io.StringIO("")))
        self.assertTrue(preprocess.is_empty_table(io.StringIO(self.BARCODE_HEADER)))

    def test_table_with_rows_is_kept(self):
        self.assertFalse(preprocess.is_empty_table(io.StringIO(
            self.BARCODE_HEADER + "PBMC,PBMC_1,,CMO301,,PBMCs 1\r\n"
        )))


class CheckLibraryTypes(unittest.TestCase):

    def test_single_library_aligners_ignore_the_library_columns(self):
        samplesheet = table({"sample": "Sample_X", "fastq_1": "x_R1.fastq.gz", "fastq_2": "x_R2.fastq.gz"})
        for aligner in ("simpleaf", "kallisto", "star", "cellranger"):
            preprocess.check_library_types(samplesheet, aligner)

    def test_arc_rows_with_their_types_and_third_reads_pass(self):
        preprocess.check_library_types(table(
            {"sample": "PBMC", "sample_type": "gex"},
            {"sample": "PBMC", "sample_type": "atac", "fastq_barcode": "atac_R3.fastq.gz"},
        ), "cellrangerarc")

    def test_arc_row_without_a_sample_type_fails_naming_its_sample(self):
        with self.assertRaisesRegex(ValueError, "sample_type on every row; check sample\\(s\\) PBMC2"):
            preprocess.check_library_types(table(
                {"sample": "PBMC1", "sample_type": "gex"},
                {"sample": "PBMC2"},
            ), "cellrangerarc")

    def test_atac_row_without_its_third_read_fails(self):
        with self.assertRaisesRegex(ValueError, "third read \\(R3\\).*PBMC"):
            preprocess.check_library_types(table(
                {"sample": "PBMC", "sample_type": "atac"},
            ), "cellrangerarc")

    def test_multi_rows_with_their_feature_types_pass(self):
        preprocess.check_library_types(table(
            {"sample": "PBMC_10K", "feature_type": "gex"},
            {"sample": "PBMC_10K", "feature_type": "vdj"},
            {"sample": "PBMC_10K", "feature_type": "ab"},
        ), "cellrangermulti")

    def test_multi_row_without_a_feature_type_fails(self):
        with self.assertRaisesRegex(ValueError, "feature_type on every row; check sample\\(s\\) PBMC_10K"):
            preprocess.check_library_types(table(
                {"sample": "PBMC_10K", "feature_type": "gex"},
                {"sample": "PBMC_10K"},
            ), "cellrangermulti")


class SimpleafRoute(unittest.TestCase):

    def test_promotes_index_and_t2g(self):
        ds, kept = resolve("simpleaf", {
            "simpleaf_genome_source": "dataset", "simpleaf_ref": REF})
        self.assertEqual(ds.params["simpleaf_index"], f"{REF}/simpleaf/index")
        self.assertEqual(ds.params["txp2gene"], f"{REF}/simpleaf/ref/t2g_3col.tsv")
        self.assertEqual(ds.params["fasta"], f"{REF}/genome.fasta")
        self.assertEqual(ds.params["gtf"], f"{REF}/genome.gtf")
        self.assertIn("simpleaf_index", kept)

    def test_no_other_aligner_index_survives(self):
        ds, _ = resolve("simpleaf", {
            "simpleaf_genome_source": "dataset", "simpleaf_ref": REF,
            "kallisto_index": "leftover", "star_index": "leftover",
            "cellranger_index": "leftover"})
        for leaked in ("kallisto_index", "star_index", "cellranger_index",
                       "cellrangerarc_reference", "cellranger_vdj_index"):
            self.assertNotIn(leaked, ds.params)


class KallistoRoute(unittest.TestCase):

    def test_standard_workflow_omits_velocity_maps(self):
        ds, _ = resolve("kallisto", {
            "kallisto_genome_source": "dataset", "kb_ref": REF,
            "kb_workflow": "standard"})
        self.assertEqual(ds.params["kallisto_index"], f"{REF}/index.idx")
        self.assertEqual(ds.params["txp2gene"], f"{REF}/t2g.txt")
        self.assertNotIn("kb_t1c", ds.params)
        self.assertNotIn("kb_t2c", ds.params)

    def test_nac_workflow_adds_velocity_maps(self):
        ds, kept = resolve("kallisto", {
            "kallisto_genome_source": "dataset", "kb_ref": REF,
            "kb_workflow": "nac"})
        self.assertEqual(ds.params["kb_t1c"], f"{REF}/cdna_t2c.txt")
        self.assertEqual(ds.params["kb_t2c"], f"{REF}/intron_t2c.txt")
        self.assertIn("kb_t1c", kept)

    def test_missing_kb_workflow_defaults_to_standard(self):
        ds, _ = resolve("kallisto", {
            "kallisto_genome_source": "dataset", "kb_ref": REF})
        self.assertNotIn("kb_t1c", ds.params)


class StarRoute(unittest.TestCase):

    def test_index_is_the_dataset_root(self):
        ds, _ = resolve("star", {
            "star_genome_source": "dataset", "star_ref": REF})
        self.assertEqual(ds.params["star_index"], REF)
        self.assertEqual(ds.params["gtf"], f"{REF}/genome.gtf")

    def test_trailing_slash_does_not_double(self):
        ds, _ = resolve("star", {
            "star_genome_source": "dataset", "star_ref": REF + "/"})
        self.assertEqual(ds.params["star_index"], REF)
        self.assertEqual(ds.params["fasta"], f"{REF}/genome.fasta")


class CellRangerRoute(unittest.TestCase):

    def test_reference_found_by_content(self):
        ds, _ = resolve("cellranger", {
            "cellranger_genome_source": "dataset", "cellranger_ref": REF},
            FakeS3(*reference_package(REF)))
        self.assertEqual(ds.params["cellranger_index"], f"{REF}/cellranger_reference")

    def test_custom_directory_name_is_honoured(self):
        ds, _ = resolve("cellranger", {
            "cellranger_genome_source": "dataset", "cellranger_ref": REF},
            FakeS3(*reference_package(REF, "my_own_name")))
        self.assertEqual(ds.params["cellranger_index"], f"{REF}/my_own_name")

    def test_picking_the_package_folder_itself_works(self):
        package = f"{REF}/my_own_name"
        ds, _ = resolve("cellranger", {
            "cellranger_genome_source": "dataset", "cellranger_ref": package},
            FakeS3(*reference_package(REF, "my_own_name")))
        self.assertEqual(ds.params["cellranger_index"], package)

    def test_fasta_and_gtf_not_passed_with_a_reference_package(self):
        ds, _ = resolve("cellranger", {
            "cellranger_genome_source": "dataset", "cellranger_ref": REF},
            FakeS3(*reference_package(REF)))
        self.assertNotIn("fasta", ds.params)
        self.assertNotIn("gtf", ds.params)

    def test_prebuilt_source_uses_the_10X_package_verbatim(self):
        package = "s3://test-references/cellranger/refdata-gex-GRCh38-2024-A"
        ds, _ = resolve("cellranger", {
            "cellranger_genome_source": "prebuilt",
            "cellranger_prebuilt": package})
        self.assertEqual(ds.params["cellranger_index"], package)
        self.assertNotIn("igenomes_base", ds.params)


class CellRangerArcRoute(unittest.TestCase):

    def test_promotes_arc_reference_as_a_prebuilt_index(self):
        # cellrangerarc_reference names a genome in a mkref config, so a pre-built
        # package travels as cellranger_index. scrnaseq 4.2.0 cannot consume it -- see
        # the note in resolve_reference -- but that is the param to set.
        ds, kept = resolve("cellrangerarc", {
            "cellrangerarc_genome_source": "dataset", "cellrangerarc_ref": REF},
            FakeS3(*reference_package(REF, "cellrangerarc_reference")))
        self.assertEqual(ds.params["cellranger_index"], f"{REF}/cellrangerarc_reference")
        self.assertEqual(kept, ["cellranger_index"])
        self.assertNotIn("cellrangerarc_reference", ds.params)
        # The reference package carries the annotation, so no loose fasta/gtf is passed.
        self.assertNotIn("fasta", ds.params)
        self.assertNotIn("gtf", ds.params)


class CellRangerMultiRoute(unittest.TestCase):

    def test_both_references_promoted(self):
        ds, kept = resolve("cellrangermulti", {
            "cellrangermulti_genome_source": "dataset",
            "cellrangermulti_gex_ref": REF,
            "cellrangermulti_vdj_ref": VDJ_REF},
            FakeS3(*reference_package(REF), *reference_package(VDJ_REF)))
        self.assertEqual(ds.params["cellranger_index"], f"{REF}/cellranger_reference")
        self.assertEqual(ds.params["cellranger_vdj_index"], f"{VDJ_REF}/cellranger_reference")
        self.assertIn("cellranger_vdj_index", kept)

    def test_vdj_is_optional(self):
        for unset in ({}, {"cellrangermulti_vdj_ref": None}):
            ds, kept = resolve("cellrangermulti", {
                "cellrangermulti_genome_source": "dataset",
                "cellrangermulti_gex_ref": REF, **unset},
                FakeS3(*reference_package(REF)))
            self.assertIn("cellranger_index", ds.params)
            self.assertNotIn("cellranger_vdj_index", ds.params)
            self.assertNotIn("cellranger_vdj_index", kept)


class IgenomesRoute(unittest.TestCase):

    def test_promotes_genome_and_keeps_igenomes_base(self):
        ds, kept = resolve("star", {
            "star_genome_source": "igenomes", "star_genome": "GRCm38",
            "igenomes_base": "s3://test-references/igenomes/"})
        self.assertEqual(ds.params["genome"], "GRCm38")
        self.assertEqual(ds.params["igenomes_base"], "s3://test-references/igenomes/")
        self.assertEqual(kept, [])

    def test_no_index_params_are_passed(self):
        ds, _ = resolve("cellranger", {
            "cellranger_genome_source": "igenomes", "cellranger_genome": "GRCh38",
            "cellranger_index": "leftover"})
        self.assertNotIn("cellranger_index", ds.params)


class Failures(unittest.TestCase):

    def test_dataset_source_without_a_reference_raises(self):
        with self.assertRaisesRegex(ValueError, "simpleaf_ref"):
            resolve("simpleaf", {"simpleaf_genome_source": "dataset"})

    def test_igenomes_without_a_genome_raises(self):
        with self.assertRaises(ValueError):
            resolve("star", {"star_genome_source": "igenomes"})

    def test_unknown_genome_source_raises(self):
        with self.assertRaisesRegex(ValueError, "somewhere_else"):
            resolve("star", {"star_genome_source": "somewhere_else"})

    def test_unknown_aligner_raises(self):
        with self.assertRaises(ValueError):
            resolve("bowtie2", {"bowtie2_genome_source": "dataset"})


class ScaffoldingCleanup(unittest.TestCase):

    def test_every_selection_param_is_removed(self):
        params = {name: "x" for name in preprocess.SELECTION_PARAMS}
        params.update({"star_genome_source": "dataset", "star_ref": REF})
        ds, _ = resolve("star", params)
        for name in preprocess.SELECTION_PARAMS:
            self.assertNotIn(name, ds.params, f"{name} leaked into the pipeline params")

    def test_selection_params_cover_every_aligner(self):
        for aligner in ("simpleaf", "kallisto", "star", "cellranger",
                        "cellrangerarc", "cellrangermulti"):
            self.assertIn(f"{aligner}_genome_source", preprocess.SELECTION_PARAMS)
            self.assertIn(f"{aligner}_genome", preprocess.SELECTION_PARAMS)


class ReferenceScan(unittest.TestCase):

    def test_finds_the_package_directory(self):
        s3 = FakeS3(f"{REF}/genome.fasta", *reference_package(REF, "my_ref"))
        self.assertEqual(
            preprocess.find_reference_dir(FakeDataset({}), s3, REF), f"{REF}/my_ref")

    def test_absent_reference_json_raises(self):
        s3 = FakeS3(f"{REF}/genome.fasta")
        with self.assertRaisesRegex(ValueError, "No reference.json"):
            preprocess.find_reference_dir(FakeDataset({}), s3, REF)

    def test_scan_stays_inside_the_picked_folder(self):
        s3 = FakeS3(*reference_package(REF + "_old"))
        with self.assertRaisesRegex(ValueError, "No reference.json"):
            preprocess.find_reference_dir(FakeDataset({}), s3, REF)


if __name__ == "__main__":
    unittest.main()
