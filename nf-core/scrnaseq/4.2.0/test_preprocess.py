"""Offline tests for the nf-core/scrnaseq reference resolution.

Runs without the Cirro SDK or boto3: both are stubbed before preprocess is
loaded, so this is a plain `python -m unittest`.

The point of these tests is the part that is easy to get wrong and expensive to
discover on the platform -- that each aligner promotes the right index params,
that no other branch's values leak through, and that a Cell Ranger reference is
found by content rather than by an assumed directory name.
"""
import sys
import types
import unittest
from pathlib import Path


def _load_preprocess():
    helpers = types.ModuleType("cirro.helpers.preprocess_dataset")
    helpers.PreprocessDataset = object
    s3_path_mod = types.ModuleType("cirro.api.models.s3_path")

    class S3Path:
        def __init__(self, path):
            rest = path.replace("s3://", "", 1)
            self.bucket, _, self.key = rest.partition("/")

    s3_path_mod.S3Path = S3Path

    for name, mod in [
        ("cirro", types.ModuleType("cirro")),
        ("cirro.helpers", types.ModuleType("cirro.helpers")),
        ("cirro.api", types.ModuleType("cirro.api")),
        ("cirro.api.models", types.ModuleType("cirro.api.models")),
    ]:
        sys.modules.setdefault(name, mod)
    sys.modules["cirro.helpers.preprocess_dataset"] = helpers
    sys.modules["cirro.api.models.s3_path"] = s3_path_mod
    sys.modules.setdefault("boto3", types.ModuleType("boto3"))

    import importlib.util
    path = Path(__file__).with_name("preprocess.py")
    spec = importlib.util.spec_from_file_location("scrnaseq_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()

REF = "s3://bucket/datasets/abc/data"


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


def resolve(aligner, params):
    """Run resolve_reference plus the cleanup the entrypoint does after it."""
    ds = FakeDataset(params)
    kept = preprocess.resolve_reference(ds, aligner)
    for param in preprocess.INDEX_PARAMS:
        if param not in kept:
            ds.remove_param(param, force=True)
    for param in preprocess.SELECTION_PARAMS:
        ds.remove_param(param, force=True)
    return ds, kept


def _stub_reference_scan(found="cellranger_reference"):
    """Replace the S3 scan with one that reports a package at <base>/<found>."""
    def fake(ds, base):
        if found is None:
            raise ValueError(f"No reference.json found under {base}")
        return f"{base.rstrip('/')}/{found}"
    return fake


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
        ds, kept = resolve("kallisto", {
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

    def setUp(self):
        self._real = preprocess.find_reference_dir
        preprocess.find_reference_dir = _stub_reference_scan()

    def tearDown(self):
        preprocess.find_reference_dir = self._real

    def test_reference_found_by_content(self):
        ds, _ = resolve("cellranger", {
            "cellranger_genome_source": "dataset", "cellranger_ref": REF})
        self.assertEqual(ds.params["cellranger_index"], f"{REF}/cellranger_reference")

    def test_custom_directory_name_is_honoured(self):
        preprocess.find_reference_dir = _stub_reference_scan("my_own_name")
        ds, _ = resolve("cellranger", {
            "cellranger_genome_source": "dataset", "cellranger_ref": REF})
        self.assertEqual(ds.params["cellranger_index"], f"{REF}/my_own_name")

    def test_fasta_and_gtf_not_passed_with_a_reference_package(self):
        ds, _ = resolve("cellranger", {
            "cellranger_genome_source": "dataset", "cellranger_ref": REF})
        self.assertNotIn("fasta", ds.params)
        self.assertNotIn("gtf", ds.params)

    def test_prebuilt_source_uses_the_10X_package_verbatim(self):
        package = "s3://pubweb-references/cellranger/refdata-gex-GRCh38-2024-A"
        ds, _ = resolve("cellranger", {
            "cellranger_genome_source": "prebuilt",
            "cellranger_prebuilt": package})
        self.assertEqual(ds.params["cellranger_index"], package)
        self.assertNotIn("igenomes_base", ds.params)


class CellRangerArcRoute(unittest.TestCase):

    def setUp(self):
        self._real = preprocess.find_reference_dir
        preprocess.find_reference_dir = _stub_reference_scan("cellrangerarc_reference")

    def tearDown(self):
        preprocess.find_reference_dir = self._real

    def test_promotes_arc_reference_as_a_prebuilt_index(self):
        # cellrangerarc_reference names a genome in a mkref config, so a pre-built
        # package travels as cellranger_index. scrnaseq 4.2.0 cannot consume it -- see
        # the note in resolve_reference -- but that is the param to set.
        ds, kept = resolve("cellrangerarc", {
            "cellrangerarc_genome_source": "dataset", "cellrangerarc_ref": REF})
        self.assertEqual(ds.params["cellranger_index"],
                         f"{REF}/cellrangerarc_reference")
        self.assertEqual(kept, ["cellranger_index"])
        self.assertNotIn("cellrangerarc_reference", ds.params)
        # The reference package carries the annotation, so no loose fasta/gtf is passed.
        self.assertNotIn("fasta", ds.params)
        self.assertNotIn("gtf", ds.params)


class CellRangerMultiRoute(unittest.TestCase):

    def setUp(self):
        self._real = preprocess.find_reference_dir
        preprocess.find_reference_dir = _stub_reference_scan()

    def tearDown(self):
        preprocess.find_reference_dir = self._real

    def test_both_references_promoted(self):
        ds, kept = resolve("cellrangermulti", {
            "cellrangermulti_genome_source": "dataset",
            "cellrangermulti_gex_ref": REF,
            "cellrangermulti_vdj_ref": "s3://bucket/vdj/data"})
        self.assertEqual(ds.params["cellranger_index"], f"{REF}/cellranger_reference")
        self.assertEqual(ds.params["cellranger_vdj_index"],
                         "s3://bucket/vdj/data/cellranger_reference")
        self.assertIn("cellranger_vdj_index", kept)

    def test_vdj_is_optional(self):
        ds, kept = resolve("cellrangermulti", {
            "cellrangermulti_genome_source": "dataset",
            "cellrangermulti_gex_ref": REF})
        self.assertIn("cellranger_index", ds.params)
        self.assertNotIn("cellranger_vdj_index", ds.params)
        self.assertNotIn("cellranger_vdj_index", kept)


class IgenomesRoute(unittest.TestCase):

    def test_promotes_genome_and_keeps_igenomes_base(self):
        ds, kept = resolve("star", {
            "star_genome_source": "igenomes", "star_genome": "GRCm38",
            "igenomes_base": "s3://pubweb-references/igenomes/"})
        self.assertEqual(ds.params["genome"], "GRCm38")
        self.assertEqual(ds.params["igenomes_base"], "s3://pubweb-references/igenomes/")
        self.assertEqual(kept, [])

    def test_no_index_params_are_passed(self):
        ds, _ = resolve("cellranger", {
            "cellranger_genome_source": "igenomes", "cellranger_genome": "GRCh38",
            "cellranger_index": "leftover"})
        self.assertNotIn("cellranger_index", ds.params)


class Failures(unittest.TestCase):

    def test_dataset_source_without_a_reference_raises(self):
        with self.assertRaises(ValueError) as caught:
            resolve("simpleaf", {"simpleaf_genome_source": "dataset"})
        self.assertIn("simpleaf_ref", str(caught.exception))

    def test_igenomes_without_a_genome_raises(self):
        with self.assertRaises(ValueError):
            resolve("star", {"star_genome_source": "igenomes"})

    def test_unknown_genome_source_raises(self):
        with self.assertRaises(ValueError) as caught:
            resolve("star", {"star_genome_source": "somewhere_else"})
        self.assertIn("somewhere_else", str(caught.exception))

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
    """find_reference_dir walks the real pagination shape boto3 returns."""

    def _paginated(self, keys):
        class Paginator:
            def paginate(self, Bucket, Prefix):
                yield {"Contents": [{"Key": k} for k in keys]}

        class Client:
            def get_paginator(self, _name):
                return Paginator()

        preprocess.boto3.client = lambda _service: Client()

    def test_finds_the_package_directory(self):
        self._paginated([
            "datasets/abc/data/genome.fasta",
            "datasets/abc/data/my_ref/star/SAindex",
            "datasets/abc/data/my_ref/reference.json",
        ])
        ds = FakeDataset({})
        self.assertEqual(
            preprocess.find_reference_dir(ds, "s3://bucket/datasets/abc/data"),
            "s3://bucket/datasets/abc/data/my_ref")

    def test_absent_reference_json_raises(self):
        self._paginated(["datasets/abc/data/genome.fasta"])
        ds = FakeDataset({})
        with self.assertRaises(ValueError) as caught:
            preprocess.find_reference_dir(ds, "s3://bucket/datasets/abc/data")
        self.assertIn("No reference.json", str(caught.exception))


if __name__ == "__main__":
    unittest.main()


class SamplesheetBuilder(unittest.TestCase):
    """The multi-library aligners need per-row metadata Cirro holds per sample."""

    def _dataset(self, params, metadata_rows, files_rows):
        import pandas as pd

        class Pivoting(FakeDataset):
            """Stands in for the SDK: rows keyed on sample, metadata merged on sample."""

            def __init__(self, params, samplesheet, files):
                super().__init__(params)
                self.samplesheet = samplesheet
                self.files = files

            def pivot_samplesheet(self, metadata_columns, file_filter_predicate):
                wide = (
                    self.files.query(file_filter_predicate)
                    .pivot_table(index="sample", columns="read", values="file",
                                 aggfunc="first")
                    .rename(columns=lambda read: f"fastq_{read}")
                    .reset_index()
                )
                combined = wide.merge(self.samplesheet, on="sample", how="inner")
                keep = ["sample"] + [c for c in combined.columns
                                     if c.startswith("fastq_")]
                keep += [c for c in metadata_columns if c in combined.columns]
                return combined[keep]

        return Pivoting(params, pd.DataFrame(metadata_rows), pd.DataFrame(files_rows))

    def _files(self, samples):
        """samples maps a sample name to how many reads its library sequenced."""
        return [
            {"sample": s, "read": read, "readType": "R", "file": f"{s}_R{read}.fastq.gz"}
            for s, n_reads in samples.items() for read in range(1, n_reads + 1)
        ]

    def test_arc_libraries_collapse_onto_the_grouping_column(self):
        ds = self._dataset(
            {},
            [{"sample": "PBMC_gex", "sample_type": "gex", "sample_group": "PBMC"},
             {"sample": "PBMC_atac", "sample_type": "atac", "sample_group": "PBMC"}],
            self._files({"PBMC_gex": 2, "PBMC_atac": 3}),
        )
        sheet = preprocess.make_samplesheet(ds, "cellrangerarc")
        self.assertEqual(set(sheet["sample"]), {"PBMC"})
        self.assertEqual(set(sheet["sample_type"]), {"gex", "atac"})
        self.assertNotIn("sample_group", sheet.columns)

    def test_multi_libraries_keep_their_feature_types(self):
        ds = self._dataset(
            {},
            [{"sample": "P_gex", "feature_type": "gex", "sample_group": "PBMC_10K"},
             {"sample": "P_vdj", "feature_type": "vdj", "sample_group": "PBMC_10K"},
             {"sample": "P_ab", "feature_type": "ab", "sample_group": "PBMC_10K"}],
            self._files({"P_gex": 2, "P_vdj": 2, "P_ab": 2}),
        )
        sheet = preprocess.make_samplesheet(ds, "cellrangermulti")
        self.assertEqual(set(sheet["sample"]), {"PBMC_10K"})
        self.assertEqual(set(sheet["feature_type"]), {"gex", "vdj", "ab"})
        self.assertEqual(len(sheet), 3)

    def test_arc_without_the_metadata_raises_naming_both_columns(self):
        ds = self._dataset({}, [{"sample": "PBMC_gex"}], self._files({"PBMC_gex": 3}))
        with self.assertRaises(ValueError) as caught:
            preprocess.make_samplesheet(ds, "cellrangerarc")
        self.assertIn("sample_type", str(caught.exception))
        self.assertIn("sample_group", str(caught.exception))

    def test_multi_without_the_grouping_column_raises(self):
        ds = self._dataset(
            {}, [{"sample": "P_gex", "feature_type": "gex"}], self._files({"P_gex": 2}))
        with self.assertRaises(ValueError) as caught:
            preprocess.make_samplesheet(ds, "cellrangermulti")
        self.assertIn("sample_group", str(caught.exception))

    def test_single_library_aligner_needs_neither_column(self):
        ds = self._dataset({}, [{"sample": "Sample_X"}], self._files({"Sample_X": 2}))
        sheet = preprocess.make_samplesheet(ds, "simpleaf")
        self.assertEqual(list(sheet["sample"]), ["Sample_X"])
        self.assertEqual(sheet["fastq_1"].iloc[0], "Sample_X_R1.fastq.gz")

    def test_optional_columns_match_the_pipeline_input_schema(self):
        # assets/schema_input.json for 4.2.0; sub_sample is not among them.
        self.assertEqual(
            sorted(preprocess.OPTIONAL_SAMPLESHEET_COLUMNS),
            sorted(["expected_cells", "seq_center", "feature_type",
                    "sample_type", "fastq_barcode"]),
        )

    def test_atac_barcode_read_moves_out_of_fastq_2(self):
        ds = self._dataset(
            {},
            [{"sample": "PBMC_atac", "sample_type": "atac", "sample_group": "PBMC"},
             {"sample": "PBMC_gex", "sample_type": "gex", "sample_group": "PBMC"}],
            self._files({"PBMC_atac": 3, "PBMC_gex": 2}),
        )
        sheet = preprocess.make_samplesheet(ds, "cellrangerarc").set_index("sample_type")

        # ATAC: R1 and R3 are the genomic pair, R2 is the cell barcode.
        self.assertEqual(sheet.loc["atac", "fastq_1"], "PBMC_atac_R1.fastq.gz")
        self.assertEqual(sheet.loc["atac", "fastq_2"], "PBMC_atac_R3.fastq.gz")
        self.assertEqual(sheet.loc["atac", "fastq_barcode"], "PBMC_atac_R2.fastq.gz")

        # GEX sequenced two reads and carries no barcode file.
        self.assertEqual(sheet.loc["gex", "fastq_1"], "PBMC_gex_R1.fastq.gz")
        self.assertEqual(sheet.loc["gex", "fastq_2"], "PBMC_gex_R2.fastq.gz")
        self.assertEqual(sheet.loc["gex", "fastq_barcode"], "")

        self.assertNotIn("fastq_3", sheet.columns)

    def test_arc_without_a_third_read_raises(self):
        ds = self._dataset(
            {},
            [{"sample": "PBMC_atac", "sample_type": "atac", "sample_group": "PBMC"}],
            self._files({"PBMC_atac": 2}),
        )
        with self.assertRaises(ValueError) as caught:
            preprocess.make_samplesheet(ds, "cellrangerarc")
        self.assertIn("three reads", str(caught.exception))

    def test_multi_route_leaves_read_columns_alone(self):
        ds = self._dataset(
            {},
            [{"sample": "P_gex", "feature_type": "gex", "sample_group": "PBMC_10K"}],
            self._files({"P_gex": 2}),
        )
        sheet = preprocess.make_samplesheet(ds, "cellrangermulti")
        self.assertEqual(sheet["fastq_2"].iloc[0], "P_gex_R2.fastq.gz")
        self.assertNotIn("fastq_barcode", sheet.columns)
