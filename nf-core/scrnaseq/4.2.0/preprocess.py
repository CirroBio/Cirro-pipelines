#!/usr/bin/env python3

from __future__ import annotations

import json
import urllib.request
from typing import TYPE_CHECKING

import pandas as pd

# The SDK is imported in __main__, so test_preprocess.py runs without it
if TYPE_CHECKING:
    from cirro.helpers.preprocess_dataset import PreprocessDataset


def resolve_references(ds: PreprocessDataset, *params: str):
    """Make reference params absolute against the references bucket.

    Form values name a path under that bucket rather than a full URI, so the
    configuration is not tied to one deployment. A value that is already absolute
    is left alone, which keeps datasets created before that change re-runnable.
    """
    for param in params:
        value = ds.params.get(param)
        if not isinstance(value, str) or not value or value.startswith("s3://"):
            continue
        ds.add_param(param, f"{ds.references_base}/{value}", overwrite=True)


# Form params that exist only to drive the reference picker. None of them are
# nf-core/scrnaseq parameters, so every one is removed before launch.
SELECTION_PARAMS = [
    f"{aligner}_{suffix}"
    for aligner in ("simpleaf", "kallisto", "star", "cellranger", "cellrangerarc", "cellrangermulti")
    for suffix in ("genome_source", "genome")
] + [
    "simpleaf_ref", "kb_ref", "star_ref", "cellranger_ref", "cellranger_prebuilt",
    "cellrangerarc_ref", "cellrangermulti_gex_ref", "cellrangermulti_vdj_ref",
    "cellrangermulti_prebuilt",
]

# Index params the pipeline understands. Any that the active aligner does not
# set are removed, so a stale value from another branch cannot leak through.
INDEX_PARAMS = [
    "simpleaf_index", "kallisto_index", "star_index", "cellranger_index",
    "cellrangerarc_reference", "cellranger_vdj_index", "txp2gene", "kb_t1c", "kb_t2c",
    "fasta", "gtf",
]

# Cell Ranger ARC and Cell Ranger multi analyze the libraries of a sample together, and
# tell them apart by this samplesheet column, which the form cannot require for these
# aligners alone.
LIBRARY_TYPE_COLUMN = {
    "cellrangerarc": "sample_type",
    "cellrangermulti": "feature_type",
}


def read_samplesheet(path) -> pd.DataFrame:
    """Read the samplesheet filled out in the form, keeping every value as text."""

    # The form field is marked required, but an empty table can still be submitted,
    # in which case Cirro either omits the param or writes an empty file
    if not path:
        raise ValueError("No samplesheet was provided: add at least one row of FASTQ files")
    try:
        # Read every value as text, so that sample names like "01" are kept as-is
        samplesheet = pd.read_csv(path, dtype=str, keep_default_na=False)
    except pd.errors.EmptyDataError:
        samplesheet = pd.DataFrame()
    if samplesheet.empty:
        raise ValueError("The samplesheet is empty: add at least one row of FASTQ files")
    return samplesheet


def is_empty_table(path) -> bool:
    """Whether an optional table filled out in the form has no rows.

    Cirro omits the param when the field was never touched, and writes an empty file
    (not even a header) when the table has no rows.
    """
    if not path:
        return True
    try:
        return pd.read_csv(path, dtype=str, keep_default_na=False).empty
    except pd.errors.EmptyDataError:
        return True


def check_library_types(samplesheet: pd.DataFrame, aligner: str):
    """Fail unless every row says what kind of library it is, when the aligner needs to know.

    nf-core/scrnaseq silently drops a Cell Ranger multi row that has no feature_type.
    """

    column = LIBRARY_TYPE_COLUMN.get(aligner)
    if column is None:
        return

    untyped = samplesheet.loc[samplesheet[column] == "", "sample"]
    if not untyped.empty:
        raise ValueError(
            f"Aligner {aligner} needs a {column} on every row; "
            f"check sample(s) {', '.join(untyped.unique())}"
        )

    if aligner == "cellrangerarc":
        is_atac = samplesheet["sample_type"] == "atac"
        no_third_read = samplesheet.loc[is_atac & (samplesheet["fastq_barcode"] == ""), "sample"]
        if not no_third_read.empty:
            raise ValueError(
                "Every ATAC row needs its third read (R3) in fastq_barcode; "
                f"check sample(s) {', '.join(no_third_read.unique())}"
            )


def find_reference_dir(ds: PreprocessDataset, s3, base: str):
    """Locate the Cell Ranger reference package inside a folder.

    Cell Ranger reference directories are named at build time
    (--cellranger_reference_name), and an uploaded reference can use any name
    at all, so the directory is found by the reference.json every package
    carries rather than by assuming a name. The folder may be the whole dataset
    or the package itself.
    """

    bucket, _, key = base.removeprefix("s3://").partition("/")
    prefix = key.rstrip("/") + "/"

    for page in s3.get_paginator("list_objects_v2").paginate(Bucket=bucket, Prefix=prefix):
        for obj in page.get("Contents", []):
            if obj["Key"].endswith("/reference.json"):
                found = f"s3://{bucket}/{obj['Key'][:-len('/reference.json')]}"
                ds.logger.info(f"Found Cell Ranger reference package at {found}")
                return found

    raise ValueError(
        f"No reference.json found under {base} - the folder does not contain a "
        "Cell Ranger reference package"
    )


def resolve_reference(ds: PreprocessDataset, aligner: str, s3):
    """Promote the selected reference onto the params nf-core/scrnaseq reads."""

    source = ds.params.get(f"{aligner}_genome_source")
    ds.logger.info(f"aligner={aligner} genome_source={source}")

    if source == "igenomes":
        genome = ds.params.get(f"{aligner}_genome")
        if genome is None:
            raise ValueError(f"genome_source=igenomes but no genome selected for aligner {aligner}")
        ds.add_param("genome", genome)
        ds.logger.info(f"Using iGenomes reference {genome}")
        return []

    ds.remove_param("igenomes_base", force=True)

    if source == "prebuilt":
        prebuilt = ds.params.get(
            "cellranger_prebuilt" if aligner == "cellranger" else "cellrangermulti_prebuilt"
        )
        if prebuilt is None:
            raise ValueError(f"genome_source=prebuilt but no reference selected for aligner {aligner}")
        ds.add_param("cellranger_index", prebuilt)
        ds.logger.info(f"Using 10X curated reference {prebuilt}")
        return ["cellranger_index"]

    if source != "dataset":
        raise ValueError(f"Unrecognized genome_source {source!r} for aligner {aligner}")

    # The simpleaf, kb and STAR paths below are relative to the root of the index
    # dataset, which is what picking the whole dataset in the folder picker gives
    if aligner == "simpleaf":
        ref = require_ref(ds, "simpleaf_ref", aligner)
        ds.add_param("simpleaf_index", f"{ref}/simpleaf/index")
        ds.add_param("txp2gene", f"{ref}/simpleaf/ref/t2g_3col.tsv")
        return ["simpleaf_index", "txp2gene"] + add_genome_files(ds, ref)

    if aligner == "kallisto":
        ref = require_ref(ds, "kb_ref", aligner)
        ds.add_param("kallisto_index", f"{ref}/index.idx")
        ds.add_param("txp2gene", f"{ref}/t2g.txt")
        kept = ["kallisto_index", "txp2gene"]
        if ds.params.get("kb_workflow", "standard") != "standard":
            ds.add_param("kb_t1c", f"{ref}/cdna_t2c.txt")
            ds.add_param("kb_t2c", f"{ref}/intron_t2c.txt")
            kept += ["kb_t1c", "kb_t2c"]
        return kept + add_genome_files(ds, ref)

    if aligner == "star":
        ref = require_ref(ds, "star_ref", aligner)
        ds.add_param("star_index", ref)
        return ["star_index"] + add_genome_files(ds, ref)

    if aligner == "cellranger":
        ref = require_ref(ds, "cellranger_ref", aligner)
        ds.add_param("cellranger_index", find_reference_dir(ds, s3, ref))
        return ["cellranger_index"]

    if aligner == "cellrangerarc":
        # A pre-built ARC package travels as cellranger_index, the param every prepared
        # Cell Ranger reference uses. cellrangerarc_reference is not a path -- it names a
        # genome inside a mkref config, and only applies when the pipeline builds its own.
        #
        # This route does not run on nf-core/scrnaseq 4.2.0, for two upstream reasons.
        # CELLRANGERARC_COUNT declares `tuple val(meta), path(reference)` while the
        # workflow passes cellranger_index as a bare path, so the reference arrives null;
        # and handing over fasta/gtf instead only reaches CELLRANGERARC_MKGTF, whose
        # output mkref rejects as invalid. The correct param is set here so the route
        # works as soon as the pipeline is fixed.
        ref = require_ref(ds, "cellrangerarc_ref", aligner)
        ds.add_param("cellranger_index", find_reference_dir(ds, s3, ref))
        return ["cellranger_index"]

    if aligner == "cellrangermulti":
        ref = require_ref(ds, "cellrangermulti_gex_ref", aligner)
        ds.add_param("cellranger_index", find_reference_dir(ds, s3, ref))
        kept = ["cellranger_index"]

        vdj = ds.params.get("cellrangermulti_vdj_ref")
        if vdj:
            ds.add_param("cellranger_vdj_index", find_reference_dir(ds, s3, vdj))
            kept.append("cellranger_vdj_index")
        else:
            ds.logger.info("No V(D)J reference selected; the pipeline will build or skip its own")
        return kept

    raise ValueError(f"Unrecognized aligner {aligner!r}")


def require_ref(ds: PreprocessDataset, param: str, aligner: str) -> str:
    ref = ds.params.get(param)
    if not ref:
        raise ValueError(f"genome_source=dataset but {param} is not set for aligner {aligner}")
    return ref.rstrip("/")


def add_genome_files(ds: PreprocessDataset, ref: str):
    """Every nf-index-genome build publishes the FASTA and GTF alongside the index."""
    ds.add_param("fasta", f"{ref}/genome.fasta")
    ds.add_param("gtf", f"{ref}/genome.gtf")
    return ["fasta", "gtf"]


def filter_params_by_schema(ds: PreprocessDataset):
    """Remove any params not present in the nf-core/scrnaseq nextflow_schema.json."""

    version = ds.params.get("workflow_version", "4.2.0")
    url = f"https://raw.githubusercontent.com/nf-core/scrnaseq/{version}/nextflow_schema.json"

    ds.logger.info(f"Fetching nextflow_schema.json for nf-core/scrnaseq {version}")
    try:
        with urllib.request.urlopen(url) as response:
            schema = json.loads(response.read().decode())
    except Exception as e:
        ds.logger.warning(f"Could not fetch nextflow_schema.json: {e}")
        return

    allowed = set()
    for section in {**schema.get("$defs", {}), **schema.get("definitions", {})}.values():
        allowed.update(section.get("properties", {}).keys())

    ds.logger.info(f"Schema defines {len(allowed):,} parameters")

    for key in list(ds.params.keys()):
        if key not in allowed:
            ds.logger.info(f"Removing param not in schema: {key}")
            ds.remove_param(key, force=True)


if __name__ == "__main__":

    # Imported here so that test_preprocess.py runs without the SDK or boto3
    import boto3
    from cirro.helpers.preprocess_dataset import PreprocessDataset

    ds = PreprocessDataset.from_running()

    aligner = ds.params["aligner"]

    check_library_types(read_samplesheet(ds.params.get("input")), aligner)

    # The pipeline checks the multiplexed samples itself, but a file with no header
    # is not a barcode CSV at all
    if is_empty_table(ds.params.get("cellranger_multi_barcodes")):
        ds.remove_param("cellranger_multi_barcodes", force=True)

    resolve_references(ds, "cellranger_prebuilt", "cellrangermulti_prebuilt")

    kept = resolve_reference(ds, aligner, boto3.client("s3"))

    for param in INDEX_PARAMS:
        if param not in kept:
            ds.remove_param(param, force=True)

    for param in SELECTION_PARAMS:
        ds.remove_param(param, force=True)

    filter_params_by_schema(ds)
