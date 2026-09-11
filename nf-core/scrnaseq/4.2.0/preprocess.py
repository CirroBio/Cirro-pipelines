#!/usr/bin/env python3

import json
import urllib.request

import boto3
from cirro.models.s3_path import S3Path
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

# Samplesheet columns nf-core/scrnaseq accepts beyond sample/fastq_1/fastq_2
# (assets/schema_input.json).
OPTIONAL_SAMPLESHEET_COLUMNS = [
    "expected_cells", "seq_center", "feature_type", "sample_type", "fastq_barcode",
]

# Cell Ranger ARC and Cell Ranger multi describe one biological sample with several
# libraries, told apart per samplesheet row by sample_type / feature_type. Cirro carries
# one metadata row per sample, so each library is uploaded as its own sample and this
# column names the sample those libraries belong to.
SAMPLE_GROUP_COLUMN = "sample_group"

# The per-row column whose value distinguishes the libraries of one sample.
LIBRARY_TYPE_COLUMN = {
    "cellrangerarc": "sample_type",
    "cellrangermulti": "feature_type",
}


def make_samplesheet(ds: PreprocessDataset, aligner: str):
    """Pivot the dataset files into the samplesheet nf-core/scrnaseq expects."""

    available = ds.samplesheet.columns.values
    metadata_columns = [col for col in OPTIONAL_SAMPLESHEET_COLUMNS if col in available]

    library_column = LIBRARY_TYPE_COLUMN.get(aligner)
    if library_column:
        missing = [
            col for col in (library_column, SAMPLE_GROUP_COLUMN) if col not in available
        ]
        if missing:
            raise ValueError(
                f"aligner {aligner} needs every library to carry {missing} in its sample "
                f"metadata: {library_column} says what kind of library it is, and "
                f"{SAMPLE_GROUP_COLUMN} names the sample the libraries belong to. "
                f"The samplesheet has {sorted(available)}"
            )
        metadata_columns.append(SAMPLE_GROUP_COLUMN)

    ds.logger.info(f"Carrying metadata columns into the samplesheet: {metadata_columns}")

    samplesheet = ds.pivot_samplesheet(
        metadata_columns=metadata_columns,
        file_filter_predicate='readType == "R"',
    )

    if library_column:
        # One row per library, all naming the sample they are libraries of.
        samplesheet["sample"] = samplesheet[SAMPLE_GROUP_COLUMN]
        samplesheet = samplesheet.drop(columns=[SAMPLE_GROUP_COLUMN])

    if aligner == "cellrangerarc":
        samplesheet = place_atac_barcode_read(samplesheet)

    ds.logger.info("Samplesheet:")
    for line in samplesheet.to_csv(index=None).split("\n"):
        ds.logger.info(line)

    return samplesheet


def place_atac_barcode_read(samplesheet):
    """Move ATAC rows onto the fastq_1/fastq_2/fastq_barcode layout ARC expects.

    A 10X multiome ATAC library sequences R1 and R3 as the genomic pair and R2 as the
    cell barcode, so a pivot ordered by read number leaves the barcode in fastq_2 and
    the genomic mate in fastq_3. GEX libraries have only two reads and are already in
    the right columns, keeping fastq_barcode empty.
    """

    if "fastq_3" not in samplesheet.columns:
        raise ValueError(
            "cellrangerarc needs three reads (R1, R2, R3) for each ATAC library, but "
            "no library in this dataset has a third read"
        )

    is_atac = samplesheet["sample_type"] == "atac"
    samplesheet["fastq_barcode"] = ""
    samplesheet.loc[is_atac, "fastq_barcode"] = samplesheet.loc[is_atac, "fastq_2"]
    samplesheet.loc[is_atac, "fastq_2"] = samplesheet.loc[is_atac, "fastq_3"]

    return samplesheet.drop(columns=["fastq_3"])


def find_reference_dir(ds: PreprocessDataset, base: str):
    """Locate the Cell Ranger reference package inside a dataset.

    Cell Ranger reference directories are named at build time
    (--cellranger_reference_name), and an uploaded reference can use any name
    at all, so the directory is found by the reference.json every package
    carries rather than by assuming a name.
    """

    s3_path = S3Path(base)
    prefix = s3_path.key.rstrip("/") + "/"
    paginator = boto3.client("s3").get_paginator("list_objects_v2")

    for page in paginator.paginate(Bucket=s3_path.bucket, Prefix=prefix):
        for obj in page.get("Contents", []):
            if obj["Key"].endswith("/reference.json"):
                found = f"s3://{s3_path.bucket}/{obj['Key'][:-len('/reference.json')]}"
                ds.logger.info(f"Found Cell Ranger reference package at {found}")
                return found

    raise ValueError(
        f"No reference.json found under {base} - the dataset does not contain a "
        "Cell Ranger reference package"
    )


def resolve_reference(ds: PreprocessDataset, aligner: str):
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
        ds.add_param("cellranger_index", find_reference_dir(ds, ref))
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
        ds.add_param("cellranger_index", find_reference_dir(ds, ref))
        return ["cellranger_index"]

    if aligner == "cellrangermulti":
        ref = require_ref(ds, "cellrangermulti_gex_ref", aligner)
        ds.add_param("cellranger_index", find_reference_dir(ds, ref))
        kept = ["cellranger_index"]

        vdj = ds.params.get("cellrangermulti_vdj_ref")
        if vdj:
            ds.add_param("cellranger_vdj_index", find_reference_dir(ds, vdj))
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

    ds = PreprocessDataset.from_running()

    aligner = ds.params["aligner"]

    # Write to the dataset's config/ folder (mapped in process-input.json)
    make_samplesheet(ds, aligner).to_csv(ds.params["input"], index=False)

    resolve_references(ds, "cellranger_prebuilt", "cellrangermulti_prebuilt")

    kept = resolve_reference(ds, aligner)

    for param in INDEX_PARAMS:
        if param not in kept:
            ds.remove_param(param, force=True)

    for param in SELECTION_PARAMS:
        ds.remove_param(param, force=True)

    filter_params_by_schema(ds)
