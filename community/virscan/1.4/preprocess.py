#!/usr/bin/env python3

import posixpath
import re

import pandas as pd

# phip-flow V1.14's validate-sample-table.py requires more than one beads_only row
# whenever Z-scores are fit, which process-compute.config always turns on
MIN_BEADS_ONLY = 2

FASTQ_EXTENSION = re.compile(r"\.f(ast)?q\.gz$")


def read_samplesheet(path) -> pd.DataFrame:
    """Read the sample table filled out in the form, keeping every value as text."""

    # The form field is marked required, but an empty table can still be submitted,
    # in which case Cirro either omits the param or writes an empty file
    if not path:
        raise ValueError("No sample table was provided: add one row per FASTQ file")
    try:
        samplesheet = pd.read_csv(path, dtype=str, keep_default_na=False)
    except pd.errors.EmptyDataError:
        samplesheet = pd.DataFrame()
    if samplesheet.empty:
        raise ValueError("The sample table is empty: add one row per FASTQ file")
    return samplesheet


def public_epitopes_csv(peptide_table) -> str:
    """The public epitopes published beside the selected peptide library."""

    if not peptide_table:
        raise ValueError("No peptide library was selected")
    return peptide_table.replace("library.csv", "public_epitopes.csv")


def label_controls(samplesheet: pd.DataFrame) -> pd.DataFrame:
    """Treat every row left without a control status as an empirical sample."""

    return samplesheet.assign(control_status=samplesheet["control_status"].replace("", "empirical"))


def check_beads_only(samplesheet: pd.DataFrame):
    n_beads_only = (samplesheet["control_status"] == "beads_only").sum()
    if n_beads_only < MIN_BEADS_ONLY:
        raise ValueError(
            f"Mark at least {MIN_BEADS_ONLY} rows as beads_only in the control_status column, "
            f"since the Z-score model is fit to them (found {n_beads_only})"
        )


def name_replicates(samplesheet: pd.DataFrame) -> pd.DataFrame:
    """Name each replicate left without a name after its FASTQ file."""

    file_names = samplesheet["fastq_filepath"].map(
        lambda uri: FASTQ_EXTENSION.sub("", posixpath.basename(uri))
    )
    return samplesheet.assign(
        sample_name=samplesheet["sample_name"].where(samplesheet["sample_name"] != "", file_names)
    )


def split_reads_prefix(fastq_uris: pd.Series) -> tuple[str, pd.Series]:
    """Split the FASTQ URIs into the folder they all share and each file's path under it.

    phip-flow opens every FASTQ at `<reads_prefix>/<fastq_filepath>`, so the table's
    absolute URIs have to become paths relative to a single prefix.
    """

    if (fastq_uris == "").any():
        raise ValueError("Every row of the sample table needs a FASTQ file")

    # "s3://bucket/key" splits into ["s3:", "", "bucket", ...]
    buckets = sorted({uri.split("/")[2] for uri in fastq_uris})
    if len(buckets) > 1:
        raise ValueError(
            "All FASTQ files must be in the same S3 bucket, since the workflow reads them "
            f"from one prefix. Found files in: {', '.join(buckets)}"
        )

    shared = []
    for parts in zip(*(uri.split("/")[:-1] for uri in fastq_uris)):
        if len(set(parts)) > 1:
            break
        shared.append(parts[0])
    reads_prefix = "/".join(shared)
    return reads_prefix, fastq_uris.str[len(reads_prefix) + 1:]


if __name__ == "__main__":

    # Imported here so that test_preprocess.py runs without the SDK
    from cirro.helpers.preprocess_dataset import PreprocessDataset

    ds = PreprocessDataset.from_running()

    ds.add_param("public_epitopes_csv", public_epitopes_csv(ds.params.get("peptide_table")))

    samplesheet = label_controls(read_samplesheet(ds.params.get("sample_table")))
    check_beads_only(samplesheet)
    samplesheet = name_replicates(samplesheet)
    reads_prefix, fastq_filepath = split_reads_prefix(samplesheet["fastq_filepath"])
    samplesheet = samplesheet.assign(fastq_filepath=fastq_filepath)
    ds.add_param("reads_prefix", reads_prefix)


    # Overwrite the file Cirro wrote to the dataset's config/ folder (mapped in process-input.json)
    samplesheet.to_csv(ds.params["sample_table"], index=None)

    ds.logger.info(ds.params)

    # Force params.json to be written: the HealthOmics pre-process Lambda fails the run
    # when the file is absent, and the SDK writes it only when a parameter changes.
    ds.keep_params(list(ds.params.keys()))
