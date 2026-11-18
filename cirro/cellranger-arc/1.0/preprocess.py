#!/usr/bin/env python3

import pandas as pd

FASTQ_TYPES = ["gex", "atac"]


def absolute_reference(value, references_base: str):
    """A reference param's value, made absolute against the references bucket.

    Form values name a path under that bucket rather than a full URI, so the
    configuration is not tied to one deployment. A value that is already absolute
    is returned as is, which keeps datasets created before that change re-runnable.
    """
    if not isinstance(value, str) or not value or value.startswith("s3://"):
        return value
    return f"{references_base}/{value}"


def read_samplesheet(path) -> pd.DataFrame:
    """Read the samplesheet filled out in the form, keeping every value as text."""

    # The form field is marked required, but an empty table can still be submitted, in
    # which case Cirro either omits the param or writes an empty file. The workflow would
    # then finish successfully without running Cell Ranger ARC at all.
    if not path:
        raise ValueError("No samplesheet was provided: add at least one row of FASTQ files")
    try:
        samplesheet = pd.read_csv(path, dtype=str, keep_default_na=False)
    except pd.errors.EmptyDataError:
        samplesheet = pd.DataFrame()
    if samplesheet.empty:
        raise ValueError("The samplesheet is empty: add at least one row of FASTQ files")
    return samplesheet


def check_libraries(samplesheet: pd.DataFrame):
    """Fail unless every sample has gex and atac rows, and only atac rows have an R3 FASTQ.

    arc_count.nf skips rows of any other type, and samples missing either library,
    without an error. Cell Ranger ARC needs all three reads of an ATAC library.
    """

    unknown = samplesheet.loc[~samplesheet["fastq_type"].isin(FASTQ_TYPES), "sample"]
    if not unknown.empty:
        raise ValueError(
            f"The library type of each row must be gex or atac; check sample(s) {', '.join(unknown.unique())}"
        )

    for fastq_type in FASTQ_TYPES:
        has_type = samplesheet.loc[samplesheet["fastq_type"] == fastq_type, "sample"]
        lacking = samplesheet.loc[~samplesheet["sample"].isin(has_type), "sample"]
        if not lacking.empty:
            raise ValueError(
                f"Each sample needs both gex and atac rows; sample(s) {', '.join(lacking.unique())} "
                f"have no {fastq_type} rows"
            )

    no_r3 = samplesheet.loc[(samplesheet["fastq_type"] == "atac") & (samplesheet["fastq_3"] == ""), "sample"]
    if not no_r3.empty:
        raise ValueError(
            f"ATAC rows need all three reads, with R3 as FASTQ (R3); check sample(s) {', '.join(no_r3.unique())}"
        )

    # It would be staged with the gene expression reads
    gex_r3 = samplesheet.loc[(samplesheet["fastq_type"] == "gex") & (samplesheet["fastq_3"] != ""), "sample"]
    if not gex_r3.empty:
        raise ValueError(
            f"Only ATAC rows have an R3 read; clear FASTQ (R3) on the gex rows of sample(s) {', '.join(gex_r3.unique())}"
        )


if __name__ == "__main__":

    # Imported here so that test_preprocess.py runs without the SDK
    from cirro.helpers.preprocess_dataset import PreprocessDataset

    ds = PreprocessDataset.from_running()

    reference_dir = ds.params.get("reference_dir")
    resolved = absolute_reference(reference_dir, ds.references_base)
    if resolved != reference_dir:
        ds.add_param("reference_dir", resolved, overwrite=True)

    check_libraries(read_samplesheet(ds.params.get("samplesheet")))

    # Force params.json to be written: the HealthOmics pre-process Lambda fails the run
    # when the file is absent, and the SDK writes it only when a parameter changes.
    ds.keep_params(list(ds.params.keys()))
