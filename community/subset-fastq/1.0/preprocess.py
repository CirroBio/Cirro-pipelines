#!/usr/bin/env python3

import pandas as pd


def read_samplesheet(path) -> pd.DataFrame:
    """Read the table of FASTQ files filled out in the form, keeping every value as text."""

    # The form field is marked required, but an empty table can still be submitted,
    # in which case Cirro either omits the param or writes an empty file. The workflow
    # would then exit without subsetting anything.
    if not path:
        raise ValueError("No FASTQ files were provided: add at least one pair of FASTQ files")
    try:
        samplesheet = pd.read_csv(path, dtype=str, keep_default_na=False)
    except pd.errors.EmptyDataError:
        samplesheet = pd.DataFrame()
    if samplesheet.empty:
        raise ValueError("The table of FASTQ files is empty: add at least one pair of FASTQ files")
    return samplesheet


def fastq_list(samplesheet: pd.DataFrame) -> str:
    """The workflow takes its FASTQ files as one comma-separated list."""
    return ",".join(
        fastq
        for row in samplesheet.itertuples()
        for fastq in (row.fastq_1, row.fastq_2)
        if fastq
    )


if __name__ == "__main__":

    # Imported here so that test_preprocess.py runs without the SDK
    from cirro.helpers.preprocess_dataset import PreprocessDataset

    ds = PreprocessDataset.from_running()

    ds.add_param("input", fastq_list(read_samplesheet(ds.params.get("samplesheet"))))

    # Consumed above; not a parameter the workflow declares
    ds.remove_param("samplesheet")
