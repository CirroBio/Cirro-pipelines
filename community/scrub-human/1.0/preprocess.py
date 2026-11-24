#!/usr/bin/env python3

import pandas as pd


def read_samplesheet(path) -> pd.DataFrame:
    """Read the table of FASTQ files filled out in the form, keeping every value as text."""

    # The form field is marked required, but an empty table can still be submitted,
    # in which case Cirro either omits the param or writes an empty file
    if not path:
        raise ValueError("No FASTQ files were provided: add at least one pair of FASTQ files")
    try:
        samplesheet = pd.read_csv(path, dtype=str, keep_default_na=False)
    except pd.errors.EmptyDataError:
        samplesheet = pd.DataFrame()
    if samplesheet.empty:
        raise ValueError("The table of FASTQ files is empty: add at least one pair of FASTQ files")
    return samplesheet


def file_list(samplesheet: pd.DataFrame) -> pd.Series:
    """Every FASTQ in the table, for the one-path-per-line --input_filelist."""
    return pd.Series([
        fastq
        for row in samplesheet.itertuples()
        for fastq in (row.fastq_1, row.fastq_2)
        if fastq
    ])


if __name__ == "__main__":

    # Imported here so that test_preprocess.py runs without the SDK
    from cirro.helpers.preprocess_dataset import PreprocessDataset

    ds = PreprocessDataset.from_running()

    # Write to the dataset's config/ folder (mapped in process-input.json). The workflow
    # drops the last character of every line, so the final line needs its newline too,
    # which to_csv writes.
    file_list(read_samplesheet(ds.params.get("samplesheet"))).to_csv(
        ds.params["input_filelist"],
        header=False,
        index=False
    )

    # Consumed above; not a parameter the workflow declares
    ds.remove_param("samplesheet")

    # Force params.json to be written: the HealthOmics pre-process Lambda fails the run
    # when the file is absent, and the SDK writes it only when a parameter changes.
    ds.keep_params(list(ds.params.keys()))
