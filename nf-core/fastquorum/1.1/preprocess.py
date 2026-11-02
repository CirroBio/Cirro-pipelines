#!/usr/bin/env python3

import pandas as pd


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


def fill_read_structure(samplesheet: pd.DataFrame, read_structure) -> pd.DataFrame:
    """Give every row that leaves read_structure blank the one set in the form."""

    blank = samplesheet["read_structure"].str.strip() == ""
    if blank.any() and not read_structure:
        raise ValueError(
            "Set a read structure in the form, or in every row of the samplesheet; "
            f"missing for sample(s) {', '.join(samplesheet.loc[blank, 'sample'].unique())}"
        )
    return samplesheet.assign(read_structure=samplesheet["read_structure"].mask(blank, read_structure))


if __name__ == "__main__":

    # Imported here so that test_preprocess.py runs without the SDK
    from cirro.helpers.preprocess_dataset import PreprocessDataset

    ds = PreprocessDataset.from_running()

    samplesheet = read_samplesheet(ds.params.get("input"))
    filled = fill_read_structure(samplesheet, ds.params.get("read_structure"))

    if not filled.equals(samplesheet):
        # Overwrite the file Cirro wrote to the dataset's config/ folder (mapped in process-input.json)
        filled.to_csv(ds.params["input"], index=None)
