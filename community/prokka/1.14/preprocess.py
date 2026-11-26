#!/usr/bin/env python3

import pandas as pd


def read_samplesheet(path) -> pd.DataFrame:
    """Read the genome table filled out in the form, keeping every value as text."""

    # The form field is marked required, but an empty table can still be submitted,
    # in which case Cirro either omits the param or writes an empty file
    if not path:
        raise ValueError("No genomes were provided: add one row per genome")
    try:
        samplesheet = pd.read_csv(path, dtype=str, keep_default_na=False)
    except pd.errors.EmptyDataError:
        samplesheet = pd.DataFrame()
    if samplesheet.empty:
        raise ValueError("The genome table is empty: add one row per genome")
    return samplesheet


def check_unique_names(samplesheet: pd.DataFrame):
    """Each genome's results are written to a folder named after it."""

    duplicated = samplesheet.loc[samplesheet["name"].duplicated(), "name"].unique()
    if len(duplicated) > 0:
        raise ValueError(
            f"Only one FASTA is allowed per genome. Listed more than once: {', '.join(duplicated)}"
        )


if __name__ == "__main__":

    # Imported here so that test_preprocess.py runs without the SDK
    from cirro.helpers.preprocess_dataset import PreprocessDataset

    ds = PreprocessDataset.from_running()

    check_unique_names(read_samplesheet(ds.params.get("sample_sheet")))

    # Force params.json to be written: the HealthOmics pre-process Lambda fails the run
    # when the file is absent, and the SDK writes it only when a parameter changes.
    ds.keep_params(list(ds.params.keys()))
