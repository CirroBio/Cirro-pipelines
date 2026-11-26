#!/usr/bin/env python3

import pandas as pd


def check_params(params: dict):
    if not params.get("formula"):
        raise ValueError("Must provide formula indicating columns to test")
    if not 0 <= params["fdr_cutoff"] < 1:
        raise ValueError("fdr_cutoff must be >= 0 and < 1")


def read_samplesheet(path) -> pd.DataFrame:
    """Read the samplesheet filled out in the form, keeping every value as text."""

    # The form field is marked required, but an empty table can still be submitted,
    # in which case Cirro either omits the param or writes an empty file
    if not path:
        raise ValueError("No samplesheet was provided: add a row for each sample")
    try:
        samplesheet = pd.read_csv(path, dtype=str, keep_default_na=False)
    except pd.errors.EmptyDataError:
        samplesheet = pd.DataFrame()
    if samplesheet.empty:
        raise ValueError("The samplesheet is empty: add a row for each sample")
    return samplesheet


def add_metadata(samplesheet: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    """The samplesheet with each sample's metadata columns, which the formula and filter use.

    The metadata is a free-form table picked in the form, so its columns are checked
    here rather than by the table schema.
    """

    if "sample" not in metadata.columns:
        raise ValueError(f"The metadata table must have a 'sample' column (found: {', '.join(metadata.columns)})")
    shared = sorted(set(metadata.columns) & set(samplesheet.columns) - {"sample"})
    if shared:
        raise ValueError(f"The metadata table cannot have a column named {', '.join(shared)}")
    for name, table in [("samplesheet", samplesheet), ("metadata table", metadata)]:
        repeated = table.loc[table["sample"].duplicated(), "sample"].unique()
        if len(repeated) > 0:
            raise ValueError(f"Each sample can have only one row in the {name}, but these are repeated: {', '.join(repeated)}")
    missing = sorted(set(samplesheet["sample"]) - set(metadata["sample"]))
    if missing:
        raise ValueError(f"These samples have no row in the metadata table: {', '.join(missing)}")

    return samplesheet.merge(metadata, on="sample")


if __name__ == "__main__":

    # Imported here so that test_preprocess.py runs without the SDK
    from cirro.helpers.preprocess_dataset import PreprocessDataset

    ds = PreprocessDataset.from_running()

    check_params(ds.params)

    samplesheet = add_metadata(
        read_samplesheet(ds.params.get("samplesheet")),
        pd.read_csv(ds.params["metadata"], dtype=str, keep_default_na=False)
    )
    samplesheet.to_csv(ds.params["samplesheet"], index=None)

    # Force params.json to be written: the HealthOmics pre-process Lambda fails the run
    # when the file is absent, and the SDK writes it only when a parameter changes.
    ds.keep_params(list(ds.params.keys()))
