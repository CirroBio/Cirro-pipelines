#!/usr/bin/env python3

import pandas as pd

CATEGORY = "category"


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
    """Each genome gets a contigs database named after it."""

    duplicated = samplesheet.loc[samplesheet["name"].duplicated(), "name"].unique()
    if len(duplicated) > 0:
        raise ValueError(
            f"Each genome name can be listed only once. Listed more than once: {', '.join(duplicated)}"
        )


def drop_blank_category(samplesheet: pd.DataFrame) -> pd.DataFrame:
    """Remove the category column when no row fills it in.

    Every column after the genome path is imported into anvi'o as metadata, and the
    workflow reads a blank cell as the text 'nan'.
    """
    if (samplesheet[CATEGORY] == "").all():
        return samplesheet.drop(columns=[CATEGORY])
    if samplesheet[CATEGORY].nunique() < 2:
        raise ValueError("To compare genomes by Category, it must have at least two distinct values")
    return samplesheet


if __name__ == "__main__":

    # Imported here so that test_preprocess.py runs without the SDK
    from cirro.helpers.preprocess_dataset import PreprocessDataset

    ds = PreprocessDataset.from_running()

    samplesheet = read_samplesheet(ds.params.get("sample_sheet"))
    check_unique_names(samplesheet)
    samplesheet = drop_blank_category(samplesheet)
    samplesheet.to_csv(ds.params["sample_sheet"], index=None)

    if CATEGORY in samplesheet.columns:
        ds.add_param("category_name", CATEGORY)

    # Remove any spaces from the output name
    ds.add_param(
        "output_name",
        ds.params["output_name"].replace(" ", "_"),
        overwrite=True
    )
