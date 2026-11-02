#!/usr/bin/env python3

import pandas as pd


def read_samplesheet(path) -> pd.DataFrame:
    """Read the samplesheet filled out in the form, keeping every value as text."""

    # The form field is marked required, but an empty table can still be submitted,
    # in which case Cirro either omits the param or writes an empty file
    if not path:
        raise ValueError("No samplesheet was provided: add one row per BAM file")
    try:
        samplesheet = pd.read_csv(path, dtype=str, keep_default_na=False)
    except pd.errors.EmptyDataError:
        samplesheet = pd.DataFrame()
    if samplesheet.empty:
        raise ValueError("The samplesheet is empty: add one row per BAM file")
    return samplesheet


def check_unique_samples(samplesheet: pd.DataFrame):
    """The workflow names each sample's genotype files after it and merges them into a
    single folder, where a repeated sample name collides."""

    duplicated = samplesheet.loc[samplesheet["sample"].duplicated(), "sample"].unique()
    if len(duplicated) > 0:
        raise ValueError(
            f"Each sample can be listed only once. Listed more than once: {', '.join(duplicated)}"
        )


def format_genes(genes: list) -> str:
    """The genes selected in the form, as arcasHLA's --genes takes them."""

    if "all" in genes:
        return "all"
    if len(genes) == 0:
        raise ValueError("Must specify at least 1 gene")
    return ",".join(genes)


if __name__ == "__main__":

    # Imported here so that test_preprocess.py runs without the SDK
    from cirro.helpers.preprocess_dataset import PreprocessDataset

    ds = PreprocessDataset.from_running()

    samplesheet = read_samplesheet(ds.params.get("samplesheet"))
    check_unique_samples(samplesheet)

    ds.add_param("genes", format_genes(ds.params["genes"]), overwrite=True)

    # Force params.json to be written: the HealthOmics pre-process Lambda fails the run
    # when the file is absent, and the SDK writes it only when a parameter changes.
    ds.keep_params(list(ds.params.keys()))
