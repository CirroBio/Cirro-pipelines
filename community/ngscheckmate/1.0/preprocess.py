#!/usr/bin/env python3

import pandas as pd


def read_samplesheet(path) -> pd.DataFrame:
    """Read the samplesheet filled out in the form, keeping every value as text."""

    # The form field is marked required, but an empty table can still be submitted,
    # in which case Cirro either omits the param or writes an empty file
    if not path:
        raise ValueError("No samplesheet was provided: add one row of FASTQ files per sample")
    try:
        samplesheet = pd.read_csv(path, dtype=str, keep_default_na=False)
    except pd.errors.EmptyDataError:
        samplesheet = pd.DataFrame()
    if samplesheet.empty:
        raise ValueError("The samplesheet is empty: add one row of FASTQ files per sample")
    return samplesheet


def first_pair_per_sample(samplesheet: pd.DataFrame) -> pd.DataFrame:
    """The workflow writes one result per sample name, so keep a single pair for each."""
    return samplesheet.drop_duplicates(subset="sample", keep="first")


if __name__ == "__main__":

    # Imported here so that test_preprocess.py runs without the SDK
    from cirro.helpers.preprocess_dataset import PreprocessDataset

    ds = PreprocessDataset.from_running()

    samplesheet = read_samplesheet(ds.params.get("paired_fastq_samplesheet"))
    deduplicated = first_pair_per_sample(samplesheet)
    if deduplicated.shape[0] < samplesheet.shape[0]:
        ds.logger.info("Only keeping the first pair of files for each sample")
        deduplicated.to_csv(ds.params["paired_fastq_samplesheet"], index=None)

    # Force params.json to be written: the HealthOmics pre-process Lambda fails the run
    # when the file is absent, and the SDK writes it only when a parameter changes.
    ds.keep_params(list(ds.params.keys()))
