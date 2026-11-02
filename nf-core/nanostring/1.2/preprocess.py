#!/usr/bin/env python3

import pandas as pd


def read_samplesheet(path) -> pd.DataFrame:
    """Read the samplesheet filled out in the form, keeping every value as text."""

    # The form field is marked required, but an empty table can still be submitted,
    # in which case Cirro either omits the param or writes an empty file
    if not path:
        raise ValueError("No samplesheet was provided: add at least one RCC file")
    try:
        samplesheet = pd.read_csv(path, dtype=str, keep_default_na=False)
    except pd.errors.EmptyDataError:
        samplesheet = pd.DataFrame()
    if samplesheet.empty:
        raise ValueError("The samplesheet is empty: add at least one RCC file")
    return samplesheet


def fill_heatmap_ids(samplesheet: pd.DataFrame) -> pd.DataFrame:
    """Label the heatmap row of every RCC file that has no label with its sample name.

    process-input.json points heatmap_id_column at HEATMAP_ID, so a blank cell would
    leave that row of the heatmap unlabelled.
    """

    labels = samplesheet["HEATMAP_ID"].mask(samplesheet["HEATMAP_ID"] == "", samplesheet["SAMPLE_ID"])
    return samplesheet.assign(HEATMAP_ID=labels)


if __name__ == "__main__":

    # Imported here so that test_preprocess.py runs without the SDK
    from cirro.helpers.preprocess_dataset import PreprocessDataset

    ds = PreprocessDataset.from_running()

    samplesheet = fill_heatmap_ids(read_samplesheet(ds.params.get("input")))
    ds.logger.info("Samplesheet:")
    ds.logger.info(samplesheet.to_csv(index=False))

    # Overwrite the file Cirro wrote to the dataset's config/ folder (mapped in process-input.json)
    samplesheet.to_csv(ds.params["input"], index=False)

    # Force params.json to be written: the HealthOmics pre-process Lambda fails the run
    # when the file is absent, and the SDK writes it only when a parameter changes.
    ds.keep_params(list(ds.params.keys()))
