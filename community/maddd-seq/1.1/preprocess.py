#!/usr/bin/env python3

import pandas as pd


def read_samplesheet(path) -> pd.DataFrame:
    """Read the samplesheet filled out in the form, keeping every value as text."""

    # The form field is marked required, but an empty table can still be submitted,
    # in which case Cirro either omits the param or writes an empty file
    if not path:
        raise ValueError("No samplesheet was provided: add one row per specimen")
    try:
        samplesheet = pd.read_csv(path, dtype=str, keep_default_na=False)
    except pd.errors.EmptyDataError:
        samplesheet = pd.DataFrame()
    if samplesheet.empty:
        raise ValueError("The samplesheet is empty: add one row per specimen")
    return samplesheet


def check_unique_specimens(samplesheet: pd.DataFrame):
    """The workflow takes a single pair of FASTQ files per specimen (FredHutch/maddd-seq
    main.nf), and repeated specimens collide in its per-specimen outputs."""

    duplicated = samplesheet.loc[samplesheet["specimen"].duplicated(), "specimen"].unique()
    if len(duplicated) > 0:
        raise ValueError(
            "Each specimen can be listed only once, with a single pair of FASTQ files "
            "(combine the files of a specimen sequenced on several lanes first). "
            f"Listed more than once: {', '.join(duplicated)}"
        )


if __name__ == "__main__":

    # Imported here so that test_preprocess.py runs without the SDK
    from cirro.helpers.preprocess_dataset import PreprocessDataset

    ds = PreprocessDataset.from_running()

    samplesheet = read_samplesheet(ds.params.get("sample_sheet"))
    ds.logger.info("Samplesheet:")
    ds.logger.info(samplesheet.to_csv(index=None))
    check_unique_specimens(samplesheet)

    ds.logger.info(ds.params)

    # Force params.json to be written: the HealthOmics pre-process Lambda fails the run
    # when the file is absent, and the SDK writes it only when a parameter changes.
    ds.keep_params(list(ds.params.keys()))
