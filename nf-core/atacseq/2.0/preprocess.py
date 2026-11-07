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


def check_replicates(samplesheet: pd.DataFrame):
    """Fail unless the replicates of each sample are numbered 1 to n, as the pipeline requires."""

    not_numbered = samplesheet.loc[~samplesheet["replicate"].str.isdecimal(), "sample"]
    if not not_numbered.empty:
        raise ValueError(
            f"Every row needs a whole-number replicate; check sample(s) {', '.join(not_numbered.unique())}"
        )

    for sample, replicates in samplesheet["replicate"].astype(int).groupby(samplesheet["sample"]):
        numbers = sorted(replicates.unique())
        if numbers != list(range(1, len(numbers) + 1)):
            raise ValueError(
                f"The replicates of sample {sample} must be numbered 1 to {len(numbers)}, "
                f"not {', '.join(map(str, numbers))}"
            )


if __name__ == "__main__":

    # Imported here so that test_preprocess.py runs without the SDK
    from cirro.helpers.preprocess_dataset import PreprocessDataset

    ds = PreprocessDataset.from_running()

    samplesheet = read_samplesheet(ds.params.get("input"))
    check_replicates(samplesheet)

    # log
    ds.logger.info(ds.params)

    # Force params.json to be written: the HealthOmics pre-process Lambda fails the run
    # when the file is absent, and the SDK writes it only when a parameter changes.
    ds.keep_params(list(ds.params.keys()))
