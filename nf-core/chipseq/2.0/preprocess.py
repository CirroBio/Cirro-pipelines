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


def check_controls(samplesheet: pd.DataFrame):
    """Fail unless each row with a control names an antibody and a replicate of another sample.

    Run after check_replicates, which makes sure every replicate is a whole number.
    """

    has_antibody = samplesheet["antibody"] != ""
    has_control = samplesheet["control"] != ""
    unpaired = samplesheet.loc[has_antibody != has_control, "sample"]
    if not unpaired.empty:
        raise ValueError(
            f"Set both antibody and control, or neither; check sample(s) {', '.join(unpaired.unique())}"
        )

    chip = samplesheet.loc[has_control]
    no_control_replicate = chip.loc[~chip["control_replicate"].str.isdecimal(), "sample"]
    if not no_control_replicate.empty:
        raise ValueError(
            "Every row with a control needs a whole-number control replicate; "
            f"check sample(s) {', '.join(no_control_replicate.unique())}"
        )

    replicates = set(zip(samplesheet["sample"], samplesheet["replicate"].astype(int)))
    for sample, control, control_replicate in zip(
        chip["sample"], chip["control"], chip["control_replicate"].astype(int)
    ):
        if control == sample:
            raise ValueError(f"Sample {sample} names itself as its control")
        if (control, control_replicate) not in replicates:
            raise ValueError(
                f"Sample {sample} names replicate {control_replicate} of {control} as its control, "
                "but the samplesheet has no such row"
            )


if __name__ == "__main__":

    # Imported here so that test_preprocess.py runs without the SDK
    from cirro.helpers.preprocess_dataset import PreprocessDataset

    ds = PreprocessDataset.from_running()

    samplesheet = read_samplesheet(ds.params.get("input"))
    ds.logger.info("Samplesheet:")
    ds.logger.info(samplesheet.to_csv(index=None))
    check_replicates(samplesheet)
    check_controls(samplesheet)

    # log
    ds.logger.info(ds.params)
