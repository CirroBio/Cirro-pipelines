#!/usr/bin/env python3

import pandas as pd


def read_samplesheet(path) -> pd.DataFrame:
    """Read the samplesheet filled out in the form, keeping every value as text."""

    # The form field is marked required, but an empty table can still be submitted,
    # in which case Cirro either omits the param or writes an empty file
    if not path:
        raise ValueError("No samplesheet was provided: add at least one row of FASTQ files")
    try:
        # Read every value as text, so that group names like "01" are kept as-is
        samplesheet = pd.read_csv(path, dtype=str, keep_default_na=False)
    except pd.errors.EmptyDataError:
        samplesheet = pd.DataFrame()
    if samplesheet.empty:
        raise ValueError("The samplesheet is empty: add at least one row of FASTQ files")
    return samplesheet


def check_replicates(samplesheet: pd.DataFrame):
    """Fail unless the replicates of each group are numbered 1 to n, as the pipeline requires."""

    not_numbered = samplesheet.loc[~samplesheet["replicate"].str.isdecimal(), "group"]
    if not not_numbered.empty:
        raise ValueError(
            f"Every row needs a whole-number replicate; check group(s) {', '.join(not_numbered.unique())}"
        )

    for group, replicates in samplesheet["replicate"].astype(int).groupby(samplesheet["group"]):
        numbers = sorted(replicates.unique())
        if numbers != list(range(1, len(numbers) + 1)):
            raise ValueError(
                f"The replicates of group {group} must be numbered 1 to {len(numbers)}, "
                f"not {', '.join(map(str, numbers))}"
            )


def check_controls(samplesheet: pd.DataFrame, use_control: bool):
    """Fail unless each control names another group and, with use_control on, every target group has one."""

    has_control = samplesheet["control"] != ""
    controls = samplesheet.loc[has_control, "control"]

    unknown = controls.loc[~controls.isin(samplesheet["group"])]
    if not unknown.empty:
        raise ValueError(f"Control(s) {', '.join(unknown.unique())} are not a group in the samplesheet")

    own = samplesheet.loc[has_control & (samplesheet["control"] == samplesheet["group"]), "group"]
    if not own.empty:
        raise ValueError(f"Group(s) {', '.join(own.unique())} name themselves as their control")

    if use_control:
        no_control = samplesheet.loc[~has_control & ~samplesheet["group"].isin(controls), "group"]
        if not no_control.empty:
            raise ValueError(
                "Use Control is on, so every group that is not itself a control needs one; "
                f"set Control for group(s) {', '.join(no_control.unique())}"
            )


if __name__ == "__main__":

    # Imported here so that test_preprocess.py runs without the SDK
    from cirro.helpers.preprocess_dataset import PreprocessDataset

    ds = PreprocessDataset.from_running()

    samplesheet = read_samplesheet(ds.params.get("input"))
    check_replicates(samplesheet)
    check_controls(samplesheet, ds.params["use_control"])

    # macs_gsize is around 2.7e9. A whole number that large is a Long by the time
    # Cirro builds the HealthOmics parameter list, and only String, Integer, Double
    # and Boolean are handled there. A float is a Double either way.
    gsize = ds.params.get("macs_gsize")
    if gsize is not None:
        ds.add_param("macs_gsize", float(gsize), overwrite=True)

    # log
    ds.logger.info(ds.params)

    # Force params.json to be written: the HealthOmics pre-process Lambda fails the run
    # when the file is absent, and the SDK writes it only when a parameter changes.
    ds.keep_params(list(ds.params.keys()))
