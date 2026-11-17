#!/usr/bin/env python3

import pandas as pd


def reference_dir(value, references_base: str):
    """A CellRanger reference directory param, made absolute against the references bucket.

    Form values name a path under that bucket rather than a full URI, so the
    configuration is not tied to one deployment. A value that is already absolute
    is left alone, which keeps datasets created before that change re-runnable.
    HealthOmics reads an S3 value with no trailing separator as an object key and
    rejects the run, so the result always ends with one.
    """
    if not value:
        return value
    if not value.startswith("s3://"):
        value = f"{references_base}/{value}"
    return value.rstrip("/") + "/"


def read_table(path, name: str, rows: str) -> pd.DataFrame:
    """Read a table filled out in the form, keeping every value as text."""

    # The form field is marked required, but an empty table can still be submitted,
    # in which case Cirro either omits the param or writes an empty file
    if not path:
        raise ValueError(f"No {name} was provided: add {rows}")
    try:
        table = pd.read_csv(path, dtype=str, keep_default_na=False)
    except pd.errors.EmptyDataError:
        table = pd.DataFrame()
    if table.empty:
        raise ValueError(f"The {name} is empty: add {rows}")
    return table


def drop_repeated_rows(samplesheet: pd.DataFrame) -> pd.DataFrame:
    """Rows filled from a dataset's file names repeat each library once per lane and
    read, and each row becomes a library in the CellRanger config."""
    return samplesheet.drop_duplicates(ignore_index=True)


def unset_umi_thresholds(params: dict) -> list:
    """The UMI thresholds left blank, which are dropped so the CellRanger default applies.

    CellRanger added emptydrops-minimum-umis in 7.1 and global-minimum-umis in 8.0,
    so a threshold that is set needs at least that version.
    """
    selected = params["cellranger_version"]
    version = tuple(int(n) for n in selected.split(".")[:2])
    unset = []
    for param, min_version in [("emptydrops_minimum_umis", (7, 1)), ("global_minimum_umis", (8, 0))]:
        if params.get(param) is None:
            unset.append(param)
        elif version < min_version:
            raise ValueError(
                f"{param} requires CellRanger {min_version[0]}.{min_version[1]} or later "
                f"(selected: {selected})"
            )
    return unset


if __name__ == "__main__":

    # Imported here so that test_preprocess.py runs without the SDK
    from cirro.helpers.preprocess_dataset import PreprocessDataset

    ds = PreprocessDataset.from_running()

    for param in ("transcriptome_dir", "vdj_dir"):
        resolved = reference_dir(ds.params.get(param), ds.references_base)
        if resolved != ds.params.get(param):
            ds.add_param(param, resolved, overwrite=True)

    samplesheet = read_table(ds.params.get("grouping"), "samplesheet", "one row per library")
    libraries = drop_repeated_rows(samplesheet)
    if len(libraries) < len(samplesheet):
        # Overwrite the file Cirro wrote to the dataset's config/ folder (mapped in process-input.json)
        libraries.to_csv(ds.params["grouping"], index=None)

    # Read only to fail clearly when the table is missing or empty
    read_table(
        ds.params.get("multiplexing"), "sample hashtag table", "one row per sample with its hashtag"
    )

    for param in unset_umi_thresholds(ds.params):
        ds.remove_param(param, force=True)

    # Left blank, the workflow's empty feature reference is used
    if ds.params.get("feature_csv") is None:
        ds.remove_param("feature_csv", force=True)

    for k, v in ds.params.items():
        ds.logger.info(f"{k}: {v}")
