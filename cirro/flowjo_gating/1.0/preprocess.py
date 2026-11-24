#!/usr/bin/env python3

import pandas as pd


def read_fcs_files(path) -> pd.DataFrame:
    """Read the FCS file table filled out in the form, keeping every value as text."""

    # The form field is marked required, but an empty table can still be submitted,
    # in which case Cirro either omits the param or writes an empty file
    if not path:
        raise ValueError("No FCS files were provided: add one row per FCS file")
    try:
        fcs_files = pd.read_csv(path, dtype=str, keep_default_na=False)
    except pd.errors.EmptyDataError:
        fcs_files = pd.DataFrame()
    if fcs_files.empty:
        raise ValueError("The FCS file table is empty: add one row per FCS file")
    return fcs_files


def deduplicate_fcs(fcs_files: list, logger) -> str:
    """Comma-separated FCS files, keeping the first of any that share a file name.

    The workflow stages every file into one folder, so files with the same name
    collide even when they come from different folders.
    """
    seen = set()
    output = []
    for file in fcs_files:
        fn = file.split("/")[-1]
        if fn in seen:
            logger.warning(f"Removing duplicate file: {file}")
        else:
            seen.add(fn)
            output.append(file)
    return ",".join(output)


if __name__ == "__main__":

    # Imported here so that test_preprocess.py runs without the SDK
    from cirro.helpers.preprocess_dataset import PreprocessDataset

    ds = PreprocessDataset.from_running()

    fcs_files = read_fcs_files(ds.params.get("fcs_files"))
    ds.logger.info(f"Found {fcs_files.shape[0]:,} FCS files to analyze")
    ds.add_param("input_fcs", deduplicate_fcs(fcs_files["fcs"].tolist(), ds.logger))

    # Consumed above; not a parameter the workflow declares.
    ds.remove_param("fcs_files", force=True)

    # An optional pick that was left empty
    if not ds.params.get("input_anndata"):
        ds.remove_param("input_anndata", force=True)
