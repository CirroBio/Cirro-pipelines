#!/usr/bin/env python3

import pandas as pd


def read_table(path) -> pd.DataFrame:
    """Read a tab-separated table filled out in the form, keeping every value as text.

    An empty table reaches here either as no param at all or as an empty file.
    """
    if not path:
        return pd.DataFrame()
    try:
        return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def check_samples(samples: pd.DataFrame):
    """Fail before the run on an empty or repeated sample list: QIIME 2 imports the
    manifest only when every sample ID is unique."""

    if samples.empty:
        raise ValueError("The samples table is empty: add at least one row of HiFi FASTQ files")
    repeated = samples.loc[samples["sample-id"].duplicated(), "sample-id"].unique()
    if len(repeated) > 0:
        raise ValueError(f"Each sample can have only one row, but these are repeated: {', '.join(repeated)}")


def complete_metadata(metadata: pd.DataFrame, samples: pd.DataFrame) -> pd.DataFrame:
    """The metadata the workflow reads, which needs a sample_name and condition for
    every sample.

    Left empty, every sample gets a placeholder condition. A pool column with no values
    is dropped, since its presence alone makes the workflow split samples by pool.
    """
    if metadata.empty:
        return pd.DataFrame({
            "sample_name": samples["sample-id"],
            "condition": [f"cat_{i % 2}" for i in range(samples.shape[0])],
        })

    missing = sorted(set(samples["sample-id"]) - set(metadata["sample_name"]))
    if missing:
        raise ValueError(f"The metadata table has no row for these samples: {', '.join(missing)}")
    if "pool" in metadata.columns and (metadata["pool"] == "").all():
        metadata = metadata.drop(columns="pool")
    return metadata


if __name__ == "__main__":

    # Imported here so that test_preprocess.py runs without the SDK
    from cirro.helpers.preprocess_dataset import PreprocessDataset

    ds = PreprocessDataset.from_running()

    samples = read_table(ds.params.get("input"))
    check_samples(samples)
    metadata = complete_metadata(read_table(ds.params.get("metadata")), samples)

    # An untouched metadata table has no param, so it is written beside the samples
    # table in the dataset's config/ folder
    if not ds.params.get("metadata"):
        ds.add_param("metadata", ds.params["input"].rsplit("/", 1)[0] + "/pipeline_metadata.tsv")

    # Overwrite the files Cirro wrote to the dataset's config/ folder (mapped in
    # process-input.json). Rewriting also drops their CRLF line endings, which the
    # workflow's shell steps would keep as part of the last column.
    samples.to_csv(ds.params["input"], index=None, sep="\t")
    metadata.to_csv(ds.params["metadata"], index=None, sep="\t")

    # Force params.json to be written: the HealthOmics pre-process Lambda fails the run
    # when the file is absent, and the SDK writes it only when a parameter changes.
    ds.keep_params(list(ds.params.keys()))
