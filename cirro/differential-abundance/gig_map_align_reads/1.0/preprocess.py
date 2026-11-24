#!/usr/bin/env python3

from cirro.helpers.preprocess_dataset import PreprocessDataset
import pandas as pd

ds = PreprocessDataset.from_running()


def log_table(title: str, df: pd.DataFrame, **kwargs):
    ds.logger.info(title)
    msg = "No rows found in table -- stopping"
    assert df.shape[0] > 0, msg
    for line in df.to_csv(**kwargs).split("\n"):
        ds.logger.info(line)


def make_metadata() -> pd.DataFrame:

    log_table("Sample metadata", ds.samplesheet, index=None)

    # Write to the dataset's config/ folder (mapped in process-input.json)
    ds.samplesheet.to_csv(ds.params["metadata"], index=None)


def main():
    # Get the input data
    make_metadata()

    msg = "Must provide formula indicating columns to test"
    assert len(ds.params.get("formula", "")) > 0, msg

    for kw in ["min_abund", "fdr_cutoff", "min_coverage"]:
        assert ds.params[kw] >= 0, f"{kw} must be >= 0"
        assert ds.params[kw] < 1, f"{kw} must be < 1"


if __name__ == "__main__":
    main()

    # Consumed above; not a parameter the workflow declares.
    ds.remove_param("min_abund", force=True)

    # Force params.json to be written: the HealthOmics pre-process Lambda fails the run
    # when the file is absent, and the SDK writes it only when a parameter changes.
    ds.keep_params(list(ds.params.keys()))
