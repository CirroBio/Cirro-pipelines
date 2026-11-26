#!/usr/bin/env python3


def check_params(params: dict):
    if not params.get("formula"):
        raise ValueError("Must provide formula indicating columns to test")
    if not 0 <= params["fdr_cutoff"] < 1:
        raise ValueError("fdr_cutoff must be >= 0 and < 1")


if __name__ == "__main__":

    # Imported here so that test_preprocess.py runs without the SDK
    from cirro.helpers.preprocess_dataset import PreprocessDataset

    ds = PreprocessDataset.from_running()

    check_params(ds.params)

    # Force params.json to be written: the HealthOmics pre-process Lambda fails the run
    # when the file is absent, and the SDK writes it only when a parameter changes.
    ds.keep_params(list(ds.params.keys()))
