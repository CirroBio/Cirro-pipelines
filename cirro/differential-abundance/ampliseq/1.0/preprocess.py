#!/usr/bin/env python3


def check_params(params: dict):
    if not params.get("formula"):
        raise ValueError("Must provide formula indicating columns to test")
    if not 0 <= params["fdr_cutoff"] < 1:
        raise ValueError("fdr_cutoff must be >= 0 and < 1")


def ampliseq_tables(ampliseq_results: str) -> dict:
    """The counts and taxonomy tables inside the picked nf-core/ampliseq output.

    A whole dataset picked in the folder picker arrives as its data/ folder with a
    trailing slash, so the paths are joined here rather than with a `|/...` suffix.
    """
    folder = ampliseq_results.rstrip("/")
    return {
        "counts": f"{folder}/qiime2/abundance_tables/feature-table.tsv",
        "taxonomy": f"{folder}/qiime2/taxonomy/taxonomy.tsv",
    }


if __name__ == "__main__":

    # Imported here so that test_preprocess.py runs without the SDK
    from cirro.helpers.preprocess_dataset import PreprocessDataset

    ds = PreprocessDataset.from_running()

    check_params(ds.params)

    for param, value in ampliseq_tables(ds.params["ampliseq_results"]).items():
        ds.add_param(param, value, overwrite=True)

    # Consumed above; not a parameter the workflow declares
    ds.remove_param("ampliseq_results")
