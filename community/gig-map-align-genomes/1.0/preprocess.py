def catalog_genes(catalog: str) -> str:
    """The deduplicated genes inside the picked gene catalog dataset.

    Built here rather than with a `|/...` suffix in process-input.json: a whole dataset
    picked in the folder picker arrives as its data/ folder with a trailing slash.
    """
    return f"{catalog.rstrip('/')}/centroids.faa.gz"


if __name__ == "__main__":

    # Imported here so that test_preprocess.py runs without the SDK
    from cirro.helpers.preprocess_dataset import PreprocessDataset

    ds = PreprocessDataset.from_running()

    ds.add_param("genes", catalog_genes(ds.params["genes"]), overwrite=True)

    # Force params.json to be written: the HealthOmics pre-process Lambda fails the run
    # when the file is absent, and the SDK writes it only when a parameter changes.
    ds.keep_params(list(ds.params.keys()))
