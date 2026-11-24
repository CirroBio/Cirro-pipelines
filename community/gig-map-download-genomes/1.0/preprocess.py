def genome_table_globs(genome_tables: str) -> dict:
    """The CSV and TSV tables of NCBI genomes in the picked folders, as gig-map's
    comma-separated glob params.

    Built here rather than with a `|/...` suffix in process-input.json: a whole dataset
    picked in the folder picker arrives as its data/ folder with a trailing slash.
    """
    folders = [folder.rstrip("/") for folder in genome_tables.split(",")]
    return {
        f"genome_{ext}": ",".join(f"{folder}/*.{ext}" for folder in folders)
        for ext in ["csv", "tsv"]
    }


if __name__ == "__main__":

    # Imported here so that test_preprocess.py runs without the SDK
    from cirro.helpers.preprocess_dataset import PreprocessDataset

    ds = PreprocessDataset.from_running()

    for param, value in genome_table_globs(ds.params["genome_tables"]).items():
        ds.add_param(param, value)
    ds.remove_param("genome_tables")
