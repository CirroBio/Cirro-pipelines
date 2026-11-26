#!/usr/bin/env python3

import pandas as pd


def read_genome_table(path) -> pd.DataFrame:
    """Read the genome table filled out in the form, keeping every value as text."""

    # The form field is marked required, but an empty table can still be submitted,
    # in which case Cirro either omits the param or writes an empty file
    if not path:
        raise ValueError("No genomes were provided: add one row per genome")
    try:
        genomes = pd.read_csv(path, dtype=str, keep_default_na=False)
    except pd.errors.EmptyDataError:
        genomes = pd.DataFrame()
    if genomes.empty:
        raise ValueError("The genome table is empty: add one row per genome")
    return genomes


def format_manifest(genomes: pd.DataFrame) -> pd.DataFrame:
    """The genome manifest in the shape the workflow reads.

    A row with no GenBank FTP folder is read from its uri rather than downloaded.
    """
    return genomes.rename(columns={"organism": "#Organism Name"}).assign(**{"GenBank FTP": ""})


if __name__ == "__main__":

    # Imported here so that test_preprocess.py runs without the SDK
    from cirro.helpers.preprocess_dataset import PreprocessDataset

    ds = PreprocessDataset.from_running()

    genomes = read_genome_table(ds.params.get("genome_table"))
    # Write to the dataset's config/ folder (mapped in process-input.json)
    format_manifest(genomes).to_csv(ds.params["genomes"], index=None)

    # Consumed above; not a parameter the workflow declares.
    ds.remove_param("genome_table", force=True)
