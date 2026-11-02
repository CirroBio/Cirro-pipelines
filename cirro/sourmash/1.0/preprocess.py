#!/usr/bin/env python3

import json

import pandas as pd

# Map the user-selected databases to their corresponding paths
db_map = {
    "1.1M: NCBI GenBank (2022) - bacteria": "genbank/genbank-2022.03-bacteria-k",
    "48k: NCBI GenBank (2022) - viruses": "genbank/genbank-2022.03-viral-k",
    "9k: NCBI GenBank (2022) - archaea": "genbank/genbank-2022.03-archaea-k",
    "1k: NCBI GenBank (2022) - protozoa": "genbank/genbank-2022.03-protozoa-k",
    "10k: NCBI GenBank (2022) - fungi": "genbank/genbank-2022.03-fungi-k",
    "85k: GTDB R08-RS214 bacterial representatives": "sourmash/gtdb-rs214-reps.k",
    "403k: GTDB R08-RS214 all bacteria": "sourmash/gtdb-rs214-k",
}


def read_samplesheet(path) -> pd.DataFrame:
    """Read the samplesheet filled out in the form, keeping every value as text."""

    # The form field is marked required, but an empty table can still be submitted,
    # in which case Cirro either omits the param or writes an empty file
    if not path:
        raise ValueError("No samplesheet was provided: add at least one row of sequence files")
    try:
        samplesheet = pd.read_csv(path, dtype=str, keep_default_na=False)
    except pd.errors.EmptyDataError:
        samplesheet = pd.DataFrame()
    if samplesheet.empty:
        raise ValueError("The samplesheet is empty: add at least one row of sequence files")
    return samplesheet


def format_db(
    bact_db: str,
    viral_db: str,
    archaea_db: str,
    protozoa_db: str,
    fungi_db: str,
    ksize: int,
    prefix: str,
    suffix: str = ".zip",
    **kwargs
) -> str:
    """Format the database path based on the selected database and ksize."""
    dbs = [
        db_map[db] + str(ksize)
        for db in [bact_db, viral_db, archaea_db, protozoa_db, fungi_db]
        if db in db_map
    ]
    return prefix + "{" + ",".join(dbs) + "}" + suffix


def check_threshold(threshold_bp: int):
    """sourmash gather accepts a threshold of 0, which keeps every match."""

    if threshold_bp < 0:
        raise ValueError(f"Minimum threshold cannot be negative ({threshold_bp})")


if __name__ == "__main__":

    # Imported here so that test_preprocess.py runs without the SDK
    from cirro.helpers.preprocess_dataset import PreprocessDataset

    ds = PreprocessDataset.from_running()

    samplesheet = read_samplesheet(ds.params.get("samplesheet"))
    ds.logger.info("Samplesheet:")
    ds.logger.info(samplesheet.to_csv(index=None))

    ds.add_param(
        "db",
        format_db(**ds.params, prefix=f"{ds.references_base}/")
    )

    check_threshold(ds.params["threshold_bp"])

    ds.logger.info(json.dumps(ds.params, indent=4))

    # Consumed above; not parameters the workflow declares.
    for param in ("archaea_db", "bact_db", "fungi_db", "protozoa_db", "viral_db"):
        ds.remove_param(param, force=True)
