#!/usr/bin/env python3

import json
from cirro.helpers.preprocess_dataset import PreprocessDataset

ds = PreprocessDataset.from_running()

# Write out the list of files
ds.logger.info(ds.files.to_csv(index=None))
ds.files.to_csv("samplesheet.csv", index=None)
ds.add_param("samplesheet", "samplesheet.csv")

# Map the user-selected databases to their corresponding paths
db_map = {
    "1.1M: NCBI GenBank (2022) - bacteria": "genbank/genbank-2022.03-bacteria-k",
    "48k: NCBI GenBank (2022) - viruses": "genbank/genbank-2022.03-viral-k",
    "9k: NCBI GenBank (2022) - archaea": "genbank/genbank-2022.03-archaea-k",
    "1k: NCBI GenBank (2022) - protozoa": "genbank/genbank-2022.03-protozoa-k",
    "10k: NCBI GenBank (2022) - fungi": "genbank/genbank-2022.03-fungi-k",
    "85k: GTDB R08-RS214 bacterial genomic representatives": "sourmash/gtdb-rs214-reps.k",
    "403k: GTDB R08-RS214 all bacterial genomes": "sourmash/gtdb-rs214-k",
}


def format_db(
    bact_db: str,
    viral_db: str,
    archaea_db: str,
    protozoa_db: str,
    fungi_db: str,
    ksize: int,
    prefix: str = "s3://pubweb-references/",
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


ds.add_param(
    "db",
    format_db(**ds.params)
)

# Make sure that the user does not select a negative threshold value
msg = f"Minimum threshold cannot be negative ({ds.params['threshold_bp']})"
assert ds.params["threshold_bp"] > 0

ds.logger.info(json.dumps(ds.params, indent=4))
