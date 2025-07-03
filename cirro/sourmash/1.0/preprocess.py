#!/usr/bin/env python3

import json
from cirro.helpers.preprocess_dataset import PreprocessDataset

ds = PreprocessDataset.from_running()

# Write out the list of files
ds.logger.info(ds.files.to_csv(index=None))
ds.files.to_csv("samplesheet.csv", index=None)
ds.add_param("samplesheet", "samplesheet.csv")

# Set up the database based on the size and k selected by the user
prefix = "s3://pubweb-references/"

db_size = ds.params["db_size"]

db_map = {
    "2M: NCBI GenBank (2022) - bacteria, viruses, archaea, protozoa, and fungi": "genbank/genbank-2022.03-*-k",
    "1.1M: NCBI GenBank (2022) - bacteria": "genbank/genbank-2022.03-bacteria-k",
    "48k: NCBI GenBank (2022) - viruses": "genbank/genbank-2022.03-viral-k",
    "9k: NCBI GenBank (2022) - archaea": "genbank/genbank-2022.03-archaea-k",
    "1k: NCBI GenBank (2022) - protozoa": "genbank/genbank-2022.03-protozoa-k",
    "10k: NCBI GenBank (2022) - fungi": "genbank/genbank-2022.03-fungi-k",
    "85k: GTDB R08-RS214 bacterial genomic representatives": "sourmash/gtdb-rs214-reps.k",
    "403k: GTDB R08-RS214 all bacterial genomes": "sourmash/gtdb-rs214-k",
}
if db_size not in db_map:
    raise ValueError(f"Unrecognized option: {db_size}")

if "GTDB R08-RS214" in db_size:
    ds.add_param("tax_db", "s3://pubweb-references/sourmash/gtdb-rs214.tax.db")
else:
    ds.add_param("tax_db", "s3://pubweb-references/genbank/genbank-2022.03-tax.db")

ds.add_param(
    "db",
    f"{prefix}{db_map[db_size]}{ds.params['ksize']}.zip"
)

# Make sure that the user does not select a negative threshold value
msg = f"Minimum threshold cannot be negative ({ds.params['threshold_bp']})"
assert ds.params["threshold_bp"] > 0

ds.logger.info(json.dumps(ds.params, indent=4))