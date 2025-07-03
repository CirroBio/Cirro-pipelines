#!/usr/bin/env python3

import json
from cirro.helpers.preprocess_dataset import PreprocessDataset

ds = PreprocessDataset.from_running()

# Write out the list of files
ds.logger.info(ds.files.to_csv(index=None))
ds.files.to_csv("samplesheet.csv", index=None)
ds.add_param("samplesheet", "samplesheet.csv")

# Set up the database based on the size and k selected by the user
prefix = "s3://pubweb-references/sourmash/"
if ds.params["db_size"] == "85k: GTDB R08-RS214 genomic representatives":
    ds.add_param(
        "db",
        f"{prefix}gtdb-rs214-reps.k{ds.params['ksize']}.zip"
    )
else:
    msg = f"Unrecognized option: {ds.params['db_size']}"
    assert ds.params["db_size"] == "403k: GTDB R08-RS214 all genomes", msg
    ds.add_param(
        "db",
        f"{prefix}gtdb-rs214-k{ds.params['ksize']}.zip"
    )

# Make sure that the user does not select a negative threshold value
msg = f"Minimum threshold cannot be negative ({ds.params['threshold_bp']})"
assert ds.params["threshold_bp"] > 0

ds.logger.info(json.dumps(ds.params, indent=4))
