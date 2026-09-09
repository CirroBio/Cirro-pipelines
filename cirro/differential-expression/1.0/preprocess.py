#!/usr/bin/env python3
from cirro.helpers.preprocess_dataset import PreprocessDataset


ds = PreprocessDataset.from_running()
samplesheet = ds.samplesheet

ds.logger.info(f"Read in samplesheet with {samplesheet.shape[0]:,} rows and {samplesheet.shape[1]:,} columns")
assert samplesheet.shape[0] > 0, "No files detected -- there may be an error with data ingest"

ds.logger.info(samplesheet.head())

# Write to the dataset's config/ folder (mapped in process-input.json)
samplesheet.to_csv(ds.params["manifest"], index=None)
