#!/usr/bin/env python3

from cirro.helpers.preprocess_dataset import PreprocessDataset

ds = PreprocessDataset.from_running()

# Make a wide manifest
manifest = ds.wide_samplesheet(
    index=["sampleIndex", "sample", "lane"],
    columns="read",
    values="file",
    column_prefix="fastq_"
).sort_values(
    by="sample"
)
ds.logger.info(manifest)
assert manifest.shape[0] > 0, "No files detected -- there may be an error with data ingest"

# Cast comprehensive from int 0/1 to boolean
comprehensive = bool(ds.params.get('comprehensive'))
ds.add_param('comprehensive', comprehensive, overwrite=True)

manifest.to_csv("samplesheet.csv", index=None)
