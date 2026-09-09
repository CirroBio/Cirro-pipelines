#!/usr/bin/env python3

from cirro.helpers.preprocess_dataset import PreprocessDataset

ds = PreprocessDataset.from_running()

ds.logger.info("Files annotated in the dataset:")
ds.logger.info(ds.files.to_csv(index=None))

# Filter out any index files that may have been uploaded
ds.files = ds.files.loc[
    ds.files.apply(
        lambda r: r.get('readType', 'R') == 'R',
        axis=1
    )
]

# Make a wide samplesheet with the columns
# sample, R1, R1
samplesheet = (
    ds.files
    .reindex(columns=["dataset", "sampleIndex", "sample", "lane", "read", "file"])
    .pivot(
        index=["dataset", "sampleIndex", "sample", "lane"],
        columns="read",
        values="file"
    )
    .rename(columns=lambda i: f"R{int(i)}")
    .reset_index()
    .rename(columns=dict(sample="specimen"))
    .reindex(columns=["specimen", "R1", "R2"])
)

ds.logger.info("Formatted samplesheet:")
ds.logger.info(samplesheet.to_csv(index=None))
assert samplesheet.shape[0] > 0, "No files detected -- there may be an error with data ingest"

# Write to the dataset's config/ folder (mapped in process-input.json)
samplesheet.to_csv(ds.params["sample_sheet"], index=None)

# log
ds.logger.info(ds.params)
