#!/usr/bin/env python3

from cirro.helpers.preprocess_dataset import PreprocessDataset


def resolve_references(ds: PreprocessDataset, *params: str):
    """Make reference params absolute against the references bucket.

    Form values name a path under that bucket rather than a full URI, so the
    configuration is not tied to one deployment. A value that is already absolute
    is left alone, which keeps datasets created before that change re-runnable.
    """
    for param in params:
        value = ds.params.get(param)
        if not isinstance(value, str) or not value or value.startswith("s3://"):
            continue
        ds.add_param(param, f"{ds.references_base}/{value}", overwrite=True)

ds = PreprocessDataset.from_running()

resolve_references(ds, "db")

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
# sample, fastq_1, fastq_1
samplesheet = (
    ds.files
    .reindex(columns=["dataset", "sampleIndex", "sample", "lane", "read", "file"])
    .pivot(
        index=["dataset", "sampleIndex", "sample", "lane"],
        columns="read",
        values="file"
    )
    .rename(columns=lambda i: f"fastq_{int(i)}")
    .reset_index()
    .reindex(columns=["sample", "fastq_1", "fastq_2"])
)

ds.logger.info("Formatted samplesheet:")
ds.logger.info(samplesheet.to_csv(index=None))
assert samplesheet.shape[0] > 0, "No files detected -- there may be an error with data ingest"

# Write to the dataset's config/ folder (mapped in process-input.json)
samplesheet.to_csv(ds.params["samplesheet"], index=None)

# log
ds.logger.info(ds.params)
