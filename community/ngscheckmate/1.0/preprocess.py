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
# sample, fastq_1, fastq_1
samplesheet = ds.files.pivot(
    index=["sampleIndex", "sample"],
    columns="read",
    values="file"
).rename(
    columns=lambda i: f"fastq_{int(i)}"
).reset_index(
).reindex(
    columns=["sample", "fastq_1", "fastq_2"]
)

ds.logger.info("Formatted samplesheet:")
ds.logger.info(samplesheet.to_csv(index=None))
assert samplesheet.shape[0] > 0, "No files detected -- there may be an error with data ingest"

# Drop any rows with missing data
total_n = samplesheet.shape[0]
samplesheet = samplesheet.dropna()
if samplesheet.shape[0] < total_n:
    ds.logger.info(f"After dropping lines with missing data, {samplesheet.shape[0]:,} remain")
    assert samplesheet.shape[0] > 0, "No files detected -- there may be an error with data ingest"

# Only keep a single file pair per sample
if samplesheet['sample'].value_counts().max() > 1:
    ds.logger.info("Only keeping a single pair of files per sample")
    samplesheet = samplesheet.groupby('sample').head(1)
    ds.logger.info(samplesheet.to_csv(index=None))

# Write to the dataset's config/ folder (mapped in process-input.json)
samplesheet.to_csv(ds.params["paired_fastq_samplesheet"], index=None)

# log
ds.logger.info(ds.params)

# Force params.json to be written: the HealthOmics pre-process Lambda fails the run
# when the file is absent, and the SDK writes it only when a parameter changes.
ds.keep_params(list(ds.params.keys()))
