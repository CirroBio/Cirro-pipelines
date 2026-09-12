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

resolve_references(ds, "transcriptome_dir")

# Make a wide samplesheet with the columns
# sample, fastq_1, fastq_1
samplesheet = (
    ds.files
    .assign(
        readType=lambda d: d.apply(
            lambda r: r.get("readType", "R"),
            axis=1
        ),
        lane=lambda d: d.apply(
            lambda r: r.get("lane", 1),
            axis=1
        ),
    )
    .query("readType == 'R'")
    # A row with no read number cannot be paired, and would pivot
    # into a column the rename below cannot name.
    .loc[lambda d: d["read"].notna()]
    .pivot(
        index=["sampleIndex", "lane", "sample", "dataset"],
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

# Force params.json to be written: the HealthOmics pre-process Lambda fails the run
# when the file is absent, and the SDK writes it only when a parameter changes.
ds.keep_params(list(ds.params.keys()))
