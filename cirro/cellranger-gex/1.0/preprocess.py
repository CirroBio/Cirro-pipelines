from cirro.helpers.preprocess_dataset import PreprocessDataset

ds = PreprocessDataset.from_running()

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

# Write out to a file
samplesheet.to_csv("samplesheet.csv", index=None)

# Add the param for the samplesheet
ds.add_param("samplesheet", "samplesheet.csv")
