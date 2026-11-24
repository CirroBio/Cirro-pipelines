
from cirro.helpers.preprocess_dataset import PreprocessDataset

# Instantiate the Cirro dataset object
ds = PreprocessDataset.from_running()

# Filter out any index files that may have been uploaded
ds.files = ds.files.loc[
    ds.files.apply(
        lambda r: r.get('readType', 'R') != 'I',
        axis=1
    )
]

# Make a wide samplesheet
samplesheet = (
    ds
    .files
    .reindex(columns=["sample", "sampleIndex", "dataset", "lane", "read", "file"])
    # A row with no read number cannot be paired, and would pivot
    # into a column the rename below cannot name.
    .loc[lambda d: d["read"].notna()]
    .pivot(
        index=["sampleIndex", "sample", "dataset", "lane"],
        columns="read",
        values="file"
    )
    .rename(columns=lambda i: f"R{int(i)}")
    .reset_index()
    .reindex(columns=["sample", "R1", "R2"])
)
assert samplesheet.shape[0] > 0, "No files detected -- there may be an error with data ingest"

ds.logger.info("Samplesheet:")
ds.logger.info(samplesheet.to_csv(index=None))

# Write out the sample sheet
ds.logger.info(f"Writing out {samplesheet.shape[0]:,} lines to samplesheet.csv")
# Write to the dataset's config/ folder (mapped in process-input.json)
samplesheet.to_csv(
    ds.params["samplesheet"],
    index=None
)

# Provide the metadata table
ds.samplesheet.to_csv(ds.params["metadata"], index=None)

# Log the parameters present
for k, v in ds.params.items():
    ds.logger.info(f"{k}: {v}")

# Force params.json to be written: the HealthOmics pre-process Lambda fails the run
# when the file is absent, and the SDK writes it only when a parameter changes.
ds.keep_params(list(ds.params.keys()))
