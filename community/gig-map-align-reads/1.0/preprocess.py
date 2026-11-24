
from cirro.helpers.preprocess_dataset import PreprocessDataset

# Instantiate the Cirro dataset object
ds = PreprocessDataset.from_running()

# Filter out any index files that may have been uploaded
ds.files = ds.files.loc[
    ds.files.apply(
        lambda r: r.get('readType', 'R') == 'R',
        axis=1
    )
]

# Make a wide samplesheet
samplesheet = ds.wide_samplesheet(
    index=["sample", "sampleIndex", "dataset"],
    columns="read",
    values="file",
    column_prefix="R"
).reindex(
    columns=["sample", "R1", "R2"]
)
assert samplesheet.shape[0] > 0, "No files detected -- there may be an error with data ingest"

ds.logger.info("Samplesheet:")
ds.logger.info(samplesheet.to_csv(index=None))

# Write to the dataset's config/ folder (mapped in process-input.json)
ds.logger.info(f"Writing out {samplesheet.shape[0]:,} lines to {ds.params['samplesheet']}")
samplesheet.to_csv(
    ds.params["samplesheet"],
    index=None
)

# Log the parameters present
for k, v in ds.params.items():
    ds.logger.info(f"{k}: {v}")

# Force params.json to be written: the HealthOmics pre-process Lambda fails the run
# when the file is absent, and the SDK writes it only when a parameter changes.
ds.keep_params(list(ds.params.keys()))
