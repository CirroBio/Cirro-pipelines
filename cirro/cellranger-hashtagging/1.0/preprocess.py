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

# Instantiate the Cirro dataset object
ds = PreprocessDataset.from_running()

resolve_references(ds, "transcriptome_dir", "vdj_dir")

# These name CellRanger reference directories. HealthOmics reads an S3 value with no
# trailing separator as an object key and rejects the run, so mark them as folders.
for ref_dir in ("transcriptome_dir", "vdj_dir"):
    ref_path = ds.params.get(ref_dir)
    if ref_path:
        ds.add_param(ref_dir, ref_path.rstrip("/") + "/", overwrite=True)

# The user must provide a column for `feature_types`
for cname in ['feature_types']:
    msg = f"The user must annotate the '{cname}' for each sample"
    assert cname in ds.samplesheet.columns.values, msg

# If the `grouping` column is provided, remove it
if "grouping" in ds.samplesheet.columns.values:
    ds.logger.info("Removing the 'grouping' column from the sample sheet")
    ds.samplesheet.drop(columns=["grouping"], inplace=True)

# Build fastq_dir as a comma-delimited list of all input dataset paths
data_paths = [dataset['dataPath'] for dataset in ds.metadata['inputs']]
assert len(data_paths) > 0, "No input datasets found"
ds.add_param("fastq_dir", ",".join(data_paths))
ds.logger.info(f"fastq_dir: {ds.params['fastq_dir']}")

ds.logger.info("Sample sheet provided by the user:")
ds.logger.info(ds.samplesheet)
assert ds.samplesheet.shape[0] > 0, "No files detected -- there may be an error with data ingest"

# Write to the dataset's config/ folder (mapped in process-input.json)
ds.logger.info(f"Writing out {ds.samplesheet.shape[0]:,} lines to {ds.params['grouping']}")
ds.samplesheet.to_csv(ds.params["grouping"], index=None)

# If the feature_csv was not provided
if "feature_csv" in ds.params and ds.params["feature_csv"] is None:

    # Remove it from the dict (so that the workflow default is used)
    ds.remove_param("feature_csv")

# Log the parameters present
for k, v in ds.params.items():
    ds.logger.info(f"{k}: {v}")
