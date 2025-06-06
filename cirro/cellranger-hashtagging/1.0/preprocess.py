#!/usr/bin/env python3

from cirro.helpers.preprocess_dataset import PreprocessDataset

# Instantiate the Cirro dataset object
ds = PreprocessDataset.from_running()

# The user must provide a column for `feature_types`
for cname in ['feature_types']:
    msg = f"The user must annotate the '{cname}' for each sample"
    assert cname in ds.samplesheet.columns.values, msg

# If the `grouping` column is provided, remove it
if "grouping" in ds.samplesheet.columns.values:
    ds.logger.info("Removing the 'grouping' column from the sample sheet")
    ds.samplesheet.drop(columns=["grouping"], inplace=True)

ds.logger.info("Sample sheet provided by the user:")
ds.logger.info(ds.samplesheet)
assert ds.samplesheet.shape[0] > 0, "No files detected -- there may be an error with data ingest"

# Write out the sample sheet
ds.logger.info(f"Writing out {ds.samplesheet.shape[0]:,} lines to sample.grouping.csv")
ds.samplesheet.to_csv("sample.grouping.csv", index=None)

# Add it to the params
ds.add_param(
    "grouping",
    "sample.grouping.csv"
)

# If the feature_csv was not provided
if "feature_csv" in ds.params and ds.params["feature_csv"] is None:

    # Remove it from the dict (so that the workflow default is used)
    ds.remove_param("feature_csv")

# Log the parameters present
for k, v in ds.params.items():
    ds.logger.info(f"{k}: {v}")
