#!/usr/bin/env python3

import pandas as pd
from io import StringIO
from cirro.helpers.preprocess_dataset import PreprocessDataset

# Instantiate the Cirro dataset object
ds = PreprocessDataset.from_running()

# The user must provide columns for `grouping` and `feature_types`
for cname in ['grouping', 'feature_types']:
    msg = f"The user must annotate the '{cname}' for each sample"
    assert cname in ds.samplesheet.columns.values, msg

groupings = ds.samplesheet[["sample", "grouping", "feature_types"]]
ds.logger.info("Sample sheet provided by the user:")
ds.logger.info(groupings)
assert groupings.shape[0] > 0, "No files detected -- there may be an error with data ingest"

# Write out the sample sheet
ds.logger.info(f"Writing out {groupings.shape[0]:,} lines to sample.grouping.csv")
groupings.to_csv("sample.grouping.csv", index=None)

# Add it to the params
ds.add_param(
    "grouping",
    "sample.grouping.csv"
)

# Find the folder which contains the FASTQ files (to account for subfolders)
ds.add_param(
    "fastq_dir",
    ds.files["file"][0].rsplit("/", 1)[0]
)

# If either the feature_csv was not provided
for kw in ["feature_csv"]:

    # If the user did not provide the keyword
    if kw in ds.params and ds.params[kw] is None:

        # Remove it from the dict (so that the workflow default is used)
        ds.remove_param(kw)

# If the user indicated that this is fixed RNA profiling
if ds.params.get("is_frp"):
    ds.logger.info("User indicated that this is fixed RNA profiling")

    # Add the appropriate probe set
    if "GRCh38" in ds.params["transcriptome_dir"]:
        ds.logger.info("Adding human reference probe set")
        ds.add_param(
            "probes_csv",
            "s3://pubweb-references/cellranger/flex/Chromium_Human_Transcriptome_Probe_Set_v1.0.1_GRCh38-2020-A.csv"
        )
    else:
        ds.logger.info("Adding mouse reference probe set")
        ds.add_param(
            "probes_csv",
            "s3://pubweb-references/cellranger/flex/Chromium_Mouse_Transcriptome_Probe_Set_v1.0.1_mm10-2020-A.csv"
        )

    # Parse the samples table provided by the user
    if ds.params.get("frp_samples") is None or len(ds.params["frp_samples"]) == 0:
        ds.logger.info("User did not provide a FRP samples table")
    else:
        ds.logger.info("Parsing the FRP samples table provided by the user")
        ds.logger.info(ds.params["frp_samples"])
        probe_barcodes = pd.read_table(StringIO(ds.params["frp_samples"]), sep=",")
        for line in probe_barcodes.to_csv(index=None).split("\n"):
            ds.logger.info(line)

        if probe_barcodes.shape[0] == 0:
            ds.logger.info("No samples detected in the FRP samples table")
        else:
            ds.logger.info(f"Detected {probe_barcodes.shape[0]:,} samples in the FRP samples table")
            probe_barcodes.to_csv("probe_barcodes.csv", index=None)
            ds.add_param("probe_barcodes", "probe_barcodes.csv")

# Log the parameters present
for k, v in ds.params.items():
    ds.logger.info(f"{k}: {v}")
