#!/usr/bin/env python3

import pandas as pd
from io import StringIO
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

# The user must provide columns for `grouping` and `feature_types`
for cname in ['grouping', 'feature_types']:
    msg = f"The user must annotate the '{cname}' for each sample"
    assert cname in ds.samplesheet.columns.values, msg

groupings = ds.samplesheet[["sample", "grouping", "feature_types"]]
ds.logger.info("Sample sheet provided by the user:")
ds.logger.info(groupings)
assert groupings.shape[0] > 0, "No files detected -- there may be an error with data ingest"

# Write to the dataset's config/ folder (mapped in process-input.json)
ds.logger.info(f"Writing out {groupings.shape[0]:,} lines to {ds.params['grouping']}")
groupings.to_csv(ds.params["grouping"], index=None)

# Build fastq_dir as a comma-delimited list of all input dataset paths
data_paths = [dataset['dataPath'] for dataset in ds.metadata['inputs']]
assert len(data_paths) > 0, "No input datasets found"
ds.add_param("fastq_dir", ",".join(data_paths))

# If either the feature_csv was not provided
for kw in ["feature_csv"]:

    # If the user did not provide the keyword
    if kw in ds.params and ds.params[kw] is None:

        # Remove it from the dict (so that the workflow default is used)
        ds.remove_param(kw)

# Set below, only when the user supplies a probe barcode table
probe_barcodes = None

# If the user indicated that this is fixed RNA profiling
if ds.params.get("is_frp"):
    ds.logger.info("User indicated that this is fixed RNA profiling")

    # Add the appropriate probe set
    if "GRCh38" in ds.params["transcriptome_dir"]:
        ds.logger.info("Adding human reference probe set")
        ds.add_param(
            "probes_csv",
            f"{ds.references_base}/cellranger/flex/Chromium_Human_Transcriptome_Probe_Set_v1.0.1_GRCh38-2020-A.csv"
        )
    else:
        ds.logger.info("Adding mouse reference probe set")
        ds.add_param(
            "probes_csv",
            f"{ds.references_base}/cellranger/flex/Chromium_Mouse_Transcriptome_Probe_Set_v1.0.1_mm10-2020-A.csv"
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
            probe_barcodes = None
        else:
            ds.logger.info(f"Detected {probe_barcodes.shape[0]:,} samples in the FRP samples table")
            # Write to the dataset's config/ folder (mapped in process-input.json)
            probe_barcodes.to_csv(ds.params["probe_barcodes"], index=None)

# No table was written, so the mapped location holds no object
if probe_barcodes is None:
    ds.remove_param("probe_barcodes")

# Log the parameters present
for k, v in ds.params.items():
    ds.logger.info(f"{k}: {v}")

# Consumed above; not parameters the workflow declares.
for param in ("frp_samples", "is_frp"):
    ds.remove_param(param, force=True)
