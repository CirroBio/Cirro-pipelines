#!/usr/bin/env python3

from cirro.helpers.preprocess_dataset import PreprocessDataset
import pandas as pd

# Instantiate the Cirro dataset object
ds = PreprocessDataset.from_running()

# Set up the GEX reference
ds.add_param(
    "transcriptome_dir",
    f"{ds.references_base}/cellranger/" + {
        "Homo sapiens (GRCh38-2024)": "refdata-gex-GRCh38-2024-A",
        "Homo sapiens (GRCh38-2020)": "refdata-gex-GRCh38-2020-A",
        "Mus musculus (GRCm39-2024)": "refdata-gex-GRCm39-2024-A",
        "Mus musculus (mm10-2020)": "refdata-gex-mm10-2020-A"
    }[
        ds.params["reference"]
    ] + "/"   # a reference directory, which HealthOmics reads as an object without it
)

# If the user did not select a custom probe set
if ds.params.get("probe_set") is None or ds.params.get("probe_set") == "":

    # Use the default probes for the genome
    ds.add_param(
        "probe_set",
        f"{ds.references_base}/cellranger/flex/" + (
            "Chromium_Human_Transcriptome_Probe_Set_v1.0.1_GRCh38-2020-A.csv"
            if ds.params["reference"].startswith("Homo sapiens")
            else "Chromium_Mouse_Transcriptome_Probe_Set_v1.0.1_mm10-2020-A.csv" # noqa
        ),
        overwrite=True
    )

# Get the sample names used for each barcode. The form groups these under
# samples, but process-input.json maps each to a flat bcNNN key, and that is the
# name ds.params carries -- ds.params["samples"] is the samplesheet this script
# writes. A barcode the user left blank names no sample.
barcode_params = sorted(
    param for param in ds.params if param.startswith("bc0")
)
probe_barcodes = pd.DataFrame([
    dict(
        sample_id=ds.params[barcode],
        barcode=barcode.upper()
    )
    for barcode in barcode_params
    if ds.params[barcode]
])
msg = "User must specify at least one sample barcode used"
assert probe_barcodes.shape[0] > 0, msg

# If multiple Probe Barcodes were used for a sample,
# separate IDs with a pipe (e.g., BC001|BC002)
probe_barcodes = pd.DataFrame(dict(
    probe_barcode_ids=probe_barcodes.groupby(
        'sample_id'
    ).apply(
        lambda d: '|'.join(d['barcode'].tolist())
    )
)).reset_index()

# Save the sample barcode spreadsheet to the dataset's config/ folder
# (mapped in process-input.json)
ds.logger.info("Sample probe barcodes specified:")
ds.logger.info(probe_barcodes.to_csv(index=None))
probe_barcodes.to_csv(ds.params["probe_barcodes"], index=None)

ds.logger.info("Samples provided by the user:")
ds.logger.info(ds.samplesheet)
assert ds.samplesheet.shape[0] > 0, "No files detected -- there may be an error with data ingest"

# Write to the dataset's config/ folder (mapped in process-input.json)
ds.logger.info(f"Writing out {ds.samplesheet.shape[0]:,} lines to {ds.params['samples']}")
ds.samplesheet.to_csv(ds.params["samples"], index=None)

# Log the parameters present
for k, v in ds.params.items():
    ds.logger.info(f"{k}: {v}")

# Consumed above; not parameters the workflow declares.
for param in ["reference"] + barcode_params:
    ds.remove_param(param, force=True)
