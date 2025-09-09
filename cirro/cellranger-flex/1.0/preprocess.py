#!/usr/bin/env python3

from cirro.helpers.preprocess_dataset import PreprocessDataset
import pandas as pd

# Instantiate the Cirro dataset object
ds = PreprocessDataset.from_running()

# Set up the GEX reference
ds.add_param(
    "transcriptome_dir",
    "s3://pubweb-references/cellranger/" + {
        "Homo sapiens (GRCh38-2024)": "refdata-gex-GRCh38-2024-A",
        "Homo sapiens (GRCh38-2020)": "refdata-gex-GRCh38-2020-A",
        "Mus musculus (GRCm39-2024)": "refdata-gex-GRCm39-2024-A",
        "Mus musculus (mm10-2020)": "refdata-gex-mm10-2020-A"
    }[
        ds.params["reference"]
    ]
)

# If the user did not select a custom probe set
if ds.params.get("probe_set") is None or ds.params.get("probe_set") == "":

    # Use the default probes for the genome
    ds.add_param(
        "probe_set",
        "s3://pubweb-references/cellranger/flex/" + (
            "Chromium_Human_Transcriptome_Probe_Set_v1.0.1_GRCh38-2020-A.csv"
            if ds.params["reference"].startswith("Homo sapiens")
            else "Chromium_Mouse_Transcriptome_Probe_Set_v1.0.1_mm10-2020-A.csv" # noqa
        ),
        overwrite=True
    )

# Get the sample names used for each barcode
probe_barcodes = pd.DataFrame([
    dict(
        sample_id=sample,
        barcode=barcode.upper()
    )
    for barcode, sample in ds.params.items()
    if barcode.startswith("bc0")
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

# Save the sample barcode spreadsheet
ds.logger.info("Sample probe barcodes specified:")
ds.logger.info(probe_barcodes.to_csv(index=None))
probe_barcodes.to_csv("probe_barcodes.csv", index=None)
ds.add_param("probe_barcodes", "probe_barcodes.csv")

ds.logger.info("Samples provided by the user:")
ds.logger.info(ds.samplesheet)
assert ds.samplesheet.shape[0] > 0, "No files detected -- there may be an error with data ingest"

# Write out the sample sheet
ds.logger.info(f"Writing out {ds.samplesheet.shape[0]:,} lines to samples.csv")
ds.samplesheet.to_csv("samples.csv", index=None)

# Add it to the params
ds.add_param(
    "samples",
    "samples.csv"
)

# Log the parameters present
for k, v in ds.params.items():
    ds.logger.info(f"{k}: {v}")
