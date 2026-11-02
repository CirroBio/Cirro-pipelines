#!/usr/bin/env python3

from cirro.helpers.preprocess_dataset import PreprocessDataset
from io import StringIO
import pandas as pd

# Instantiate the Cirro dataset object
ds = PreprocessDataset.from_running()

# Datasets launched before the form offered a choice carry no flex_version
flex_version = ds.params.get("flex_version") or "v1"
ds.logger.info(f"Flex chemistry: {flex_version}")

# Flex v2 chemistries were added in CellRanger 10.0
if flex_version == "v2":
    msg = f"Flex v2 requires CellRanger 10.0 or later (selected: {ds.params['cellranger_version']})"
    assert int(ds.params["cellranger_version"].split(".")[0]) >= 10, msg

# CellRanger added emptydrops-minimum-umis in 7.1 and global-minimum-umis in 8.0.
# A threshold left blank is dropped so the CellRanger default applies.
cellranger_version = tuple(int(n) for n in ds.params["cellranger_version"].split(".")[:2])
for param, min_version in [("emptydrops_minimum_umis", (7, 1)), ("global_minimum_umis", (8, 0))]:
    if ds.params.get(param) is None:
        ds.remove_param(param, force=True)
    else:
        msg = f"{param} requires CellRanger {min_version[0]}.{min_version[1]} or later (selected: {ds.params['cellranger_version']})"
        assert cellranger_version >= min_version, msg

# The GEX reference and the 10X probe sets built against it, by Flex version.
# CellRanger refuses a probe set whose reference version differs from the GEX
# reference, and the v2 probe sets exist for the 2024-A references only.
transcriptome, probe_sets = {
    "Homo sapiens (GRCh38-2024)": (
        "refdata-gex-GRCh38-2024-A",
        {
            "v1": "Chromium_Human_Transcriptome_Probe_Set_v1.1.0_GRCh38-2024-A.csv",
            "v2": "Chromium_Human_Transcriptome_Probe_Set_v2.0.0_GRCh38-2024-A.csv"
        }
    ),
    "Homo sapiens (GRCh38-2020)": (
        "refdata-gex-GRCh38-2020-A",
        {"v1": "Chromium_Human_Transcriptome_Probe_Set_v1.0.1_GRCh38-2020-A.csv"}
    ),
    "Mus musculus (GRCm39-2024)": (
        "refdata-gex-GRCm39-2024-A",
        {
            "v1": "Chromium_Mouse_Transcriptome_Probe_Set_v1.1.0_GRCm39-2024-A.csv",
            "v2": "Chromium_Mouse_Transcriptome_Probe_Set_v2.0.0_GRCm39-2024-A.csv"
        }
    ),
    "Mus musculus (mm10-2020)": (
        "refdata-gex-mm10-2020-A",
        {"v1": "Chromium_Mouse_Transcriptome_Probe_Set_v1.0.1_mm10-2020-A.csv"}
    )
}[ds.params["reference"]]
msg = f"No Flex {flex_version} probe set is available for {ds.params['reference']}"
assert flex_version in probe_sets, msg
default_probe_set = probe_sets[flex_version]

ds.add_param(
    "transcriptome_dir",
    # a reference directory, which HealthOmics reads as an object without it
    f"{ds.references_base}/cellranger/{transcriptome}/"
)

# If the user did not select a custom probe set
if ds.params.get("probe_set") is None or ds.params.get("probe_set") == "":
    ds.add_param(
        "probe_set",
        f"{ds.references_base}/cellranger/flex/{default_probe_set}",
        overwrite=True
    )

# Get the sample names used for each barcode. The form groups these under
# samples, but process-input.json maps each to a flat bcNNN key, and that is the
# name ds.params carries -- ds.params["samples"] is the samplesheet this script
# writes. A barcode the user left blank names no sample.
barcode_params = sorted(
    param for param in ds.params if param.startswith("bc0")
)

# Datasets launched before the form offered a choice carry no barcode_source
barcode_source = ds.params.get("barcode_source") or "manual"
ds.logger.info(f"Probe barcode source: {barcode_source}")

if barcode_source == "manual":
    probe_barcodes = pd.DataFrame([
        dict(
            sample_id=ds.params[barcode],
            probe_barcode_ids=barcode.upper()
        )
        for barcode in barcode_params
        if ds.params.get(barcode)
    ])
else:
    # A table pasted into the form, or a Spreadsheet (CSV) reference
    if barcode_source == "table":
        probe_barcodes = pd.read_csv(StringIO(ds.params["barcode_table"]), dtype=str, skipinitialspace=True)
    else:
        probe_barcodes = pd.read_csv(ds.params["barcode_reference"], dtype=str, skipinitialspace=True)

    # CellRanger 10.1 also names this column sample_barcode_ids; 10.0 accepts only
    # probe_barcode_ids, and so does the workflow's config builder
    probe_barcodes = probe_barcodes.rename(columns={"sample_barcode_ids": "probe_barcode_ids"})

    for cname in ["sample_id", "probe_barcode_ids"]:
        msg = f"Barcode table must have a column named {cname} (found: {', '.join(probe_barcodes.columns)})"
        assert cname in probe_barcodes.columns, msg
    probe_barcodes = probe_barcodes[["sample_id", "probe_barcode_ids"]].dropna()

msg = "User must specify at least one sample barcode used"
assert probe_barcodes.shape[0] > 0, msg

# If multiple Probe Barcodes were used for a sample,
# separate IDs with a pipe (e.g., BC001|BC002)
probe_barcodes = pd.DataFrame(dict(
    probe_barcode_ids=probe_barcodes.groupby(
        'sample_id'
    ).apply(
        lambda d: '|'.join(d['probe_barcode_ids'].tolist())
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
for param in ["flex_version", "reference", "barcode_source", "barcode_table", "barcode_reference"] + barcode_params:
    ds.remove_param(param, force=True)
