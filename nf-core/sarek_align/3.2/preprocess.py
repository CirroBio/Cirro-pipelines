#!/usr/bin/env python3

from cirro.helpers.preprocess_dataset import PreprocessDataset
import os
import pandas as pd


def make_manifest(ds: PreprocessDataset) -> pd.DataFrame:

    # Filter out any index files that may have been uploaded
    ds.files = ds.files.loc[
        ds.files.apply(
            lambda r: r.get('readType', 'R') == 'R',
            axis=1
        )
    ]

    # Make a wide manifest
    manifest: pd.DataFrame = ds.wide_samplesheet(
        index=["sampleIndex", "sample", "lane", "dataset"],
        columns="read",
        values="file",
        column_prefix="fastq_"
    )
    assert manifest.shape[0] > 0, "No files detected -- error with data ingest"

    # Get the sample metadata (if any)
    # Populate the 'patient' column with the provided value,
    # falling back to the sample ID if missing
    samples = (
        ds.samplesheet
        .reindex(columns=["sample", "patient"])
        .assign(patient=lambda d: d['patient'].fillna(d['sample']))
        .set_index("sample")
    )

    # 1. Use the 'patient' column if provided, falling back to the 'sample'
    # 2. Order the columns
    # 3. Overwrite the 'lane' column to provide a unique value per-row
    # (This is necessary to account for datasets which merge flowcells)
    manifest = (
        manifest
        .set_index("sample")
        .assign(patient=samples["patient"])
        .reset_index()
        .reindex(columns=['patient', 'sample', 'lane', 'fastq_1', 'fastq_2'])
        .assign(lane=[str(i)for i in range(manifest.shape[0])])
    )

    return manifest


if __name__ == "__main__":

    # Load the information for this dataset
    ds = PreprocessDataset.from_running()

    # Make the samplesheet
    manifest = make_manifest(ds)

    # Log the manifest
    ds.logger.info(manifest.to_csv(index=None))

    # Write manifest
    manifest.to_csv("manifest.csv", index=None)

    # Dynamic Resource Usage
    # `compute_multiplier` == 2 for WGS and 1 for WES
    ds.add_param(
        "compute_multiplier",
        int(2 - int(ds.params.get('wes')))
    )

    # If an intervals file was not selected, use --no_intervals
    if not ds.params.get("intervals"):
        ds.add_param("no_intervals", True)

    # log all params
    ds.logger.info(ds.params)
