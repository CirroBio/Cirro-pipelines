#!/usr/bin/env python3

import pandas as pd
from cirro.helpers.preprocess_dataset import PreprocessDataset

def make_manifest(dataset: PreprocessDataset) -> pd.DataFrame:

    # Log the parameters present
    for k, v in dataset.params.items():
        print(f"{k}: {v}")

    # PEPTIDE TABLE

    # Get the `peptide_table`, which was used to select the reference library
    peptide_table = dataset.params.get("peptide_table")

    assert peptide_table is not None, "ERROR: parameter peptide_table was not defined"

    # Set the `public_epitopes_csv` parameter to use the same library
    dataset.add_param(
        "public_epitopes_csv",
        peptide_table.replace("library.csv", "public_epitopes.csv")
    )
    print(f"Added public_epitopes_csv: {dataset.params['public_epitopes_csv']}")


    # READ SAMPLE SHEET

    # Filter out any index files that may have been uploaded
    ds.files = ds.files.loc[
        ds.files.apply(
            lambda r: r.get('readType', 'R') == 'R',
            axis=1
        )
    ]

    # Read the list of files available in this dataset
    print(f"Number of files associated with the dataset: {dataset.files.shape[0]:,}")
    assert dataset.files.shape[0] > 0, "No files detected -- there may be an error with data ingest"

    # Read the annotations of samples available in this dataset
    print(f"Number of samples associated with the dataset: {dataset.samplesheet.shape[0]:,}")
    assert dataset.samplesheet.shape[0] > 0, "No samples detected -- there may be an error with data ingest"


    # ANNOTATE CONTROLS
    # There must be a column for `control_status` defined
    msg = "ERROR: Must tag at least one sample as control_status = beads_only (column not found)"
    assert 'control_status' in dataset.samplesheet.columns.values, msg

    # Assign a Series linking each sample to its 'control_status'
    control_status = dataset.samplesheet.set_index("sample")["control_status"]

    # At least one sample must be labeled 'beads_only'
    print("Sample labels upon ingest")
    print(control_status.value_counts())
    msg = "ERROR: Must tag at least one sample as control_status = beads_only (value not found)"
    assert "beads_only" in control_status.values, msg

    # Mark all of the other samples as 'empirical'
    control_status = control_status.apply(
        lambda v: v if v == 'beads_only' else 'empirical'
    )
    print("Sample labels prior to analysis")
    print(control_status.value_counts())

    # Filter out any index reads, if any were provided
    files = dataset.files.loc[
        dataset.files.reindex(columns=['readType'])['readType'].fillna("R") == "R"
    ]

    # FORMAT SAMPLESHEET
    # to conform to the expectations of phipflow, annotating files appropriately
    samplesheet = files.assign(
        control_status=lambda d: d["sample"].apply(control_status.get)
    ).assign(
        fastq_filepath=lambda d: d["file"].apply(
            lambda uri: uri.rsplit("/", 1)[-1]
        ),
        sample_name=lambda d: d["fastq_filepath"].apply(
            lambda s: s.replace(".fastq.gz", "")
        )
    ).reindex(
        columns=["sample_name", "sample", "control_status", "fastq_filepath"]
    )

    return samplesheet


if __name__ == '__main__':
    ds = PreprocessDataset.from_running()

    samplesheet = make_manifest(ds)

    # Save the manifest
    samplesheet.to_csv("samplesheet.csv", index=None)
    ds.add_param('sample_table', 'samplesheet.csv')
