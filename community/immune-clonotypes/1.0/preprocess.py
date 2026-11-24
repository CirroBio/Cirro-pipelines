#!/usr/bin/env python3

import logging
import pandas as pd
import os
from cirro.helpers.preprocess_dataset import PreprocessDataset

# Set up logging
logFormatter = logging.Formatter(
    '%(asctime)s %(levelname)-8s [process-prepare] %(message)s'
)
logger = logging.getLogger()
logger.setLevel(logging.INFO)

# Also write to STDOUT
consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
logger.addHandler(consoleHandler)

class Dataset:
    """
    Generate all of the information available for the dataset being executed:
    files: DataFrame of all files being input
    samples: DataFrame of all sample metadata associated with the files in this dataset
    """

    def __init__(self):

        # Read dataset information from the environment
        self.s3_dataset = os.getenv("PW_S3_DATASET")

        # Raise an error if the variable is not set
        assert self.s3_dataset is not None, "Environment variable not found: PW_S3_DATASET"

        # Read the file-level information
        self.files = self.read_csv(
            "config/files.csv",
            required_columns=["sample", "file"]
        )

        # Make sure that the expected columns are present
        assert "sample" in self.files.columns.values, "Did not find expected column in "

        # Read the sample-level information
        self.samplesheet = self.read_csv(
            "config/samplesheet.csv",
            required_columns=["sample"]
        )

    def read_csv(self, suffix: str, required_columns=[]):
        """Read a CSV from the dataset and check for any required columns."""

        df = pd.read_csv(f"{self.s3_dataset}/{suffix}")
        for col in required_columns:
            assert col in df.columns.values, f"Did not find expected columns {col} in {self.s3_dataset}/{suffix}"
        return df

    def log(self):
        """Print logging messages about the dataset."""

        print(f"Storage location for dataset: {self.s3_dataset}")
        print(f"Number of files in dataset: {self.files.shape[0]:,}")
        print(f"Number of samples in dataset: {self.samplesheet.shape[0]:,}")


# Get the dataset information
dataset = Dataset()

# Log information about the dataset
dataset.log()

# Construct a manifest with the paired FASTQ files in the input

# If the files were added with a samplesheet.csv, they will include
# an index indicating which line of the samplesheet.csv they were in
if "sampleIndex" in dataset.files.columns.values:

    # Reconstruct the manifest
    manifest = dataset.files.pivot(
        index=["sampleIndex", "sample"],
        columns="read",
        values="file"
    ).rename(
        columns=lambda i: f"fastq_{int(i)}"
    ).sort_index(
    ).sort_index(
        axis=1
    ).reset_index(
    ).drop(
        columns=["sampleIndex"]
    )

# If the files weren't added with a samplesheet.csv, then we will
# use different logic to construct the sample sheet.
# This is intended to capture the scenario when there are multiple
# pairs of FASTQs for any samples.
else:

    manifest = []

    # Iterate over each file
    for sample_name, sample_files in dataset.files.groupby("sample"):

        # Make a sorted list of the files available for this sample
        file_list = sample_files["file"].sort_values().tolist()

        # Add pairs of files to the manifest
        while len(file_list) > 0:

            assert len(file_list) >= 2, f"Unexpected odd number of files found for sample {sample_name}"

            # Add the first two files to the manifest
            manifest.append(
                dict(
                    sample=sample_name,
                    fastq_1=file_list[0],
                    fastq_2=file_list[1]
                )
            )

            # Remove those files from the list
            file_list = file_list[2:]

    manifest = pd.DataFrame(manifest)

# Write to the dataset's config/ folder, the location process-input.json maps
# paired_manifest to. This script predates the Cirro SDK and has no ds.params,
# so the path is built from the dataset root it already reads.
manifest_path = f"{dataset.s3_dataset}/config/manifest.csv"

manifest.reindex(
    columns=[
        "sample",
        "fastq_1",
        "fastq_2"
    ]
).to_csv(
    manifest_path,
    index=None
)
print(f"Wrote out {manifest.shape[0]:,} lines to {manifest_path}")

# This script predates the Cirro SDK and changes no parameters, but the HealthOmics
# pre-process Lambda fails the run when params.json is absent and only the SDK writes
# that file, so write it unchanged.
sdk_ds = PreprocessDataset.from_running()
sdk_ds.keep_params(list(sdk_ds.params.keys()))
