#!/usr/bin/env python3

import pandas as pd
import os

# Get the environment variable with the S3 path to the dataset samplesheet
pw_s3_dataset = os.getenv("PW_S3_DATASET")

# Read the samplesheet for the dataset
samplesheet = pd.read_csv(f"{pw_s3_dataset}/config/samplesheet.csv")

print(f"Read in samplesheet with {samplesheet.shape[0]:,} rows and {samplesheet.shape[1]:,} columns")
assert samplesheet.shape[0] > 0, "No files detected -- there may be an error with data ingest"

print(samplesheet.head())

# Write out the samplesheet to a local file
samplesheet.to_csv("samplesheet.csv", index=None)
