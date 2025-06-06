#!/usr/bin/env python3

from cirro.api.models.s3_path import S3Path
from cirro.helpers.preprocess_dataset import PreprocessDataset
from logging import Logger
from pathlib import Path
from typing import Union
import boto3
import json
import pandas as pd


def read_json(path):
    """Read a JSON from S3."""

    # Make the full S3 path
    s3_path = S3Path(path)

    if s3_path.valid:
        s3 = boto3.client('s3')
        retr = s3.get_object(Bucket=s3_path.bucket, Key=s3_path.key)
        text = retr['Body'].read().decode()
    else:
        with Path(path).open() as handle:
            text = handle.read()

    # Parse JSON
    return json.loads(text)


def read_input_dataset(ds: PreprocessDataset):

    # Get the list of files in the parent dataset
    manifest_fp = ds.params["input_dataset"] + "/web/manifest.json"
    ds.logger.info(f"Reading files from input_dataset ({manifest_fp})")
    manifest = read_json(manifest_fp)

    # Add the full file path
    for i in manifest['files']:
        i['file'] = ds.params["input_dataset"] + "/" + i['file']

    return manifest


def parse_fastq(r: pd.Series, suffix=".fastq.gz") -> Union[str, None]:
    """
    Check that the filename conforms to the expected structure.
    e.g. SRX3815835_SRR6860803_1.fastq.gz
    """
    fn: str = r['file'].rsplit("/", 1)[-1]

    if not fn.endswith(suffix):
        return
    elif len(fn[:-len(suffix)].split("_")) != 3:
        return
    else:
        return dict(
            zip(
                ['sample', 'run', 'spot'],
                fn[:-len(suffix)].split("_")
            ),
            **r
        )


def is_fastq(r: pd.Series) -> bool:
    return parse_fastq(r) is not None


def main():

    # Instantiate the Cirro dataset object
    ds = PreprocessDataset.from_running()

    # Read the list of files from the input dataset
    manifest = read_input_dataset(ds)

    # Parse the files, filtering to FASTQ files and parsing
    # the sample, run, and read index
    files = pd.DataFrame([
        parse_fastq(r)
        for r in manifest["files"]
        if is_fastq(r)
    ])

    # Assign the R1 and R2 based on the median filesize,
    # and use that information to format a samplesheet
    samplesheet = format_samplesheet(files, ds.logger)

    # Write to disk
    samplesheet.to_csv("samplesheet.csv", index=None)

    # Point the workflow to the spreadsheet
    ds.add_param("samplesheet", "samplesheet.csv")

    # Log the parameters present
    for k, v in ds.params.items():
        ds.logger.info(f"{k}: {v}")


def format_samplesheet(files: pd.DataFrame, logger: Logger) -> pd.DataFrame:
    """
    The approach for identifying the appropriate spots was derived from
    https://kb.10xgenomics.com/hc/en-us/articles/115003802691-How-do-I-prepare-Sequence-Read-Archive-SRA-data-from-NCBI-for-Cell-Ranger-
    """

    # Get the median filesize for each of the spot
    logger.info("Summarizing file size per spot")
    logger.info(files.head().to_csv(index=None))
    median_size = (
        files
        .groupby('spot')
        .apply(lambda d: d['size'].median())
        .sort_values(ascending=False)
    )

    msg = "Only 1 spot found -- cellranger requires paired-end reads"
    assert median_size.shape[0] >= 2, msg

    logger.info("Median file size - per spot")
    logger.info(median_size)

    logger.info("Pivoting to wide samplesheet")

    # R2 is the longest read, and R1 is the second-longest
    samplesheet = (
        files
        .pivot(
            index=["sample", "run"],
            columns="spot",
            values="file"
        )
        .rename(
            columns={
                median_size.index.values[0]: "fastq_2",
                median_size.index.values[1]: "fastq_1"
            }
        )
        .reset_index()
        .reindex(columns=["sample", "fastq_1", "fastq_2"])
        .dropna()
    )

    assert samplesheet.shape[0] > 0, "No paired-end reads found"

    return samplesheet


if __name__ == "__main__":
    main()
