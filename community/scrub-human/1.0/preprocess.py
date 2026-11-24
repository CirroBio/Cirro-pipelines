#!/usr/bin/env python3

import boto3
from cirro.helpers.preprocess_dataset import PreprocessDataset
from cirro.models.s3_path import S3Path


def make_file_list(ds: PreprocessDataset):

    ds.logger.info("Files annotated in the dataset:")
    ds.logger.info(ds.files.to_csv(index=None))

    # Filter out any index files that may have been uploaded
    ds.files = ds.files.loc[
        ds.files.apply(
            lambda r: r.get('readType', 'R') == 'R',
            axis=1
        )
    ]

    # Make a simple table with each of the files listed on single row
    # Write to the dataset's config/ folder (mapped in process-input.json)
    ds.files.reindex(
        columns=["file"]
    ).to_csv(
        ds.params["input_filelist"],
        header=False,
        index=False
    )


def copy_file_s3(s3, source_file, dest_folder):
    if not dest_folder.endswith("/"):
        dest_folder = dest_folder + "/"
    dest_file = dest_folder + source_file.split("/")[-1]

    source = S3Path(source_file)
    dest = S3Path(dest_file)

    s3.meta.client.copy(
        {
            "Bucket": source.bucket,
            "Key": source.key,
        },
        dest.bucket,
        dest.key
    )


def copy_samplesheet(ds: PreprocessDataset):
    """Copy the samplesheet from the input dataset, if it exists."""

    s3 = boto3.resource('s3')
    try:
        copy_file_s3(
            s3,
            ds.params["input_folder"] + "/samplesheet.csv",
            ds.params["output_folder"]
        )
    except Exception as e:
        ds.logger.info("Did not find samplesheet.csv from the input dataset")
        ds.logger.info(str(e))
        return

    ds.logger.info("Successfully copied samplesheet.csv to output directory")


if __name__ == "__main__":
    # Set up the Cirro Dataset
    ds = PreprocessDataset.from_running()

    # Set up the list of files to process
    make_file_list(ds)

    # Copy over the input samplesheet.csv, if it exists
    copy_samplesheet(ds)

    # log
    ds.logger.info(ds.params)

    # Force params.json to be written: the HealthOmics pre-process Lambda fails the run
    # when the file is absent, and the SDK writes it only when a parameter changes.
    ds.keep_params(list(ds.params.keys()))
