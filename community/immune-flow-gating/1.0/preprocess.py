#!/usr/bin/env python3

import json
import boto3
from cirro.helpers.preprocess_dataset import PreprocessDataset
from cirro.models.s3_path import S3Path


def write_text(uri: str, content: str):
    """Write a text file to the S3 location given by uri."""
    path = S3Path(uri)
    assert path.valid, f"Not an S3 URI: {uri}"
    boto3.client("s3").put_object(
        Bucket=path.bucket,
        Key=path.key,
        Body=content.encode()
    )


def build_channel_map(ds: PreprocessDataset):

    # Build a mapping of the expected channels for each marker
    marker_map = dict(
        CCR7="PerCP-Cy5-5-A",
        CD3="Alexa Fluor 700-A",
        CD4="APC-eFluor 780-A",
        CD8="BV650-A",
        CD14="APC-A",
        CD19="PE-Cy7-A",
        CD25="FITC-A",
        CD45RA="eFluor 450-A",
        CD56="PE-A",
        LD="V500-A"
    )

    # Build the channel_map.json file required
    channel_map = []

    for marker, cname in marker_map.items():

        # If the user selected a non-default channel name
        if ds.params[marker] != cname:

            # Add the mapping
            channel_map.append(
                dict(
                    channels=ds.params[marker],
                    alias=cname
                )
            )
            ds.logger.info(
                f"Mapping {cname} -> {ds.params[marker]} (for {marker})"
            )

    # Write to the dataset's config/ folder (mapped in process-input.json)
    write_text(ds.params["channel_map"], json.dumps(channel_map, indent=4))


def build_samplesheet(ds: PreprocessDataset):

    # Make a table with the columns:
    # sample,file_1
    samplesheet = ds.files.assign(
        file_1=ds.files["file"].apply(
            lambda s: s.rsplit("/", 1)[-1]
        )
    ).drop(
        columns=["file"]
    ).set_index(
        "sample"
    )

    if "sampleIndex" in samplesheet.columns.values:
        samplesheet = samplesheet.drop(columns=["sampleIndex"])

    # Add in any metadata provided for the samples
    for cname, cvals in ds.samplesheet.set_index("sample").items():
        samplesheet = samplesheet.assign(**{
            cname: cvals
        })

    # Log the samplesheet
    ds.logger.info(samplesheet.to_csv())

    # Write to the dataset's config/ folder (mapped in process-input.json)
    samplesheet.to_csv(ds.params["samplesheet"])


if __name__ == "__main__":
    ds = PreprocessDataset.from_running()

    build_channel_map(ds)
    build_samplesheet(ds)

    # Consumed above; not parameters the workflow declares.
    for param in ("CCR7", "CD3", "CD4", "CD8", "CD14", "CD19", "CD25", "CD45RA", "CD56", "LD"):
        ds.remove_param(param, force=True)

    # Force params.json to be written: the HealthOmics pre-process Lambda fails the run
    # when the file is absent, and the SDK writes it only when a parameter changes.
    ds.keep_params(list(ds.params.keys()))
