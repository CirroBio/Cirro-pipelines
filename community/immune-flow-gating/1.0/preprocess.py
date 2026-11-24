#!/usr/bin/env python3

import json

import pandas as pd

# The channel the workflow expects for each marker
MARKER_CHANNELS = dict(
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


def build_channel_map(params: dict) -> list:
    """Alias each marker's channel, as chosen in the form, to the channel the workflow expects."""
    return [
        dict(channels=params[marker], alias=channel)
        for marker, channel in MARKER_CHANNELS.items()
        if params[marker] != channel
    ]


def read_samplesheet(path) -> pd.DataFrame:
    """Read the samplesheet filled out in the form, keeping every value as text."""

    # The form field is marked required, but an empty table can still be submitted,
    # in which case Cirro either omits the param or writes an empty file
    if not path:
        raise ValueError("No samplesheet was provided: add one row per FCS file")
    try:
        samplesheet = pd.read_csv(path, dtype=str, keep_default_na=False)
    except pd.errors.EmptyDataError:
        samplesheet = pd.DataFrame()
    if samplesheet.empty:
        raise ValueError("The samplesheet is empty: add one row per FCS file")
    return samplesheet


def build_workflow_samplesheet(samplesheet: pd.DataFrame) -> pd.DataFrame:
    """The samplesheet in the shape the workflow reads: each FCS file by its name alone.

    The workflow stages every FCS file into one folder and finds each row's file there,
    so file names must be unique. A batch column is only kept when filled in, since the
    workflow gates by batch whenever the column is present.
    """
    file_names = samplesheet["fcs"].str.rsplit("/", n=1).str[-1]
    duplicated = file_names[file_names.duplicated()].unique()
    if len(duplicated) > 0:
        raise ValueError(
            f"Each FCS file name can be listed only once. Listed more than once: {', '.join(duplicated)}"
        )

    workflow_samplesheet = pd.DataFrame(dict(sample=samplesheet["sample"], file_1=file_names))
    if (samplesheet["batch"] != "").any():
        workflow_samplesheet["batch"] = samplesheet["batch"]
    return workflow_samplesheet


if __name__ == "__main__":

    # Imported here so that test_preprocess.py runs without the SDK
    import boto3
    from cirro.helpers.preprocess_dataset import PreprocessDataset
    from cirro.models.s3_path import S3Path

    ds = PreprocessDataset.from_running()

    channel_map = build_channel_map(ds.params)
    ds.logger.info(f"Channel map: {channel_map}")
    # Write to the dataset's config/ folder (mapped in process-input.json)
    channel_map_path = S3Path(ds.params["channel_map"])
    boto3.client("s3").put_object(
        Bucket=channel_map_path.bucket,
        Key=channel_map_path.key,
        Body=json.dumps(channel_map, indent=4).encode()
    )

    samplesheet = read_samplesheet(ds.params.get("fcs_samplesheet"))
    # Write to the dataset's config/ folder (mapped in process-input.json)
    build_workflow_samplesheet(samplesheet).to_csv(ds.params["samplesheet"], index=None)
    ds.add_param("fcs", ",".join(samplesheet["fcs"]))

    # Consumed above; not parameters the workflow declares.
    for param in ("fcs_samplesheet", *MARKER_CHANNELS):
        ds.remove_param(param, force=True)
