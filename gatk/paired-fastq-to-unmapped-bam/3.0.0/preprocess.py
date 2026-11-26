#!/usr/bin/env python3

from __future__ import annotations

import json
from datetime import datetime
from typing import TYPE_CHECKING

import pandas as pd

# The SDK is imported in __main__, so test_preprocess.py runs without it
if TYPE_CHECKING:
    from cirro.helpers.preprocess_dataset import PreprocessDataset

WORKFLOW = "ConvertPairedFastQsToUnmappedBamWf"

# Read group fields that take the sample name when left blank
DEFAULT_TO_SAMPLE = ["readgroup_name", "library_name", "platform_unit"]


def read_samplesheet(path) -> pd.DataFrame:
    """Read the table of read groups filled out in the form, keeping every value as text."""

    # The form field is marked required, but an empty table can still be submitted,
    # in which case Cirro either omits the param or writes an empty file
    if not path:
        raise ValueError("No samplesheet was provided: add at least one row of FASTQ files")
    try:
        samplesheet = pd.read_csv(path, dtype=str, keep_default_na=False)
    except pd.errors.EmptyDataError:
        samplesheet = pd.DataFrame()
    if samplesheet.empty:
        raise ValueError("The samplesheet is empty: add at least one row of FASTQ files")
    return samplesheet


def fill_read_groups(samplesheet: pd.DataFrame) -> pd.DataFrame:
    """Fill blank read group fields with the sample name."""

    samplesheet = samplesheet.assign(**{
        col: samplesheet[col].where(samplesheet[col] != "", samplesheet["sample"])
        for col in DEFAULT_TO_SAMPLE
    })

    # Each unmapped BAM is named after its read group, so two rows with the same
    # read group would write to the same output file
    duplicated = samplesheet["readgroup_name"][samplesheet["readgroup_name"].duplicated()].unique()
    if len(duplicated) > 0:
        raise ValueError(
            f"Read group names must be unique, but these are used more than once: {', '.join(duplicated)}. "
            "Give each pair of FASTQ files of a sample its own read group"
        )
    return samplesheet


def shared_inputs(params: dict, today: str) -> dict:
    """The workflow inputs which are the same for every read group."""

    inputs = {
        kw: val
        for kw, val in params.items()
        if kw.startswith(WORKFLOW)
    }

    # The workflow passes the sequencing center to FastqToSam unquoted
    inputs[f"{WORKFLOW}.sequencing_center"] = (
        inputs.get(f"{WORKFLOW}.sequencing_center") or "Sequencing_Center"
    ).replace(" ", "_")

    # If run_date was not provided, use today's date
    if not inputs.get(f"{WORKFLOW}.run_date"):
        inputs[f"{WORKFLOW}.run_date"] = today

    return inputs


def build_inputs(samplesheet: pd.DataFrame, shared: dict) -> list:
    """One set of workflow inputs per read group."""

    return [
        {
            f"{WORKFLOW}.sample_name": row["sample"],
            **{
                f"{WORKFLOW}.{kw}": row[kw]
                for kw in ["readgroup_name", "library_name", "platform_unit", "fastq_1", "fastq_2"]
            },
            **shared
        }
        for _, row in samplesheet.iterrows()
    ]


def write_json(fp, obj, indent=4):

    with open(fp, "wt") as handle:
        json.dump(obj, handle, indent=indent)


def setup_options(ds: PreprocessDataset):

    from cirro.models.s3_path import S3Path

    # Add default params
    ds.add_param("memory_retry_multiplier", 2.0)
    ds.add_param("use_relative_output_paths", True)

    # Set up the scriptBucketName
    ds.add_param(
        "scriptBucketName",
        S3Path(ds.params['final_workflow_outputs_dir']).bucket
    )

    # Isolate the options arguments for the workflow
    options = {
        kw: val
        for kw, val in ds.params.items()
        if not kw.startswith(WORKFLOW)
    }

    # Write out
    write_json("options.json", options)


if __name__ == "__main__":

    from cirro.helpers.preprocess_dataset import PreprocessDataset

    ds = PreprocessDataset.from_running()

    samplesheet = fill_read_groups(read_samplesheet(ds.params.get("samplesheet")))
    all_inputs = build_inputs(samplesheet, shared_inputs(ds.params, datetime.now().isoformat()))

    # The samplesheet is only read here, so keep it out of the Cromwell options
    ds.remove_param("samplesheet")

    # Write out the complete set of inputs
    write_json("inputs.json", all_inputs)

    # Write out each individual file pair
    for i, single_input in enumerate(all_inputs):
        write_json(f"inputs.{i}.json", single_input)

    setup_options(ds)
