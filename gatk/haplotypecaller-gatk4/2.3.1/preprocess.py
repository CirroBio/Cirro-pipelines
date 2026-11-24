#!/usr/bin/env python3

from __future__ import annotations

import json
from typing import TYPE_CHECKING

import pandas as pd

# The SDK is imported in __main__, so test_preprocess.py runs without it
if TYPE_CHECKING:
    from cirro.helpers.preprocess_dataset import PreprocessDataset

WORKFLOW = "HaplotypeCallerGvcf_GATK4"


def read_samplesheet(path) -> pd.DataFrame:
    """Read the table of alignments filled out in the form, keeping every value as text."""

    # The form field is marked required, but an empty table can still be submitted,
    # in which case Cirro either omits the param or writes an empty file
    if not path:
        raise ValueError("No samplesheet was provided: add at least one row with a BAM and its index")
    try:
        samplesheet = pd.read_csv(path, dtype=str, keep_default_na=False)
    except pd.errors.EmptyDataError:
        samplesheet = pd.DataFrame()
    if samplesheet.empty:
        raise ValueError("The samplesheet is empty: add at least one row with a BAM and its index")
    return samplesheet


def build_inputs(samplesheet: pd.DataFrame, params: dict) -> list:
    """One set of workflow inputs per row of the samplesheet."""

    workflow_params = {
        kw: val
        for kw, val in params.items()
        if kw.startswith(WORKFLOW)
    }

    all_inputs = []
    for _, row in samplesheet.iterrows():

        dat = dict(
            input_bam=row["input_bam"],
            input_bam_index=row["input_bam_index"]
        )

        # Set inputs at the task level
        for task in ["HaplotypeCaller", "CramToBamTask", "MergeGVCFs"]:

            # Hardcode the reference size
            dat[f"{task}.disk_space_gb"] = 100

            # Set the preemptible attempts
            dat[f"{task}.preemptible_attempts"] = 3

        all_inputs.append({
            **{f"{WORKFLOW}.{kw}": val for kw, val in dat.items()},
            **workflow_params
        })

    return all_inputs


def write_json(fp, obj, indent=4) -> None:

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

    all_inputs = build_inputs(read_samplesheet(ds.params.get("samplesheet")), ds.params)

    # The samplesheet is only read here, so keep it out of the Cromwell options
    ds.remove_param("samplesheet")

    # Write out the complete set of inputs
    write_json("inputs.json", all_inputs)

    # Write out each individual sample
    for i, single_input in enumerate(all_inputs):
        write_json(f"inputs.{i}.json", single_input)

    setup_options(ds)
