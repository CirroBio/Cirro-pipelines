#!/usr/bin/env python3

import json
from cirro.helpers.preprocess_dataset import PreprocessDataset
import pandas as pd
from datetime import datetime
from cirro.models.s3_path import S3Path


def setup_inputs(ds: PreprocessDataset):

    # Make sure that the sequencing center field does not contain spaces
    ds.params[
        "ConvertPairedFastQsToUnmappedBamWf.sequencing_center"
    ] = ds.params.get(
        "ConvertPairedFastQsToUnmappedBamWf.sequencing_center",
        "Sequencing Center"
    ).replace(" ", "_")

    # If run_date was not provided, use today's date
    apply_run_date(ds)

    # Make a combined set of inputs with all of the R1/R2 file pairs
    all_inputs = [
        {
            **single_input,
            **{
                kw: val
                for kw, val in ds.params.items()
                if kw.startswith("ConvertPairedFastQsToUnmappedBamWf")
            }
        }
        for single_input in yield_single_inputs(ds)
    ]

    # Write out the complete set of inputs
    write_json("inputs.json", all_inputs)

    # Write out each individual file pair
    for i, input in enumerate(all_inputs):
        write_json(f"inputs.{i}.json", input)


def write_json(fp, obj, indent=4):

    with open(fp, "wt") as handle:
        json.dump(obj, handle, indent=indent)


def yield_single_inputs(ds: PreprocessDataset):

    # Make a wide DataFrame with columns:
    # sample, lane, fastq_1, fastq_2
    wide_df = make_wide_df(ds)

    for _, row in wide_df.iterrows():

        dat = dict(
            readgroup_name=row['sample'],
            sample_name=row['sample'],
            library_name=row['sample'],
            platform_unit=row['sample'],
            fastq_1=row['fastq_1'],
            fastq_2=row['fastq_2'],
        )

        yield {
            f"ConvertPairedFastQsToUnmappedBamWf.{kw}": val
            for kw, val in dat.items()
        }


def apply_run_date(ds: PreprocessDataset):

    # Get the run_date value provided by the user
    run_date = ds.params.get("ConvertPairedFastQsToUnmappedBamWf.run_date")

    # If the user did not provide a value
    if not run_date:

        # Use today's date
        ds.add_param(
            "ConvertPairedFastQsToUnmappedBamWf.run_date",
            datetime.now().isoformat(),
            overwrite=True
        )


def setup_options(ds: PreprocessDataset):

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
        if not kw.startswith("ConvertPairedFastQsToUnmappedBamWf")
    }

    # Write out
    with open("options.json", "wt") as handle:
        json.dump(options, handle, indent=4)


def make_wide_df(ds: PreprocessDataset) -> pd.DataFrame:
    """Format the FASTQ inputs in wide format."""

    # Format as a wide dataset
    ds.logger.info("Creating paired table of inputs")
    wide_df = ds.files.assign(
        readType=ds.files.reindex(columns=["readType"])["readType"].fillna("R")
    ).query(
        "readType == 'R'"
    ).reindex(
        columns=["sampleIndex", "sample", "lane", "read", "file"]
    ).pivot(
        index=["sampleIndex", "sample", "lane"],
        columns="read",
        values="file"
    ).rename(
        columns=lambda i: f"fastq_{int(i)}"
    ).reset_index()
    ds.logger.info("Creating paired table of inputs - DONE")
    ds.logger.info(wide_df.to_csv(index=None))

    return wide_df


if __name__ == "__main__":
    ds = PreprocessDataset.from_running()
    setup_inputs(ds)
    setup_options(ds)
