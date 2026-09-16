#!/usr/bin/env python3

import json
from cirro.helpers.preprocess_dataset import PreprocessDataset
from cirro.models.s3_path import S3Path


def setup_inputs(ds: PreprocessDataset):

    # Make a combined set of inputs with each of the BAM files
    all_inputs = [
        {
            **single_input,
            **{
                kw: val
                for kw, val in ds.params.items()
                if kw.startswith("HaplotypeCallerGvcf_GATK4")
            }
        }
        for single_input in yield_single_inputs(ds)
    ]

    # Raise an error if no inputs are found
    assert len(all_inputs) > 0, "No inputs found -- stopping execution"

    # Write out the complete set of inputs
    write_json("inputs.json", all_inputs)

    # Write out each individual file pair
    for i, input in enumerate(all_inputs):
        write_json(f"inputs.{i}.json", input)


def write_json(fp, obj, indent=4) -> None:

    with open(fp, "wt") as handle:
        json.dump(obj, handle, indent=indent)


def mark_filetype(fp: str) -> str:
    if fp.endswith(".bam"):
        return "bam"
    elif fp.endswith(".bai"):
        return "bai"
    else:
        return "other"


def yield_single_inputs(ds: PreprocessDataset) -> dict:
    """The input files will have both a BAM and a BAI for each sample."""

    # Get the BAM and BAI for each sample
    for sample, files in ds.files.groupby("sample"):

        file_dict = {
            mark_filetype(fp): fp
            for fp in files['file'].values
        }

        bam_fp = file_dict.get("bam")
        bai_fp = file_dict.get("bai")

        ds.logger.info(f"Sample: {sample}")
        if bam_fp is None:
            ds.logger.info("No BAM file found, skipping")
            continue
        ds.logger.info(f"BAM: {bam_fp}")
        if bai_fp is None:
            ds.logger.info("No BAM Index file found, skipping")
            continue

        # Set the params, which will all have the
        # workflow-level prefix appended (below)
        dat = dict(
            input_bam=bam_fp,
            input_bam_index=bai_fp
        )

        # Set inputs at the task level
        for task in ["HaplotypeCaller", "CramToBamTask", "MergeGVCFs"]:

            # Hardcode the reference size
            dat[f"{task}.disk_space_gb"] = 100

            # Set the preemptible attempts
            dat[f"{task}.preemptible_attempts"] = 3

        yield {
            f"HaplotypeCallerGvcf_GATK4.{kw}": val
            for kw, val in dat.items()
        }


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
        if not kw.startswith("HaplotypeCallerGvcf_GATK4")
    }

    # Write out
    with open("options.json", "wt") as handle:
        json.dump(options, handle, indent=4)


if __name__ == "__main__":
    ds = PreprocessDataset.from_running()
    setup_inputs(ds)
    setup_options(ds)
