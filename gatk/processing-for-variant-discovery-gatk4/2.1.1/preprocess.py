#!/usr/bin/env python3

import json
from cirro.helpers.preprocess_dataset import PreprocessDataset
from cirro.models.s3_path import S3Path


def setup_inputs(ds: PreprocessDataset):

    # There are two params which need to be set up as lists instead of strings
    for kw_suffix in ["known_indels_sites_VCFs", "known_indels_sites_indices"]:
        kw = f"PreProcessingForVariantDiscovery_GATK4.{kw_suffix}"
        ds.add_param(kw, ds.params[kw].split(","), overwrite=True)

    # Make a combined set of inputs with each of the BAM files
    all_inputs = [
        {
            **single_input,
            **{
                kw: val
                for kw, val in ds.params.items()
                if kw.startswith("PreProcessingForVariantDiscovery_GATK4")
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

    ds.logger.info("Input files:")
    ds.logger.info(ds.files.to_csv(index=None))
    assert ds.files.shape[0] > 0, "ERROR: no input files found"

    # Each of the files in the input is a BAM file for a single sample
    for sample, bam_fp in ds.files.set_index("sample")["file"].items():

        # Make sure that the file is BAM
        if not bam_fp.endswith(".bam"):
            ds.logger.info(f"Skipping non-BAM file: {bam_fp}")
            continue

        ds.logger.info(f"Sample: {sample}")
        ds.logger.info(f"BAM: {bam_fp}")

        # Make a file with the BAM file path
        with open(f"{sample}.bam_list.txt", "w") as handle:
            handle.write(bam_fp)

        dat = dict(
            sample_name=sample,
            flowcell_unmapped_bams_list=f"{sample}.bam_list.txt"
        )

        yield {
            f"PreProcessingForVariantDiscovery_GATK4.{kw}": val
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
        if not kw.startswith("PreProcessingForVariantDiscovery_GATK4")
    }

    # Write out
    with open("options.json", "wt") as handle:
        json.dump(options, handle, indent=4)


if __name__ == "__main__":
    ds = PreprocessDataset.from_running()
    setup_inputs(ds)
    setup_options(ds)
