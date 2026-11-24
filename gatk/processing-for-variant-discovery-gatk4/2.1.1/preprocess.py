#!/usr/bin/env python3

from __future__ import annotations

import json
from typing import TYPE_CHECKING

import pandas as pd

# The SDK is imported in __main__, so test_preprocess.py runs without it
if TYPE_CHECKING:
    from cirro.helpers.preprocess_dataset import PreprocessDataset

WORKFLOW = "PreProcessingForVariantDiscovery_GATK4"

KNOWN_INDELS = [
    "GATK/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
    "GATK/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz"
]


def read_samplesheet(path) -> pd.DataFrame:
    """Read the table of unmapped BAMs filled out in the form, keeping every value as text."""

    # The form field is marked required, but an empty table can still be submitted,
    # in which case Cirro either omits the param or writes an empty file
    if not path:
        raise ValueError("No samplesheet was provided: add at least one row with an unmapped BAM")
    try:
        samplesheet = pd.read_csv(path, dtype=str, keep_default_na=False)
    except pd.errors.EmptyDataError:
        samplesheet = pd.DataFrame()
    if samplesheet.empty:
        raise ValueError("The samplesheet is empty: add at least one row with an unmapped BAM")
    return samplesheet


def group_samples(samplesheet: pd.DataFrame) -> list:
    """The unmapped BAMs of each sample, as (sample_name, [BAM, ...]) in table order."""

    return list(
        samplesheet.groupby("sample_name", sort=False)["flowcell_unmapped_bams"].apply(list).items()
    )


def known_indels(references_base: str) -> dict:
    """The known indel sites, which the workflow takes as lists of files."""

    return {
        f"{WORKFLOW}.known_indels_sites_VCFs": [f"{references_base}/{fp}" for fp in KNOWN_INDELS],
        f"{WORKFLOW}.known_indels_sites_indices": [f"{references_base}/{fp}.tbi" for fp in KNOWN_INDELS]
    }


def build_inputs(samples: list, list_uris: list, params: dict) -> list:
    """One set of workflow inputs per sample, each reading its BAMs from a list file."""

    workflow_params = {
        kw: val
        for kw, val in params.items()
        if kw.startswith(WORKFLOW)
    }

    return [
        {
            f"{WORKFLOW}.sample_name": sample,
            f"{WORKFLOW}.flowcell_unmapped_bams_list": list_uri,
            **workflow_params
        }
        for (sample, _), list_uri in zip(samples, list_uris)
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

    samples = group_samples(read_samplesheet(ds.params.get("samplesheet")))

    # Stage each sample's list of BAMs beside the samplesheet, in the dataset's config/ folder
    config_dir = ds.params["samplesheet"].rsplit("/", 1)[0]
    list_uris = [f"{config_dir}/flowcell_unmapped_bams.{i}.list" for i in range(len(samples))]
    for (_, bams), list_uri in zip(samples, list_uris):
        pd.Series(bams).to_csv(list_uri, index=False, header=False)

    params = {**ds.params, **known_indels(ds.references_base)}
    all_inputs = build_inputs(samples, list_uris, params)

    # The samplesheet is only read here, so keep it out of the Cromwell options
    ds.remove_param("samplesheet")

    # Write out the complete set of inputs
    write_json("inputs.json", all_inputs)

    # Write out each individual sample
    for i, single_input in enumerate(all_inputs):
        write_json(f"inputs.{i}.json", single_input)

    setup_options(ds)
