#!/usr/bin/env python3

from __future__ import annotations

import json
from typing import TYPE_CHECKING

import pandas as pd

# The SDK is imported in __main__, so test_preprocess.py runs without it
if TYPE_CHECKING:
    from cirro.helpers.preprocess_dataset import PreprocessDataset

WORKFLOW = "WholeGenomeGermlineSingleSample"


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


def build_inputs(samplesheet: pd.DataFrame, params: dict) -> list:
    """One set of workflow inputs per sample, with all of that sample's unmapped BAMs."""

    workflow_params = {
        kw: val
        for kw, val in params.items()
        if kw.startswith(WORKFLOW)
    }

    return [
        {
            f"{WORKFLOW}.sample_and_unmapped_bams": dict(
                sample_name=sample,
                base_file_name=sample,
                final_gvcf_base_name=sample,
                unmapped_bam_suffix=".bam",
                flowcell_unmapped_bams=bams.tolist()
            ),
            **workflow_params
        }
        for sample, bams in samplesheet.groupby("sample", sort=False)["flowcell_unmapped_bams"]
    ]


def write_json(fp, obj, indent=4):

    with open(fp, "wt") as handle:
        json.dump(obj, handle, indent=indent)


def setup_options(ds: PreprocessDataset):

    from cirro.models.s3_path import S3Path

    # Add default params
    ds.add_param("memory_retry_multiplier", 2.0)
    ds.add_param("use_relative_output_paths", True)
    ds.add_param("read_from_cache", True)
    ds.add_param("write_to_cache", True)

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

    # If running in DRAGEN mode - Functional Equivalence
    if ds.params.get("dragen_mode") == "Enabled (Functional Equivalence)":

        ds.logger.info("Running in DRAGEN Mode (Functional Equivalence)")
        ds.add_param("dragen_functional_equivalence_mode", True)

    # If running in DRAGEN mode - Maximum Quality
    elif ds.params.get("dragen_mode") == "Enabled (Maximum Quality)":

        ds.logger.info("Running in DRAGEN Mode (Maximum Quality)")
        ds.add_param("dragen_maximum_quality_mode", True)

    # If running in normal mode
    else:
        msg = f"Did not recognize dragen_mode={ds.params['dragen_mode']}"
        assert ds.params["dragen_mode"] == "Disabled", msg

    ds.remove_param("dragen_mode", force=True)

    # Write out
    with open("options.json", "wt") as handle:
        json.dump(options, handle, indent=4)


def setup_reference(ds: PreprocessDataset) -> None:
    """Set up the reference genome information."""

    avail_references = {
        "Homo sapiens (GRCh38)": {
            "contamination_sites_ud": f"{ds.references_base}/GATK/hg38/v0/contamination-resources/1000g/1000g.phase3.100k.b38.vcf.gz.dat.UD",
            "contamination_sites_bed": f"{ds.references_base}/GATK/hg38/v0/contamination-resources/1000g/1000g.phase3.100k.b38.vcf.gz.dat.bed",
            "contamination_sites_mu": f"{ds.references_base}/GATK/hg38/v0/contamination-resources/1000g/1000g.phase3.100k.b38.vcf.gz.dat.mu",
            "calling_interval_list": f"{ds.references_base}/GATK/hg38/v0/wgs_calling_regions.hg38.interval_list",
            "reference_fasta" : {
                "ref_dict": f"{ds.references_base}/GATK/hg38/v0/dragen_reference/Homo_sapiens_assembly38_masked.dict",
                "ref_fasta": f"{ds.references_base}/GATK/hg38/v0/dragen_reference/Homo_sapiens_assembly38_masked.fasta",
                "ref_fasta_index": f"{ds.references_base}/GATK/hg38/v0/dragen_reference/Homo_sapiens_assembly38_masked.fasta.fai",
                "ref_sa": f"{ds.references_base}/GATK/hg38/v0/dragen_reference/Homo_sapiens_assembly38_masked.fasta.64.sa",
                "ref_alt": f"{ds.references_base}/GATK/hg38/v0/dragen_reference/Homo_sapiens_assembly38_masked.fasta.alt",
                "ref_amb": f"{ds.references_base}/GATK/hg38/v0/dragen_reference/Homo_sapiens_assembly38_masked.fasta.64.amb",
                "ref_bwt": f"{ds.references_base}/GATK/hg38/v0/dragen_reference/Homo_sapiens_assembly38_masked.fasta.64.bwt",
                "ref_ann": f"{ds.references_base}/GATK/hg38/v0/dragen_reference/Homo_sapiens_assembly38_masked.fasta.64.ann",
                "ref_pac": f"{ds.references_base}/GATK/hg38/v0/dragen_reference/Homo_sapiens_assembly38_masked.fasta.64.pac"
            },
            "known_indels_sites_vcfs": [
                f"{ds.references_base}/GATK/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
                f"{ds.references_base}/GATK/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz"
            ],
            "known_indels_sites_indices": [
                f"{ds.references_base}/GATK/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi",
                f"{ds.references_base}/GATK/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi"
            ],
            "dbsnp_vcf": f"{ds.references_base}/GATK/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf",
            "dbsnp_vcf_index": f"{ds.references_base}/GATK/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx",
            "evaluation_interval_list": f"{ds.references_base}/GATK/hg38/v0/wgs_evaluation_regions.hg38.interval_list",
            "haplotype_database_file": f"{ds.references_base}/GATK/hg38/v0/Homo_sapiens_assembly38.haplotype_database.txt"
        }
    }

    msg = f"No references defined for genome {ds.params['genome']}"
    assert ds.params["genome"] in avail_references, msg

    ds.add_param(
        "WholeGenomeGermlineSingleSample.references",
        avail_references[ds.params["genome"]]
    )

    wgs_coverage_interval_list = {
        "Homo sapiens (GRCh38)": f"{ds.references_base}/GATK/hg38/v0/wgs_coverage_regions.hg38.interval_list"
    }
    msg = f"No wgs_coverage_interval_list defined for genome {ds.params['genome']}"
    assert ds.params["genome"] in wgs_coverage_interval_list, msg

    ds.add_param(
        "WholeGenomeGermlineSingleSample.wgs_coverage_interval_list",
        wgs_coverage_interval_list[ds.params['genome']]
    )

    ds.add_param(
        "WholeGenomeGermlineSingleSample.allow_empty_ref_alt",
        True
    )

    ds.add_param(
        "WholeGenomeGermlineSingleSample.scatter_settings",
        {
            "haplotype_scatter_count": 50,
            "break_bands_at_multiples_of": 1000000
        }
    )

    ds.add_param(
        "WholeGenomeGermlineSingleSample.papi_settings",
        {
            "preemptible_tries": 3,
            "agg_preemptible_tries": 3
        }
    )

    ds.remove_param("genome")


if __name__ == "__main__":

    from cirro.helpers.preprocess_dataset import PreprocessDataset

    ds = PreprocessDataset.from_running()
    setup_reference(ds)

    all_inputs = build_inputs(read_samplesheet(ds.params.get("samplesheet")), ds.params)

    # The samplesheet is only read here, so keep it out of the Cromwell options
    ds.remove_param("samplesheet")

    # Write out the complete set of inputs
    write_json("inputs.json", all_inputs)

    # Write out each individual sample
    for i, single_input in enumerate(all_inputs):
        write_json(f"inputs.{i}.json", single_input)

    setup_options(ds)
