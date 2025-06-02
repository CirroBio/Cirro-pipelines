#!/usr/bin/env python3

import os
from cirro.helpers.preprocess_dataset import PreprocessDataset
import pandas as pd
from cirro.api.models.s3_path import S3Path
import boto3
import json


def make_manifest(ds: PreprocessDataset) -> pd.DataFrame:

    ds.logger.info("Input Files:")
    ds.logger.info(ds.files.to_csv(index=None))

    assert ds.files.shape[0] > 0, "No files detected -- error with data ingest"

    # Format a wide sample sheet
    manifest = (
        ds.files
        .assign(ext=ds.files["file"].apply(lambda s: s.split(".")[-1]))
        .pivot(
            index="sample",
            columns="ext",
            values="file"
        )
    )

    assert manifest.shape[0] > 0, "No files detected -- error with data ingest"

    ds.logger.info("Pivoted Table:")
    ds.logger.info(manifest.to_csv())

    # Each sample should have both a bam and bai file
    for ext in ["bam", "bai"]:
        assert ext in manifest.columns.values, f"Requires {ext} files"
    for sample, r in manifest.iterrows():
        for ext in ["bam", "bai"]:
            assert not pd.isnull(r[ext]), f"Missing {ext} for {sample}"

    # Reset the index to get the sample column back
    manifest = manifest.reset_index()

    # append metadata to file paths
    samples = ds.samplesheet.set_index("sample")
    manifest = manifest.assign(**{
        k: manifest["sample"].apply(v.get)
        for k, v in samples.items()
    })

    ordering = ['patient', 'sex', 'status', 'sample', 'lane', 'bam', 'bai']
    manifest: pd.DataFrame = manifest.reindex(columns=ordering)

    # Overwrite the 'lane' column to provide a unique value per-row
    # This is necessary to account for datasets which merge multiple flowcells
    manifest = manifest.assign(lane=[
        str(i)
        for i in range(manifest.shape[0])
    ])

    # Set the default value for 'sex' to be NA
    manifest = manifest.assign(
        sex=manifest['sex'].fillna('NA')
    )

    # Transform status values "Normal" -> 0 and "Tumor" -> 1
    manifest = manifest.replace(
        to_replace=dict(
            status=dict(
                normal=0,
                Normal=0,
                tumor=1,
                Tumor=1
            )
        )
    )

    # Set the default status to 0
    manifest = manifest.assign(
        status=manifest["status"].fillna(0).apply(int)
    )

    # If the user selected Germline Variant Calling
    if ds.params["analysis_type"] == "Germline Variant Calling":

        # Set the default "patient" attribute as the sample
        manifest = manifest.assign(
            patient=manifest["patient"].fillna(manifest["sample"])
        )

    # If the user selected Somatic Variant Calling
    else:

        # Run checks on each row of the manifest
        for i, row in manifest.iterrows():

            line_msg = f"\nOffending entry:\n {manifest.iloc[i].to_frame().T}"

            # Check that status is 0/1
            msg = "The column 'status' must be 0 (normal) and/or 1 (tumor)."
            msg = msg + line_msg
            assert row["status"] in [0, 1], msg

            # Check that sex is XX/XY
            msg = "ERROR: The column sex must consist of XX, XY or NA."
            msg = msg + line_msg
            assert row["sex"] in ['XX', 'XY', 'NA'], msg

            # Check that 'patient' 'sample' and 'lane' are unique
            patient = str(row['patient'])
            sample = str(row['sample'])
            lane = str(row['lane'])
            msg = "ERROR: patient, sample and lane must be unique."
            msg = msg + line_msg
            assert patient != sample and sample != lane, msg

    return manifest

if __name__ == "__main__":

    # Load the information for this dataset
    ds = PreprocessDataset.from_running()

    # Make the samplesheet
    manifest = make_manifest(ds)

    # Log the manifest
    ds.logger.info(manifest.to_csv(index=None))

    # Write manifest
    manifest.to_csv("manifest.csv", index=None)

    tools = ds.params.get('tools')
    assert tools is not None, "ERROR: You must select a variant calling tool."

    # Annotation tool is allowed to be empty, init empty list if it is
    annotation_tool = ds.params.get('annotation_tool') or []

    # Combine the two
    tools = ','.join(map(str, tools + annotation_tool))

    ds.add_param('tools', tools, overwrite=True)

    # If an intervals file was not selected, use --no_intervals
    if not ds.params.get("intervals"):
        ds.add_param("no_intervals", True)

    # vep_cache_version (attempt to resolve 'expected numeric, got string' bug).
    # this gets updated regularly - keep an eye on it / need to think
    # how to automatically update this. 
    genome = ds.params.get('genome')
    cache_key = {'GATK.GRCh37': 106, 'GATK.GRCh38': 106, 'GRCm38': 102}
    vep_cache_version = cache_key[genome]
    ds.add_param('vep_cache_version', vep_cache_version, overwrite=True)

    # if user does not select VEP/snpEff then annotation tool param does not exist.
    # script sets it as empty list, use this to toggle deleting the param to avoid error.
    if len(annotation_tool) != 0:
        ds.remove_param('annotation_tool')

    # construct logic for dbNSFP & SpliceAI
    # note reference genome selected
    # if true, construct the appropriate parameters as file paths.
    dbnsfp_param = ds.params.get('vep_dbnsfp')  # true or null
    spliceai = ds.params.get('vep_spliceai')  # true or null
    database = {'GATK.GRCh37': ['GRCh37', 'hg19'],
                'GATK.GRCh38': ['GRCh38', 'hg38']}

    # dbNSFP
    ref_prefix = f"s3://pubweb-references/VEP/{database[genome][0]}"
    if dbnsfp_param:
        dbnsfp = f"{ref_prefix}/dbNSFP4.2a_{database[genome][0].lower()}.gz"
        dbnsfp_tbi = f"{ref_prefix}/dbNSFP4.2a_{database[genome][0].lower()}.gz.tbi"
        ds.add_param('dbnsfp', dbnsfp, overwrite=True)
        ds.add_param('dbnsfp_tbi', dbnsfp_tbi, overwrite=True)
        ds.add_param('dbnsfp_consequence', 'ALL', overwrite=True)

    # Splice AI
    if spliceai:
        spliceai_snv = f"{ref_prefix}/spliceai_scores.raw.snv.{database[genome][1]}.vcf.gz"
        spliceai_snv_tbi = f"{ref_prefix}/spliceai_scores.raw.snv.{database[genome][1]}.vcf.gz.tbi"
        spliceai_indel = f"{ref_prefix}/spliceai_scores.raw.indel.{database[genome][1]}.vcf.gz"
        spliceai_indel_tbi = f"{ref_prefix}/spliceai_scores.raw.indel.{database[genome][1]}.vcf.gz.tbi"
        ds.add_param('spliceai_snv', spliceai_snv, overwrite=True)
        ds.add_param('spliceai_snv_tbi', spliceai_snv_tbi, overwrite=True)
        ds.add_param('spliceai_indel', spliceai_indel, overwrite=True)
        ds.add_param('spliceai_indel_tbi', spliceai_indel_tbi, overwrite=True)

    # workflow is failing due to 'Unknown config attribute `params.vep_version` -- check config file: /root/.nextflow/assets/nf-core/sarek/nextflow.config'
    # test by hard-coding here
    ds.add_param('vep_version', '106.1', overwrite=True)

    # PON handling

    if ds.params.get('analysis_type') == 'Somatic Variant Calling':
        if genome == 'GATK.GRCh37':
            pon = "s3://pubweb-references/igenomes/Homo_sapiens/GATK/GRCh37/Annotation/GATKBundle/Mutect2-WGS-panel-b37.vcf.gz"
            pon_tbi = "s3://pubweb-references/igenomes/Homo_sapiens/GATK/GRCh37/Annotation/GATKBundle/Mutect2-WGS-panel-b37.vcf.gz.tbi"
            ds.add_param('pon', pon, overwrite=True)
            ds.add_param('pon_tbi', pon_tbi, overwrite=True)
        if genome == "GATK.GRCh38":
            pon = "s3://pubweb-references/igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000g_pon.hg38.vcf.gz"
            pon_tbi = "s3://pubweb-references/igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000g_pon.hg38.vcf.gz.tbi"
            ds.add_param('pon', pon, overwrite=True)
            ds.add_param('pon_tbi', pon_tbi, overwrite=True)    

    ds.remove_param('analysis_type')

    # add index for germline_resource if present
    germline_resource = ds.params.get('germline_resource')
    if germline_resource:
        germline_resource_tbi = germline_resource + '.tbi'
        ds.add_param('germline_resource_tbi', germline_resource_tbi, overwrite=True)

    # Dynamic Resource Usage

    # `compute_multiplier` == 2 for WGS and 1 for WES
    ds.add_param(
        "compute_multiplier",
        int(2 - int(ds.params.get('wes')))
    )

    ds.remove_param('wes')

    # log all params
    ds.logger.info(ds.params)
