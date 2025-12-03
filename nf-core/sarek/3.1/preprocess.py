#!/usr/bin/env python3

from cirro.helpers.preprocess_dataset import PreprocessDataset
import pandas as pd
from typing import List


def make_manifest(ds: PreprocessDataset) -> pd.DataFrame:

    # Filter out any index files that may have been uploaded
    ds.files = ds.files.loc[
        ds.files.apply(
            lambda r: r.get('readType', 'R') == 'R',
            axis=1
        )
    ]

    # Make a wide manifest
    manifest = ds.wide_samplesheet(
        index=["sampleIndex", "sample", "lane", "dataset"],
        columns="read",
        values="file",
        column_prefix="fastq_"
    )
    assert manifest.shape[0] > 0, "No files detected -- there may be an error with data ingest"

    # append metadata to file paths
    samples = ds.samplesheet.set_index("sample")
    manifest = manifest.assign(**{k: manifest["sample"].apply(v.get) for k, v in samples.items()})
    ordering = ['patient', 'sex', 'status', 'sample', 'lane', 'fastq_1', 'fastq_2']
    manifest = manifest.reindex(columns=ordering)

    # Overwrite the 'lane' column to provide a unique value per-row
    # This is necessary to account for datasets which merge multiple flowcells
    manifest = manifest.assign(lane=[
        str(i)
        for i in range(manifest.shape[0])
    ])

    # Run sanity checks before writing to manifest.csv
    expected_columns = ['patient', 'sex', 'status', 'sample', 'lane', 'fastq_1', 'fastq_2']
    manifest_columns = manifest.axes[1]
    non_canonical = [i for i in manifest_columns if i not in expected_columns]

    # Drop additional columns cleanly without exit
    if len(non_canonical) != 0:
        manifest = manifest.drop(columns=non_canonical)

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

        # Do some checking on each row
        for i, row in manifest.iterrows():

            entry_msg = f"Offending entry:\n {manifest.iloc[i].to_frame().T}"

            # Check that status is 0/1
            msg = "The status must contain 0 (normal) and/or 1 (tumor).\n"
            assert row["status"] in [0, 1], msg + entry_msg

            # Check that sex is XX/XY
            msg = "The column sex must consist of XX, XY or NA.\n"
            assert row["sex"] in ['XX', 'XY', 'NA'], msg + entry_msg

            # Check that 'patient' 'sample' and 'lane' are unique
            patient = str(row['patient'])
            sample = str(row['sample'])
            lane = str(row['lane'])
            msg = "Values for patient, sample and lane must be unique.\n "
            assert patient != sample and sample != lane, msg + entry_msg

    return manifest


def delete_params(ds: PreprocessDataset, param_list: List[str]):
    """Delete multiple parameters."""
    for param in param_list:
        if param in ds.params:
            ds.remove_param(param, force=True)


def fix_umi_form_dependencies(ds: PreprocessDataset):
    """
    Make sure that the UMI processing options are only provided to the workflow
    if the appropriate section was enabled in the form. This is necessary because
    the form data is still populated even if the element is hidden by the
    deselection of the element driving its dependency.
    """

    # Group the params based on what approach the user selected
    fgbio_params = ["umi_read_structure", "group_by_umi_strategy"]
    fastp_params = ["umi_location", "umi_length", "umi_base_skip"]

    # Get the option selected by the user
    umi_tool = ds.params.get("umi_tool", "none")

    # If no UMI processing was selected
    if umi_tool == "none":
        ds.logger.info("UMIs - user selected no UMI processing")
        delete_params(ds, fgbio_params + fastp_params)

    # Use fgbio for UMI processing
    elif umi_tool == "use_fgbio":
        ds.logger.info("UMIs - use fgbio")
        delete_params(ds, fastp_params)

    # Use fastp for UMI processing
    elif umi_tool == "use_fastp":
        ds.logger.info("UMIs - use fastp")
        delete_params(ds, fgbio_params)

    else:
        raise Exception(f"Did not expect umi_tool selection: {umi_tool}")
    
    ds.remove_param("umi_tool", force=True)


if __name__ == "__main__":

    # Load the information for this dataset
    ds = PreprocessDataset.from_running()

    # Make the samplesheet
    manifest = make_manifest(ds)

    # Log the manifest
    ds.logger.info(manifest.to_csv(index=None))

    # Write manifest
    manifest.to_csv("manifest.csv", index=None)

    # JSON parameters !! revert after fix
    params = ds.params

    tools = params.get('tools')
    assert tools is not None, "ERROR: You must select a variant calling tool."

    # Annotation tool is allowed to be empty, init empty list if it is
    # not 100% sure of annotation tool behavior if neither are selected i.e empty
    annotation_tool = params.get('annotation_tool') or []

    # Combine the two
    tools = ','.join(map(str, tools + annotation_tool))

    ds.add_param('tools', tools, overwrite=True)

    genome = params.get('genome')

    # if user does not select VEP/snpEff then annotation tool param does not exist.
    # script sets it as empty list, use this to toggle deleting the param to avoid error.
    if len(annotation_tool) != 0:
        ds.remove_param('annotation_tool')

    # construct logic for dbNSFP & SpliceAI
    # note reference genome selected
    # if true, construct the appropriate parameters as file paths. 
    dbnsfp_param = params.get('vep_dbnsfp') # true or null 
    spliceai = params.get('vep_spliceai') # true or null 
    database = {'GATK.GRCh37': ['GRCh37', 'hg19'],
                'GATK.GRCh38': ['GRCh38', 'hg38']}

    # dbNSFP
    if dbnsfp_param:
        dbnsfp = f"s3://pubweb-references/VEP/{database[genome][0]}/dbNSFP4.2a_{database[genome][0].lower()}.gz"
        dbnsfp_tbi = f"s3://pubweb-references/VEP/{database[genome][0]}/dbNSFP4.2a_{database[genome][0].lower()}.gz.tbi"
        ds.add_param('dbnsfp', dbnsfp, overwrite=True)
        ds.add_param('dbnsfp_tbi', dbnsfp_tbi, overwrite=True)
        ds.add_param('dbnsfp_consequence', 'ALL', overwrite=True)

    # Splice AI
    if spliceai:
        spliceai_snv = f"s3://pubweb-references/VEP/{database[genome][0]}/spliceai_scores.raw.snv.{database[genome][1]}.vcf.gz"
        spliceai_snv_tbi = f"s3://pubweb-references/VEP/{database[genome][0]}/spliceai_scores.raw.snv.{database[genome][1]}.vcf.gz.tbi"
        spliceai_indel = f"s3://pubweb-references/VEP/{database[genome][0]}/spliceai_scores.raw.indel.{database[genome][1]}.vcf.gz"
        spliceai_indel_tbi = f"s3://pubweb-references/VEP/{database[genome][0]}/spliceai_scores.raw.indel.{database[genome][1]}.vcf.gz.tbi"
        ds.add_param('spliceai_snv', spliceai_snv, overwrite=True)
        ds.add_param('spliceai_snv_tbi', spliceai_snv_tbi, overwrite=True)
        ds.add_param('spliceai_indel', spliceai_indel, overwrite=True)
        ds.add_param('spliceai_indel_tbi', spliceai_indel_tbi, overwrite=True)

    # PON handling

    if params.get('analysis_type') == 'Somatic Variant Calling':
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
        int(2 - int(params.get('wes')))
    )

    ds.remove_param('wes')

    # If an intervals file was not selected, use --no_intervals
    if not ds.params.get("intervals"):
        ds.add_param("no_intervals", True)

    # Fix the UMI form dependencies
    fix_umi_form_dependencies(ds)

    # log all params
    ds.logger.info(ds.params)
