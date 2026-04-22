#!/usr/bin/env python3

from cirro.helpers.preprocess_dataset import PreprocessDataset
import pandas as pd
from typing import List
import urllib.request
import json


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


# Params that the extra JSON input field must never override.
# Includes both Cirro framework params (input/outdir) and Cirro-computed values
# (compute_multiplier) that must not be replaced by user input.
_PROTECTED_PARAMS = frozenset({
    "input",              # manifest.csv — built by this script
    "outdir",             # output directory — assigned by Cirro
    "igenomes_base",      # iGenomes S3 base URL
    "vep_cache",          # VEP annotation cache S3 path
    "snpeff_cache",       # snpEff annotation cache S3 path
    "monochrome_logs",    # internal logging flag
    "compute_multiplier", # WES/WGS resource multiplier — computed from wes param
})

# These Cirro-specific params are consumed by process-compute.config (not the nf-core schema)
# and must survive the schema filter even though they are not in nextflow_schema.json.
_CIRRO_PASSTHROUGH_PARAMS = frozenset({"compute_multiplier"})

# Maximum allowed size for the extra_params_json textarea to prevent resource exhaustion.
_MAX_EXTRA_JSON_BYTES = 10_000


def apply_extra_json_params(ds: PreprocessDataset):
    """Parse the extra_params_json textarea and merge user-supplied params into the workflow.

    The value is expected to be a JSON object (dict). Each key-value pair is added to
    ds.params with overwrite=True so that users can tune params set earlier in the script.
    A subsequent call to filter_params_by_schema() will discard any keys not recognized by
    the active pipeline version, so this function does not need to validate param names.

    Protected params (input, outdir, compute_multiplier, etc.) are never overridden.
    """
    extra_json_str = (ds.params.get("extra_params_json") or "").strip()

    # Always remove this Cirro-only field — it is never a valid Nextflow parameter
    ds.remove_param("extra_params_json", force=True)

    if not extra_json_str:
        ds.logger.info("extra_params_json: no extra parameters provided")
        return

    # Reject oversized payloads before parsing to prevent resource exhaustion
    if len(extra_json_str) > _MAX_EXTRA_JSON_BYTES:
        ds.logger.warning(
            f"extra_params_json: payload is {len(extra_json_str):,} chars "
            f"(limit {_MAX_EXTRA_JSON_BYTES:,}), skipping extra params"
        )
        return

    ds.logger.info("extra_params_json: parsing user-supplied JSON parameters")

    try:
        extra_params = json.loads(extra_json_str)
    except json.JSONDecodeError as e:
        # Warn but do not abort — a malformed JSON string should not block the run
        ds.logger.warning(f"extra_params_json: JSON parse error ({e}), skipping extra params")
        return

    if not isinstance(extra_params, dict):
        ds.logger.warning(
            f"extra_params_json: expected a JSON object but got {type(extra_params).__name__}, skipping"
        )
        return

    ds.logger.info(f"extra_params_json: found {len(extra_params)} parameter(s)")

    applied, skipped = 0, 0
    for key, value in extra_params.items():
        # Enforce that keys are non-empty strings (JSON requires string keys, but guard anyway)
        if not isinstance(key, str) or not key.strip():
            ds.logger.warning(f"extra_params_json: skipping invalid key {key!r}")
            skipped += 1
            continue
        if key in _PROTECTED_PARAMS:
            ds.logger.warning(
                f"extra_params_json: skipping protected parameter '{key}' (cannot be overridden via extra JSON)"
            )
            skipped += 1
            continue
        ds.logger.info(f"extra_params_json: applying {key}={value!r}")
        ds.add_param(key, value, overwrite=True)
        applied += 1

    ds.logger.info(
        f"extra_params_json: applied {applied} parameter(s), skipped {skipped} protected/invalid parameter(s)"
    )


def filter_params_by_schema(ds: PreprocessDataset):
    """Remove any params not present in the nf-core/sarek nextflow_schema.json.

    Fetches the schema for the user-selected workflow version from GitHub, then
    removes every key in ds.params that is not listed as a recognized pipeline
    parameter. This achieves two things:
      1. Cirro-only housekeeping params (workflow_version, extra_params_json) are
         cleaned up automatically without requiring explicit remove_param calls.
      2. Any extra params supplied via the JSON textarea that are not valid for
         the chosen pipeline version are silently dropped, preventing Nextflow
         schema-validation errors at runtime.

    workflow_version is always removed here (even on schema fetch failure) since it
    must never reach Nextflow. Params in _CIRRO_PASSTHROUGH_PARAMS are preserved
    even if absent from the nf-core schema, because they are consumed by Cirro's
    process-compute.config rather than by the pipeline itself.
    """
    version = ds.params.get("workflow_version", "3.8.1")

    # Remove workflow_version unconditionally — it is a Cirro-only param and must
    # not reach Nextflow regardless of whether the schema fetch succeeds.
    ds.remove_param("workflow_version", force=True)

    url = f"https://raw.githubusercontent.com/nf-core/sarek/{version}/nextflow_schema.json"
    ds.logger.info(f"filter_params_by_schema: fetching schema for nf-core/sarek {version}")

    try:
        with urllib.request.urlopen(url) as response:
            schema = json.loads(response.read().decode())
    except Exception as e:
        # If the schema cannot be fetched, log a warning and skip filtering.
        # Note: without filtering, any keys injected via extra_params_json that
        # are not valid Nextflow params will cause pipeline schema-validation errors.
        ds.logger.warning(f"filter_params_by_schema: could not fetch schema — {e}")
        ds.logger.warning("filter_params_by_schema: skipping schema-based param filtering")
        return

    # Collect all parameter names defined across the schema's top-level $defs sections
    allowed = set()
    for section in schema.get("$defs", {}).values():
        allowed.update(section.get("properties", {}).keys())

    ds.logger.info(f"filter_params_by_schema: schema defines {len(allowed):,} recognized parameters")

    # Remove params not recognized by the schema, but preserve Cirro passthrough params
    # (e.g. compute_multiplier) that are used by process-compute.config, not the pipeline schema.
    removed = [
        key for key in list(ds.params.keys()) if key not in allowed and key not in _CIRRO_PASSTHROUGH_PARAMS
    ]
    for key in removed:
        ds.remove_param(key, force=True)

    if removed:
        ds.logger.info(f"filter_params_by_schema: removed {len(removed)} unrecognized param(s): {removed}")
    else:
        ds.logger.info("filter_params_by_schema: all params are recognized by the schema")

    ds.logger.info(f"filter_params_by_schema: {len(ds.params)} param(s) remain after filtering")


if __name__ == "__main__":

    # Load the information for this dataset
    ds = PreprocessDataset.from_running()

    ds.logger.info(f"Starting sarek preprocess — workflow_version={ds.params.get('workflow_version', '3.8.1')!r}")
    ds.logger.info(f"analysis_type={ds.params.get('analysis_type')!r}, genome={ds.params.get('genome')!r}")

    # Make the samplesheet
    manifest = make_manifest(ds)

    # Log the manifest
    ds.logger.info(manifest.to_csv(index=None))

    # Write manifest
    manifest.to_csv("manifest.csv", index=None)
    ds.logger.info(f"Wrote {manifest.shape[0]} row(s) to manifest.csv")

    # JSON parameters !! revert after fix
    params = ds.params

    tools = params.get('tools')
    # Guard against None and empty list — both mean no variant caller was selected
    assert tools, "ERROR: You must select at least one variant calling tool."

    # Annotation tool is allowed to be empty, init empty list if it is
    # not 100% sure of annotation tool behavior if neither are selected i.e empty
    annotation_tool = params.get('annotation_tool') or []

    # Combine the two
    tools = ','.join(map(str, tools + annotation_tool))
    ds.logger.info(f"Combined tools string: {tools!r}")

    ds.add_param('tools', tools, overwrite=True)

    genome = params.get('genome')
    ds.logger.info(f"Reference genome: {genome!r}")

    # annotation_tool is a Cirro UI staging param — always remove it.
    # The combined tools string already incorporates annotation tool selections.
    ds.remove_param('annotation_tool', force=True)

    # construct logic for dbNSFP & SpliceAI
    # note reference genome selected
    # if true, construct the appropriate parameters as file paths.
    dbnsfp_param = params.get('vep_dbnsfp') # true or null
    spliceai = params.get('vep_spliceai') # true or null
    # Maps iGenomes genome identifiers to the (directory, ucsc) naming used in VEP reference paths
    database = {'GATK.GRCh37': ['GRCh37', 'hg19'],
                'GATK.GRCh38': ['GRCh38', 'hg38']}

    # dbNSFP — only available for human genomes
    if dbnsfp_param:
        if genome in database:
            dbnsfp = f"s3://pubweb-references/VEP/{database[genome][0]}/dbNSFP4.2a_{database[genome][0].lower()}.gz"
            dbnsfp_tbi = f"s3://pubweb-references/VEP/{database[genome][0]}/dbNSFP4.2a_{database[genome][0].lower()}.gz.tbi"
            ds.add_param('dbnsfp', dbnsfp, overwrite=True)
            ds.add_param('dbnsfp_tbi', dbnsfp_tbi, overwrite=True)
            ds.add_param('dbnsfp_consequence', 'ALL', overwrite=True)
            ds.logger.info(f"dbNSFP: enabled for {genome}")
        else:
            ds.logger.warning(f"No dbNSFP reference genome found for {genome} -- removing dbNSFP plugin")
            ds.remove_param('vep_dbnsfp')

    # SpliceAI — only available for human genomes; GRCm38 (mouse) has no reference data
    if spliceai:
        if genome not in database:
            ds.logger.warning(f"SpliceAI: no reference data available for genome {genome!r} -- skipping SpliceAI plugin")
        else:
            spliceai_snv = f"s3://pubweb-references/VEP/{database[genome][0]}/spliceai_scores.raw.snv.{database[genome][1]}.vcf.gz"
            spliceai_snv_tbi = f"s3://pubweb-references/VEP/{database[genome][0]}/spliceai_scores.raw.snv.{database[genome][1]}.vcf.gz.tbi"
            spliceai_indel = f"s3://pubweb-references/VEP/{database[genome][0]}/spliceai_scores.raw.indel.{database[genome][1]}.vcf.gz"
            spliceai_indel_tbi = f"s3://pubweb-references/VEP/{database[genome][0]}/spliceai_scores.raw.indel.{database[genome][1]}.vcf.gz.tbi"
            ds.add_param('spliceai_snv', spliceai_snv, overwrite=True)
            ds.add_param('spliceai_snv_tbi', spliceai_snv_tbi, overwrite=True)
            ds.add_param('spliceai_indel', spliceai_indel, overwrite=True)
            ds.add_param('spliceai_indel_tbi', spliceai_indel_tbi, overwrite=True)
            ds.logger.info(f"SpliceAI: enabled for {genome}")

    # PON handling

    if params.get('analysis_type') == 'Somatic Variant Calling':
        if genome == 'GATK.GRCh37':
            pon = "s3://pubweb-references/igenomes/Homo_sapiens/GATK/GRCh37/Annotation/GATKBundle/Mutect2-WGS-panel-b37.vcf.gz"
            pon_tbi = "s3://pubweb-references/igenomes/Homo_sapiens/GATK/GRCh37/Annotation/GATKBundle/Mutect2-WGS-panel-b37.vcf.gz.tbi"
            ds.add_param('pon', pon, overwrite=True)
            ds.add_param('pon_tbi', pon_tbi, overwrite=True)
            ds.logger.info("PON: added GRCh37 somatic panel of normals")
        if genome == "GATK.GRCh38":
            pon = "s3://pubweb-references/igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000g_pon.hg38.vcf.gz"
            pon_tbi = "s3://pubweb-references/igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000g_pon.hg38.vcf.gz.tbi"
            ds.add_param('pon', pon, overwrite=True)
            ds.add_param('pon_tbi', pon_tbi, overwrite=True)
            ds.logger.info("PON: added GRCh38 somatic panel of normals")

    ds.remove_param('analysis_type')

    # add index for germline_resource if present
    germline_resource = ds.params.get('germline_resource')
    if germline_resource:
        germline_resource_tbi = germline_resource + '.tbi'
        ds.add_param('germline_resource_tbi', germline_resource_tbi, overwrite=True)
        ds.logger.info("germline_resource: added .tbi index path")

    # Dynamic Resource Usage

    # `compute_multiplier` == 2 for WGS and 1 for WES; consumed by process-compute.config
    wes = params.get('wes', False)
    compute_multiplier = int(2 - int(wes))
    ds.add_param("compute_multiplier", compute_multiplier)
    ds.logger.info(f"compute_multiplier={compute_multiplier} ({'WES' if wes else 'WGS'})")

    ds.remove_param('wes')

    # If an intervals file was not selected, use --no_intervals
    if not ds.params.get("intervals"):
        ds.add_param("no_intervals", True)
        ds.logger.info("No intervals file selected — adding --no_intervals flag")

    # Fix the UMI form dependencies
    fix_umi_form_dependencies(ds)

    # Parse and apply any extra parameters the user provided as raw JSON.
    # These are added before schema filtering so that invalid keys are
    # automatically removed in the next step rather than causing Nextflow errors.
    apply_extra_json_params(ds)

    # Remove any parameters not defined in the nf-core/sarek nextflow_schema.json for
    # the selected workflow version. This also cleans up Cirro-only housekeeping params
    # (workflow_version, extra_params_json) that must not reach Nextflow.
    filter_params_by_schema(ds)

    # log all params
    ds.logger.info(f"Final params ({len(ds.params)} total): {ds.params}")
