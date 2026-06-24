#!/usr/bin/env python3

from cirro.helpers.preprocess_dataset import PreprocessDataset
import pandas as pd
from typing import List
import urllib.request
import urllib.error
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

    # fastp UMI extraction requires a UMI length and is incompatible with reading
    # UMIs from the read header (mirrors nf-core/sarek's runtime checks).
    if umi_tool == "use_fastp":
        if not ds.params.get("umi_length"):
            ds.logger.error(
                "fastp UMI extraction (UMI Location) requires a UMI Length to be specified."
            )
        if ds.params.get("umi_in_read_header"):
            ds.logger.error(
                "fastp UMI extraction and 'UMI in Read Header' cannot be used together. "
                "Please choose only one UMI handling approach."
            )

    # fgbio consensus with UMIs already in the read header requires the read
    # structure to be '+T +T'.
    if umi_tool == "use_fgbio" and ds.params.get("umi_in_read_header"):
        if ds.params.get("umi_read_structure") not in (None, "+T +T"):
            ds.logger.error(
                "When UMIs are in the read header (umi_in_read_header), the fgbio UMI "
                "read structure must be '+T +T'. Please update the UMI Read Structure."
            )


# Tool/sample/resource dependencies that the form cannot express, mirroring the
# runtime validation in nf-core/sarek's samplesheet_to_channel subworkflow. Sarek
# errors on the first group at runtime, so we fail fast here with a clear message.
_TUMOR_REQUIRING_TOOLS = ("ascat", "controlfreec", "mutect2", "msisensorpro")
_NORMAL_REQUIRING_TOOLS = ("ascat", "deepvariant", "haplotypecaller", "msisensorpro")
_SEX_REQUIRING_TOOLS = ("ascat", "controlfreec")


def validate_tool_dependencies(ds: PreprocessDataset, manifest: pd.DataFrame):
    """Check tool/sample/resource dependencies before the run starts.

    Mirrors the runtime validation in nf-core/sarek's samplesheet_to_channel
    subworkflow. All issues are logged (errors for conditions sarek treats as
    fatal, warnings for softer ones) rather than raised, so the user sees the
    problem in the preprocess logs.
    """
    selected = {str(t) for t in (ds.params.get("tools") or [])}
    if not selected:
        return

    statuses = set(manifest["status"])

    tumor_tools = sorted(selected.intersection(_TUMOR_REQUIRING_TOOLS))
    if tumor_tools and 1 not in statuses:
        ds.logger.error(
            "These selected tool(s) require at least one tumor sample (status=Tumor), "
            f"but none was found in the samplesheet: {', '.join(tumor_tools)}."
        )

    normal_tools = sorted(selected.intersection(_NORMAL_REQUIRING_TOOLS))
    if normal_tools and 0 not in statuses:
        ds.logger.error(
            "These selected tool(s) require at least one normal sample (status=Normal), "
            f"but none was found in the samplesheet: {', '.join(normal_tools)}."
        )

    sex_tools = sorted(selected.intersection(_SEX_REQUIRING_TOOLS))
    if sex_tools and (manifest["sex"] == "NA").any():
        ds.logger.error(
            f"The tool(s) {', '.join(sex_tools)} require sex (XX/XY) for every sample, "
            "but one or more samples have sex=NA. Set the sex column in the samplesheet."
        )

    # ASCAT needs allele/loci reference files that Cirro does not supply by default.
    if "ascat" in selected and not (ds.params.get("ascat_alleles") and ds.params.get("ascat_loci")):
        ds.logger.warning(
            "ASCAT was selected but ascat_alleles/ascat_loci were not provided. ASCAT will "
            "fail unless both are supplied via the Extra Parameters (JSON) field."
        )

    # Mutect2 without a germline resource runs without germline-based filtering.
    if "mutect2" in selected and not ds.params.get("germline_resource"):
        ds.logger.warning(
            "Mutect2 was selected without a germline resource; no germline-based filtering "
            "of somatic calls will be performed."
        )


def warn_custom_genome_limitations(ds: PreprocessDataset):
    """Warn about features unavailable when a custom BWA genome dataset is used.

    Custom genome datasets provide only the FASTA + BWA index — no GATK known-sites
    (dbsnp/known_indels) or annotation caches — so base recalibration and variant
    annotation cannot run against them. Must be called before resolve_reference_genome
    removes ``genome_source``.
    """
    if ds.params.get("genome_source") != "dataset":
        return
    ds.logger.warning(
        "Custom genome selected: GATK known-sites (dbsnp/known_indels) are not available, "
        "so base recalibration cannot run (it is skipped automatically — see "
        "skip_baserecalibration_without_known_sites)."
    )
    if ds.params.get("annotation_tool"):
        ds.logger.warning(
            "Custom genome selected: variant annotation (VEP/snpEff) reference data is not "
            "available for custom genomes and annotation will be skipped or fail."
        )


def skip_baserecalibration_without_known_sites(ds: PreprocessDataset, is_custom_genome: bool):
    """Skip base recalibration when no known-sites resources are available.

    GATK BaseRecalibrator requires at least one of dbsnp/known_indels. iGenomes
    references supply these via the genome config at runtime, so only custom genomes
    need this guard: when neither resource is present in the params, add
    `baserecalibrator` to `skip_tools` so the run does not fail. Call after all params
    are populated (post extra-JSON and schema filter) so user-supplied resources are
    taken into account.
    """
    if not is_custom_genome:
        return
    if ds.params.get("dbsnp") or ds.params.get("known_indels"):
        return
    existing = ds.params.get("skip_tools")
    skip = [t for t in str(existing).split(",") if t] if existing else []
    if "baserecalibrator" not in skip:
        skip.append("baserecalibrator")
    ds.add_param("skip_tools", ",".join(skip), overwrite=True)
    ds.logger.info(
        "No dbsnp/known_indels provided — adding 'baserecalibrator' to skip_tools "
        f"(skip_tools={','.join(skip)})."
    )


def resolve_reference_genome(ds: PreprocessDataset):
    """Wire up the reference based on the iGenomes vs Custom Genome selection.

    For iGenomes the curated ``genome`` key is passed through unchanged. For a
    custom genome the user selects a pre-built BWA index dataset; we point
    ``--fasta``/``--bwa`` at that dataset and drop ``--genome``/``--igenomes_base``.
    The BWA index pipeline publishes ``genome.fasta.gz`` and the flat index files
    (``genome.{amb,ann,bwt,pac,sa}``) directly into the dataset's data directory,
    so the directory itself serves as the ``--bwa`` argument (nf-core's bwa/mem
    module derives the index prefix from the ``.amb`` file).
    """
    genome_source = ds.params.get("genome_source")
    ds.remove_param("genome_source", force=True)

    bwa_index = ds.params.get("bwa_index")
    ds.remove_param("bwa_index", force=True)

    if genome_source == "dataset":
        if not bwa_index:
            raise ValueError(
                "Custom Genome selected but no BWA genome index dataset was provided."
            )
        ds.logger.info(f"genome_source=dataset: using custom BWA index at {bwa_index}")
        ds.add_param("fasta", f"{bwa_index}/genome.fasta.gz", overwrite=True)
        ds.add_param("bwa", bwa_index, overwrite=True)
        ds.remove_param("genome", force=True)
        ds.remove_param("igenomes_base", force=True)
    else:
        ds.logger.info(f"genome_source=igenomes: genome={ds.params.get('genome')!r}")


_DEFAULT_WORKFLOW_VERSION = "3.8.1"

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
    "wes",               # consumed to compute compute_multiplier before extra JSON is applied
    "intervals",         # consumed to set no_intervals before extra JSON is applied
})

# These params must survive the schema filter regardless of the sarek version.
# compute_multiplier/optical_duplicate_pixel_distance are consumed by process-compute.config
# (not the nf-core schema). igenomes_base, vep_cache, snpeff_cache, and monochrome_logs are
# standard nf-core params present in modern sarek schemas, but we pin them here as well so
# that Cirro's custom S3 paths are never accidentally removed if an unusual schema version
# omits them.
_CIRRO_PASSTHROUGH_PARAMS = frozenset({
    "compute_multiplier",
    "optical_duplicate_pixel_distance",
    "igenomes_base",
    "vep_cache",
    "snpeff_cache",
    "monochrome_logs",
})

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
    version = ds.params.get("workflow_version", _DEFAULT_WORKFLOW_VERSION)

    # Remove workflow_version unconditionally — it is a Cirro-only param and must
    # not reach Nextflow regardless of whether the schema fetch succeeds.
    ds.remove_param("workflow_version", force=True)

    url = f"https://raw.githubusercontent.com/nf-core/sarek/{version}/nextflow_schema.json"
    ds.logger.info(f"filter_params_by_schema: fetching schema for nf-core/sarek {version}")

    try:
        with urllib.request.urlopen(url, timeout=30) as response:
            schema = json.loads(response.read().decode())
    except urllib.error.HTTPError as e:
        if e.code == 404:
            raise RuntimeError(
                f"nf-core/sarek version '{version}' not found. "
                f"See https://nf-co.re/sarek/releases for available versions."
            ) from e
        ds.logger.warning(f"filter_params_by_schema: HTTP error fetching schema — {e}, skipping filter")
        return
    except Exception as e:
        ds.logger.warning(f"filter_params_by_schema: could not fetch schema — {e}, skipping filter")
        return

    # Collect all parameter names defined across the schema's top-level $defs sections
    allowed = set()
    for section in {**schema.get("$defs", {}), **schema.get("definitions", {})}.values():
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

    ds.logger.info(f"Starting sarek preprocess — workflow_version={ds.params.get('workflow_version', _DEFAULT_WORKFLOW_VERSION)!r}")
    ds.logger.info(f"analysis_type={ds.params.get('analysis_type')!r}, genome={ds.params.get('genome')!r}")

    # Make the samplesheet
    manifest = make_manifest(ds)

    # Log the manifest
    ds.logger.info(manifest.to_csv(index=None))

    # Write manifest
    manifest.to_csv("manifest.csv", index=None)
    ds.logger.info(f"Wrote {manifest.shape[0]} row(s) to manifest.csv")

    # Validate tool/sample/resource dependencies and warn about custom-genome
    # limitations while genome_source/tools/annotation_tool are still present.
    validate_tool_dependencies(ds, manifest)
    warn_custom_genome_limitations(ds)

    # Capture the genome source before resolve_reference_genome removes it.
    is_custom_genome = ds.params.get("genome_source") == "dataset"

    # Resolve the reference genome (iGenomes vs Custom BWA index) before any
    # downstream logic reads the genome param.
    resolve_reference_genome(ds)

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
        # A user-supplied panel-of-normals takes precedence; sarek generates its
        # .tbi index automatically if one is not alongside it.
        if params.get('pon'):
            ds.logger.info("PON: using user-supplied panel of normals")
        elif genome == 'GATK.GRCh37':
            pon = "s3://pubweb-references/igenomes/Homo_sapiens/GATK/GRCh37/Annotation/GATKBundle/Mutect2-WGS-panel-b37.vcf.gz"
            pon_tbi = "s3://pubweb-references/igenomes/Homo_sapiens/GATK/GRCh37/Annotation/GATKBundle/Mutect2-WGS-panel-b37.vcf.gz.tbi"
            ds.add_param('pon', pon, overwrite=True)
            ds.add_param('pon_tbi', pon_tbi, overwrite=True)
            ds.logger.info("PON: added GRCh37 somatic panel of normals")
        elif genome == "GATK.GRCh38":
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

    # With all params populated, skip base recalibration for custom genomes that
    # lack known-sites resources (dbsnp/known_indels). iGenomes supplies these.
    skip_baserecalibration_without_known_sites(ds, is_custom_genome)

    # log all params
    ds.logger.info(f"Final params ({len(ds.params)} total): {ds.params}")
