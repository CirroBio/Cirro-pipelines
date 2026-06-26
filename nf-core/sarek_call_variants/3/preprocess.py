#!/usr/bin/env python3

from cirro.helpers.preprocess_dataset import PreprocessDataset
import pandas as pd
import urllib.request
import urllib.error
import json


def make_manifest(ds: PreprocessDataset) -> pd.DataFrame:

    ds.logger.info("Input Files:")
    ds.logger.info(ds.files.to_csv(index=None))

    assert ds.files.shape[0] > 0, "No files detected -- error with data ingest"

    # Format a wide sample sheet
    def _bam_pref(filename: str) -> int:
        if filename.endswith('.recal.bam'):
            return 0
        elif filename.endswith('.sorted.bam'):
            return 1
        else:
            return 2

    manifest = (
        ds.files
        .loc[ds.files["file"].str.contains(r'\.bam(\.bai)?$', regex=True)]
        .assign(
            stem=lambda df: df["file"].str.replace(r'\.bam(\.bai)?$', '', regex=True),
            ext=lambda df: df["file"].str.extract(r'\.(bam(?:\.bai)?)$')[0]
        )
        .pivot(index=["sample", "stem"], columns="ext", values="file")
        .rename_axis(columns=None)
        .dropna(subset=["bam", "bam.bai"])
        .assign(_pref=lambda df: df["bam"].apply(_bam_pref))
        .groupby(level="sample", group_keys=False)
        .apply(lambda g: g[g["_pref"] == g["_pref"].min()])
        .drop(columns=["_pref"])
        .reset_index()
        .rename(columns={"bam.bai": "bai"})
        .drop(columns=["stem"])
    )

    assert manifest.shape[0] > 0, "No files detected -- error with data ingest"

    ds.logger.info("BAM/BAI pairs:")
    ds.logger.info(manifest.to_csv(index=None))

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

# Params set by Cirro infrastructure or computed by this script that must not
# be overridden by user-supplied extra JSON.
_PROTECTED_PARAMS = frozenset({
    "input",
    "outdir",
    "igenomes_base",
    "vep_cache",
    "snpeff_cache",
    "monochrome_logs",
    "step",             # hardcoded to "variant_calling" by process-input.json
    "compute_multiplier",  # computed from wes
    "wes",              # consumed to compute compute_multiplier before extra JSON is applied
    "intervals",        # consumed to set no_intervals before extra JSON is applied
})

# compute_multiplier and optical_duplicate_pixel_distance are consumed by process-compute.config
# (not the nf-core schema). igenomes_base, vep_cache, snpeff_cache, and monochrome_logs are
# Cirro-pinned S3 paths that must survive schema filtering regardless of sarek version.
_CIRRO_PASSTHROUGH_PARAMS = frozenset({
    "compute_multiplier",
    "optical_duplicate_pixel_distance",
    "igenomes_base",
    "vep_cache",
    "snpeff_cache",
    "monochrome_logs",
})

_MAX_EXTRA_JSON_BYTES = 10_000


def apply_extra_json_params(ds: PreprocessDataset):
    """Parse the extra_params_json textarea and merge user-supplied params into the workflow.

    Applied before filter_params_by_schema so that unrecognized keys are automatically
    removed in the subsequent schema validation step.
    """
    extra_json_str = (ds.params.get("extra_params_json") or "").strip()
    ds.remove_param("extra_params_json", force=True)

    if not extra_json_str:
        ds.logger.info("extra_params_json: no extra parameters provided")
        return

    if len(extra_json_str) > _MAX_EXTRA_JSON_BYTES:
        ds.logger.warning(
            f"extra_params_json: payload too large ({len(extra_json_str):,} chars), skipping"
        )
        return

    ds.logger.info("extra_params_json: parsing user-supplied JSON parameters")

    try:
        extra_params = json.loads(extra_json_str)
    except json.JSONDecodeError as e:
        ds.logger.warning(f"extra_params_json: JSON parse error ({e}), skipping")
        return

    if not isinstance(extra_params, dict):
        ds.logger.warning(
            f"extra_params_json: expected a JSON object but got {type(extra_params).__name__}, skipping"
        )
        return

    ds.logger.info(f"extra_params_json: found {len(extra_params)} parameter(s)")

    applied, skipped = 0, 0
    for key, value in extra_params.items():
        if not isinstance(key, str) or not key.strip():
            ds.logger.warning(f"extra_params_json: skipping invalid key {key!r}")
            skipped += 1
            continue
        if key in _PROTECTED_PARAMS:
            ds.logger.warning(f"extra_params_json: skipping protected parameter '{key}'")
            skipped += 1
            continue
        ds.logger.info(f"extra_params_json: applying {key}={value!r}")
        ds.add_param(key, value, overwrite=True)
        applied += 1

    ds.logger.info(
        f"extra_params_json: applied {applied} parameter(s), skipped {skipped} protected/invalid"
    )


def filter_params_by_schema(ds: PreprocessDataset):
    """Remove any params not present in the nf-core/sarek nextflow_schema.json.

    Fetches the schema for the user-selected version from GitHub and removes every
    ds.params key not listed there. Cleans up Cirro-only params (workflow_version,
    extra_params_json) and drops invalid extra JSON params before they reach Nextflow.
    Params in _CIRRO_PASSTHROUGH_PARAMS are preserved even if absent from the schema.
    """
    version = ds.params.get("workflow_version", _DEFAULT_WORKFLOW_VERSION)

    # Always remove — this is a Cirro-only param that must not reach Nextflow
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

    allowed = set()
    for section in {**schema.get("$defs", {}), **schema.get("definitions", {})}.values():
        allowed.update(section.get("properties", {}).keys())

    ds.logger.info(f"filter_params_by_schema: schema defines {len(allowed):,} parameters")

    removed = [
        key for key in list(ds.params.keys()) if key not in allowed and key not in _CIRRO_PASSTHROUGH_PARAMS
    ]
    for key in removed:
        ds.remove_param(key, force=True)

    if removed:
        ds.logger.info(f"filter_params_by_schema: removed {len(removed)} unrecognized param(s): {removed}")
    ds.logger.info(f"filter_params_by_schema: {len(ds.params)} param(s) remain")


if __name__ == "__main__":

    ds = PreprocessDataset.from_running()

    ds.logger.info(f"Starting sarek_call_variants preprocess — workflow_version={ds.params.get('workflow_version', _DEFAULT_WORKFLOW_VERSION)!r}")
    ds.logger.info(f"analysis_type={ds.params.get('analysis_type')!r}, genome={ds.params.get('genome')!r}")

    manifest = make_manifest(ds)
    ds.logger.info(manifest.to_csv(index=None))
    manifest.to_csv("manifest.csv", index=None)
    ds.logger.info(f"Wrote {manifest.shape[0]} row(s) to manifest.csv")

    # Validate tool/sample/resource dependencies while tools is still a list.
    validate_tool_dependencies(ds, manifest)

    tools = ds.params.get('tools')
    assert tools, "ERROR: You must select at least one variant calling tool."

    annotation_tool = ds.params.get('annotation_tool') or []

    tools = ','.join(map(str, tools + annotation_tool))
    ds.logger.info(f"Combined tools string: {tools!r}")

    ds.add_param('tools', tools, overwrite=True)

    # If an intervals file was not selected, use --no_intervals
    if not ds.params.get("intervals"):
        ds.add_param("no_intervals", True)
        ds.logger.info("No intervals file selected — adding --no_intervals flag")

    # Resolve the reference genome (iGenomes vs Custom BWA index) before any
    # downstream logic reads the genome param. Warn that custom genomes lack the
    # annotation/known-sites references the iGenomes path supplies.
    if ds.params.get("genome_source") == "dataset" and ds.params.get("annotation_tool"):
        ds.logger.warning(
            "Custom genome selected: variant annotation (VEP/snpEff) reference data is not "
            "available for custom genomes and annotation will be skipped or fail."
        )
    resolve_reference_genome(ds)

    genome = ds.params.get('genome')
    ds.logger.info(f"Reference genome: {genome!r}")

    # annotation_tool is a Cirro UI staging param — always remove it
    ds.remove_param('annotation_tool', force=True)

    # Maps iGenomes genome identifiers to the (directory, ucsc) naming used in VEP reference paths
    database = {'GATK.GRCh37': ['GRCh37', 'hg19'],
                'GATK.GRCh38': ['GRCh38', 'hg38']}

    dbnsfp_param = ds.params.get('vep_dbnsfp')  # true or null
    spliceai = ds.params.get('vep_spliceai')  # true or null

    # dbNSFP — only available for human genomes
    if dbnsfp_param:
        if genome in database:
            ref_prefix = f"s3://pubweb-references/VEP/{database[genome][0]}"
            dbnsfp = f"{ref_prefix}/dbNSFP4.2a_{database[genome][0].lower()}.gz"
            dbnsfp_tbi = f"{ref_prefix}/dbNSFP4.2a_{database[genome][0].lower()}.gz.tbi"
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
            ref_prefix = f"s3://pubweb-references/VEP/{database[genome][0]}"
            spliceai_snv = f"{ref_prefix}/spliceai_scores.raw.snv.{database[genome][1]}.vcf.gz"
            spliceai_snv_tbi = f"{ref_prefix}/spliceai_scores.raw.snv.{database[genome][1]}.vcf.gz.tbi"
            spliceai_indel = f"{ref_prefix}/spliceai_scores.raw.indel.{database[genome][1]}.vcf.gz"
            spliceai_indel_tbi = f"{ref_prefix}/spliceai_scores.raw.indel.{database[genome][1]}.vcf.gz.tbi"
            ds.add_param('spliceai_snv', spliceai_snv, overwrite=True)
            ds.add_param('spliceai_snv_tbi', spliceai_snv_tbi, overwrite=True)
            ds.add_param('spliceai_indel', spliceai_indel, overwrite=True)
            ds.add_param('spliceai_indel_tbi', spliceai_indel_tbi, overwrite=True)
            ds.logger.info(f"SpliceAI: enabled for {genome}")

    # PON handling
    if ds.params.get('analysis_type') == 'Somatic Variant Calling':
        # A user-supplied panel-of-normals takes precedence; sarek generates its
        # .tbi index automatically if one is not alongside it.
        if ds.params.get('pon'):
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

    # `compute_multiplier` == 2 for WGS and 1 for WES; consumed by process-compute.config
    wes = ds.params.get('wes', False)
    compute_multiplier = int(2 - int(wes))
    ds.add_param("compute_multiplier", compute_multiplier)
    ds.logger.info(f"compute_multiplier={compute_multiplier} ({'WES' if wes else 'WGS'})")

    ds.remove_param('wes')

    apply_extra_json_params(ds)
    filter_params_by_schema(ds)

    ds.logger.info(f"Final params ({len(ds.params)} total): {ds.params}")
