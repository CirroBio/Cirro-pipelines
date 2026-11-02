#!/usr/bin/env python3

from __future__ import annotations

import json
import urllib.error
import urllib.request
from typing import TYPE_CHECKING

import pandas as pd

# The SDK is imported in __main__, so test_preprocess.py runs without it
if TYPE_CHECKING:
    from cirro.helpers.preprocess_dataset import PreprocessDataset


# Spellings of sex seen in user samplesheets, mapped to the XX/XY encoding sarek
# requires. Anything else (including a blank cell) becomes NA.
_SEX_ALIASES = {
    'f': 'XX',
    'female': 'XX',
    'xx': 'XX',
    'm': 'XY',
    'male': 'XY',
    'xy': 'XY',
    'na': 'NA',
    'unknown': 'NA',
}


def normalize_sex(value: str) -> str:
    """Map a samplesheet sex value to sarek's XX/XY/NA encoding."""
    return _SEX_ALIASES.get(value.strip().lower(), 'NA')


# Spellings of status mapped to sarek's encoding, 0 (normal) or 1 (tumor). A blank cell
# is a normal sample, as in sarek. 0.0 and 1.0 come from sheets that stored status as a
# number with blanks.
_STATUS_ALIASES = {
    '': '0',
    '0': '0',
    '0.0': '0',
    'normal': '0',
    '1': '1',
    '1.0': '1',
    'tumor': '1',
    'tumour': '1',
}


def normalize_status(value: str) -> str:
    """Map a samplesheet status value to sarek's 0 (normal) / 1 (tumor) encoding."""
    try:
        return _STATUS_ALIASES[value.strip().lower()]
    except KeyError:
        raise ValueError(
            f"Status must be 0 (normal), 1 (tumor) or blank (normal), not {value!r}"
        ) from None


def read_samplesheet(path) -> pd.DataFrame:
    """Read the samplesheet filled out in the form, keeping every value as text."""

    # The form field is marked required, but an empty table can still be submitted,
    # in which case Cirro either omits the param or writes an empty file
    if not path:
        raise ValueError("No samplesheet was provided: add at least one row of FASTQ files")
    try:
        samplesheet = pd.read_csv(path, dtype=str, keep_default_na=False)
    except pd.errors.EmptyDataError:
        samplesheet = pd.DataFrame()
    if samplesheet.empty:
        raise ValueError("The samplesheet is empty: add at least one row of FASTQ files")
    return samplesheet


def number_blank_lanes(samplesheet: pd.DataFrame) -> pd.Series:
    """Each row's lane, with blanks numbered 1, 2, ... past the lanes its sample already uses."""

    lanes = samplesheet["lane"].copy()
    for _, rows in samplesheet.groupby("sample", sort=False):
        used = set(rows["lane"])
        number = 0
        for index in rows.index[rows["lane"] == ""]:
            number += 1
            while str(number) in used:
                number += 1
            lanes[index] = str(number)
    return lanes


def fill_samplesheet(samplesheet: pd.DataFrame) -> pd.DataFrame:
    """Fill in the values sarek needs that the form lets the user leave blank."""

    return samplesheet.assign(
        patient=samplesheet["patient"].where(samplesheet["patient"] != "", samplesheet["sample"]),
        sex=samplesheet["sex"].map(normalize_sex),
        status=samplesheet["status"].map(normalize_status),
        lane=number_blank_lanes(samplesheet),
    )


def check_one_patient_per_sample(samplesheet: pd.DataFrame):
    """Fail when the rows of one sample name different patients, which sarek rejects."""

    patients = samplesheet.groupby("sample", sort=False)["patient"].unique()
    conflicts = patients[patients.map(len) > 1]
    if not conflicts.empty:
        listing = "; ".join(f"{sample} ({', '.join(names)})" for sample, names in conflicts.items())
        raise ValueError(
            f"Each sample must belong to a single patient, but these have several: {listing}. "
            "A blank patient is filled with the sample name, so give every row of a sample "
            "the same patient, or leave it blank on all of them."
        )


def check_unique_lanes(samplesheet: pd.DataFrame):
    """Fail when two rows share a patient, sample and lane, which sarek rejects."""

    key = ["patient", "sample", "lane"]
    duplicated = samplesheet.loc[samplesheet.duplicated(key, keep=False), key].drop_duplicates()
    if not duplicated.empty:
        listing = "; ".join(
            f"patient {patient}, sample {sample}, lane {lane}"
            for patient, sample, lane in duplicated.itertuples(index=False)
        )
        raise ValueError(
            f"Each row of a sample needs its own lane, but these are on more than one row: {listing}. "
            "Correct the lanes, or leave them blank to have them numbered."
        )


def check_one_status_per_sample(samplesheet: pd.DataFrame):
    """Fail when the rows of one sample disagree on whether it is a normal or a tumor."""

    statuses = samplesheet.groupby("sample", sort=False)["status"].nunique()
    mixed = statuses.index[statuses > 1].tolist()
    if mixed:
        raise ValueError(
            f"Each sample must be either normal or tumor, but these have rows of both: {', '.join(mixed)}. "
            "A blank status counts as normal (0), so set the status on every row of a tumor sample."
        )


def check_one_normal_per_tumor_patient(samplesheet: pd.DataFrame):
    """Fail when a patient with a tumor has several normal samples, which sarek rejects."""

    samples = samplesheet.drop_duplicates(["patient", "sample"])
    tumor_patients = samples.loc[samples["status"] == "1", "patient"]
    normals = samples.loc[(samples["status"] == "0") & samples["patient"].isin(tumor_patients)]
    normals = normals.groupby("patient", sort=False)["sample"].unique()
    conflicts = normals[normals.map(len) > 1]
    if not conflicts.empty:
        listing = "; ".join(f"{patient} ({', '.join(names)})" for patient, names in conflicts.items())
        raise ValueError(
            f"A tumor is paired with the one normal sample of its patient, but these patients have several: {listing}. "
            "Give each extra normal sample a patient of its own."
        )


# Tools sarek stops on when no sample has the status they need, or when any sample's
# sex is NA (subworkflows/local/samplesheet_to_channel)
_TUMOR_REQUIRING_TOOLS = frozenset({"ascat", "controlfreec", "msisensorpro", "mutect2"})
_NORMAL_REQUIRING_TOOLS = frozenset({"ascat", "deepvariant", "haplotypecaller", "msisensorpro"})
_SEX_REQUIRING_TOOLS = frozenset({"ascat", "controlfreec"})


def check_status_requirements(samplesheet: pd.DataFrame, tools: list[str]):
    """Fail when a selected tool needs a tumor or a normal sample and the samplesheet has none."""

    statuses = set(samplesheet["status"])
    for status, label, requiring in (
        ("1", "tumor", _TUMOR_REQUIRING_TOOLS),
        ("0", "normal", _NORMAL_REQUIRING_TOOLS),
    ):
        needing = sorted(requiring.intersection(tools))
        if needing and status not in statuses:
            raise ValueError(
                f"{', '.join(needing)} needs at least one {label} sample (status {status}), "
                "but the samplesheet has none."
            )


def check_sex_requirements(samplesheet: pd.DataFrame, tools: list[str]):
    """Fail when a selected tool needs the sex of every sample and some have none."""

    needing = sorted(_SEX_REQUIRING_TOOLS.intersection(tools))
    unknown = samplesheet.loc[samplesheet["sex"] == "NA", "sample"].unique().tolist()
    if needing and unknown:
        raise ValueError(
            f"{', '.join(needing)} needs the sex (XX or XY) of every sample, but these have none: "
            f"{', '.join(unknown)}. Set the Sex column on their rows."
        )


def delete_params(ds: PreprocessDataset, param_list: list[str]):
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


def warn_missing_caller_resources(ds: PreprocessDataset, tools: list[str]):
    """Warn about resources a selected caller needs that the form does not supply."""

    # ASCAT needs allele/loci reference files that Cirro does not supply by default.
    if "ascat" in tools and not (ds.params.get("ascat_alleles") and ds.params.get("ascat_loci")):
        ds.logger.warning(
            "ASCAT was selected but ascat_alleles/ascat_loci were not provided. ASCAT will "
            "fail unless both are supplied via the Extra Parameters (JSON) field."
        )

    # Mutect2 without a germline resource runs without germline-based filtering.
    if "mutect2" in tools and not ds.params.get("germline_resource"):
        ds.logger.warning(
            "Mutect2 was selected without a germline resource; no germline-based filtering "
            "of somatic calls will be performed."
        )


def warn_custom_genome_limitations(ds: PreprocessDataset, is_custom_genome: bool):
    """Warn that base recalibration cannot run against a custom BWA genome dataset.

    Custom genome datasets provide only the FASTA + BWA index — no GATK known-sites
    (dbsnp/known_indels) — so base recalibration is skipped automatically (see
    skip_baserecalibration_without_known_sites).
    """
    if not is_custom_genome:
        return
    ds.logger.warning(
        "Custom genome selected: GATK known-sites (dbsnp/known_indels) are not available, "
        "so base recalibration cannot run (it is skipped automatically — see "
        "skip_baserecalibration_without_known_sites)."
    )


def drop_annotation_for_custom_genome(ds: PreprocessDataset, is_custom_genome: bool):
    """Drop the annotation tool selection when a custom genome is used.

    VEP and snpEff need assembly-specific caches keyed off the iGenomes genome
    (vep_genome, vep_species, snpeff_db); none of those resolve for a custom genome,
    so sarek would fail at the annotation step. Dropping the selection lets the rest
    of the run complete unannotated. Must be called before annotation_tool is merged
    into the tools string.
    """
    if not is_custom_genome:
        return
    annotation_tool = ds.params.get("annotation_tool")
    if not annotation_tool:
        return
    ds.logger.warning(
        "Custom genome selected: VEP/snpEff reference data is not available, so the "
        f"selected annotation tool(s) ({', '.join(map(str, annotation_tool))}) will be skipped."
    )
    ds.remove_param("annotation_tool", force=True)


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
    The BWA index pipeline publishes ``genome.fasta`` and the flat index files
    (``genome.{amb,ann,bwt,pac,sa}``) directly into the dataset's data directory,
    so the directory itself serves as the ``--bwa`` argument (nf-core's bwa/mem
    module derives the index prefix from the ``.amb`` file).

    Dropping ``--genome`` is not sufficient on its own: sarek's nextflow.config
    defaults ``genome`` to 'GATK.GRCh38', so every reference param Cirro leaves unset
    (dict, dbsnp, known_indels, intervals, germline_resource, pon, snpeff_db, vep_*)
    would still resolve to GRCh38 iGenomes values and clash with the custom FASTA.
    ``--igenomes_ignore`` empties ``params.genomes``, so every getGenomeAttribute
    lookup returns null and the missing references are derived from the custom FASTA
    instead.
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
        ds.add_param("fasta", f"{bwa_index}/genome.fasta", overwrite=True)
        ds.add_param("fasta_fai", f"{bwa_index}/genome.fasta.fai", overwrite=True)
        ds.add_param("bwa", bwa_index, overwrite=True)
        ds.add_param("igenomes_ignore", True, overwrite=True)
        ds.remove_param("genome", force=True)
        ds.remove_param("igenomes_base", force=True)
    else:
        ds.logger.info(f"genome_source=igenomes: genome={ds.params.get('genome')!r}")


# VCF params handed to Mutect2, each paired with the index param that must accompany
# it. Both entries can resolve to files with the same name — the Cirro reference
# library only offers germline_resource.vcf.gz as a VCF, so selecting it for the panel
# of normals as well makes --pon and --germline_resource collide.
_VCF_PARAM_PAIRS = (
    ("germline_resource", "germline_resource_tbi"),
    ("pon", "pon_tbi"),
)


def colliding_vcf_copies(params: dict) -> list[tuple[str, str, str]]:
    """The (param, source URI, copy URI) of each VCF param to copy under a name of its own.

    Nextflow stages every input of a process into a single work directory, so two
    params pointing at files with the same name abort the run:

        Process ...:MUTECT2_PAIRED input file name collision -- There are multiple
        input files for each of the following file names: germline_resource.vcf.gz

    Each offending VCF is copied into the dataset's config/ folder, beside the
    samplesheet, under a name prefixed with its param key. Its index is copied as
    ``<copy>.tbi`` because GATK requires it to sit next to the VCF under a matching
    name; where no index param is set, sarek indexes the copy itself.
    """
    paths = {
        vcf_param: params[vcf_param]
        for vcf_param, _ in _VCF_PARAM_PAIRS
        if params.get(vcf_param)
    }

    by_name = {}
    for vcf_param, path in paths.items():
        by_name.setdefault(path.rsplit("/", 1)[-1], []).append(vcf_param)

    colliding = {
        vcf_param
        for vcf_params in by_name.values() if len(vcf_params) > 1
        for vcf_param in vcf_params
    }
    if not colliding:
        return []

    config_dir = params["input"].rsplit("/", 1)[0]
    copies = []
    for vcf_param, tbi_param in _VCF_PARAM_PAIRS:
        if vcf_param not in colliding:
            continue
        copy_uri = f"{config_dir}/{vcf_param}_{paths[vcf_param].rsplit('/', 1)[-1]}"
        copies.append((vcf_param, paths[vcf_param], copy_uri))
        if params.get(tbi_param):
            copies.append((tbi_param, params[tbi_param], f"{copy_uri}.tbi"))
    return copies


def require_analysis_type_binding(ds: PreprocessDataset):
    """Fail when the launch payload did not bind to the form's analysis_type block.

    process-form.json declares ``wes`` required within ``analysis_type`` and gives it a
    default, so any launch whose paramJson matches the form supplies it. Its absence means
    the whole block resolved to nothing — and ``intervals`` lives in that same block, so it
    was dropped too, which would silently turn a targeted run into a genome-wide one.

    Nothing else catches this: the form's top-level ``required`` names only the
    samplesheet, and draft-07 permits additional properties, so a flat parameter dict
    validates cleanly and then binds none of the nested JSONPaths in process-input.json.
    """
    if "wes" in ds.params:
        return

    raise ValueError(
        "Launch parameters did not match the process form: 'wes' is absent, so the "
        "analysis_type block bound nothing and 'intervals' was dropped with it. "
        "paramJson must nest parameters exactly as process-form.json declares them "
        "(analysis_type.wes, analysis_type.intervals, analysis_type.tools); "
        "a flat parameter dict binds nothing."
    )


_DEFAULT_WORKFLOW_VERSION = "3.8.1"

# Params that the extra JSON input field must never override.
# Includes both Cirro framework params (input/outdir) and Cirro-computed values
# (compute_multiplier) that must not be replaced by user input.
_PROTECTED_PARAMS = frozenset({
    "input",              # samplesheet filled out in the form
    "outdir",             # output directory — assigned by Cirro
    "igenomes_base",      # iGenomes S3 base URL
    "igenomes_ignore",    # set by resolve_reference_genome for custom genomes
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

    import boto3
    from cirro.helpers.preprocess_dataset import PreprocessDataset
    from cirro.models.s3_path import S3Path

    # Load the information for this dataset
    ds = PreprocessDataset.from_running()

    ds.logger.info(f"Starting sarek preprocess — workflow_version={ds.params.get('workflow_version', _DEFAULT_WORKFLOW_VERSION)!r}")
    ds.logger.info(f"analysis_type={ds.params.get('analysis_type')!r}, genome={ds.params.get('genome')!r}")

    require_analysis_type_binding(ds)

    tools = ds.params.get('tools')
    # Guard against None and empty list — both mean no variant caller was selected
    assert tools, "ERROR: You must select at least one variant calling tool."

    samplesheet = read_samplesheet(ds.params.get("input"))
    filled = fill_samplesheet(samplesheet)

    for line in filled.to_csv(index=None).split("\n"):
        ds.logger.info(line)

    check_one_patient_per_sample(filled)
    check_one_status_per_sample(filled)
    check_unique_lanes(filled)
    check_one_normal_per_tumor_patient(filled)
    check_status_requirements(filled, tools)
    check_sex_requirements(filled, tools)

    if not filled.equals(samplesheet):
        # Overwrite the file Cirro wrote to the dataset's config/ folder (mapped in process-input.json)
        filled.to_csv(ds.params["input"], index=None)
        ds.logger.info(f"Wrote {filled.shape[0]} row(s) to {ds.params['input']}")

    warn_missing_caller_resources(ds, tools)

    # Capture the genome source before resolve_reference_genome removes it.
    is_custom_genome = ds.params.get("genome_source") == "dataset"

    warn_custom_genome_limitations(ds, is_custom_genome)
    drop_annotation_for_custom_genome(ds, is_custom_genome)

    # Resolve the reference genome (iGenomes vs Custom BWA index) before any
    # downstream logic reads the genome param.
    resolve_reference_genome(ds)

    params = ds.params

    # Annotation tool is allowed to be empty, init empty list if it is
    # not 100% sure of annotation tool behavior if neither are selected i.e empty
    annotation_tool = params.get('annotation_tool') or []

    # Combine the two
    combined_tools = ','.join(map(str, tools + annotation_tool))
    ds.logger.info(f"Combined tools string: {combined_tools!r}")

    ds.add_param('tools', combined_tools, overwrite=True)

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
            dbnsfp = f"{ds.references_base}/VEP/{database[genome][0]}/dbNSFP4.2a_{database[genome][0].lower()}.gz"
            dbnsfp_tbi = f"{ds.references_base}/VEP/{database[genome][0]}/dbNSFP4.2a_{database[genome][0].lower()}.gz.tbi"
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
            spliceai_snv = f"{ds.references_base}/VEP/{database[genome][0]}/spliceai_scores.raw.snv.{database[genome][1]}.vcf.gz"
            spliceai_snv_tbi = f"{ds.references_base}/VEP/{database[genome][0]}/spliceai_scores.raw.snv.{database[genome][1]}.vcf.gz.tbi"
            spliceai_indel = f"{ds.references_base}/VEP/{database[genome][0]}/spliceai_scores.raw.indel.{database[genome][1]}.vcf.gz"
            spliceai_indel_tbi = f"{ds.references_base}/VEP/{database[genome][0]}/spliceai_scores.raw.indel.{database[genome][1]}.vcf.gz.tbi"
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
            pon = f"{ds.references_base}/igenomes/Homo_sapiens/GATK/GRCh37/Annotation/GATKBundle/Mutect2-WGS-panel-b37.vcf.gz"
            pon_tbi = f"{ds.references_base}/igenomes/Homo_sapiens/GATK/GRCh37/Annotation/GATKBundle/Mutect2-WGS-panel-b37.vcf.gz.tbi"
            ds.add_param('pon', pon, overwrite=True)
            ds.add_param('pon_tbi', pon_tbi, overwrite=True)
            ds.logger.info("PON: added GRCh37 somatic panel of normals")
        elif genome == "GATK.GRCh38":
            pon = f"{ds.references_base}/igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000g_pon.hg38.vcf.gz"
            pon_tbi = f"{ds.references_base}/igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000g_pon.hg38.vcf.gz.tbi"
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

    # Give the Mutect2 VCF params distinct file names, now that every source of a
    # germline_resource/pon value (form, iGenomes defaults, extra JSON) has been applied.
    copies = colliding_vcf_copies(ds.params)
    if not copies:
        ds.logger.info("VCF inputs: no file name collisions to resolve")
    for param, source_uri, copy_uri in copies:
        source, dest = S3Path(source_uri), S3Path(copy_uri)
        assert source.valid, f"Cannot stage a copy of --{param}: {source_uri} is not an S3 path"
        ds.logger.info(f"VCF inputs: staging {source_uri} as {copy_uri}")
        boto3.client("s3").copy({"Bucket": source.bucket, "Key": source.key}, dest.bucket, dest.key)
        ds.add_param(param, copy_uri, overwrite=True)

    # Remove any parameters not defined in the nf-core/sarek nextflow_schema.json for
    # the selected workflow version. This also cleans up Cirro-only housekeeping params
    # (workflow_version, extra_params_json) that must not reach Nextflow.
    filter_params_by_schema(ds)

    # With all params populated, skip base recalibration for custom genomes that
    # lack known-sites resources (dbsnp/known_indels). iGenomes supplies these.
    skip_baserecalibration_without_known_sites(ds, is_custom_genome)

    # log all params
    ds.logger.info(f"Final params ({len(ds.params)} total): {ds.params}")
