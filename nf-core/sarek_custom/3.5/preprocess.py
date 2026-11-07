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
# is a normal sample, as in sarek.
_STATUS_ALIASES = {
    '': '0',
    '0': '0',
    'normal': '0',
    '1': '1',
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


# The offered callers sarek 3.8.1 runs on a tumor: paired callers against the normal of
# its patient, tumor-only callers when its patient has no normal
# (subworkflows/local/bam_variant_calling_somatic_all and _tumor_only_all)
_PAIRED_CALLERS = frozenset({"cnvkit", "freebayes", "manta", "mutect2", "strelka", "tiddit"})
_TUMOR_ONLY_CALLERS = frozenset({"cnvkit", "freebayes", "lofreq", "manta", "mpileup", "mutect2", "tiddit"})

# Tools sarek stops on when no sample has the status they need (samplesheet_to_channel)
_TUMOR_REQUIRING_TOOLS = frozenset({"mutect2"})
_NORMAL_REQUIRING_TOOLS = frozenset({"deepvariant", "haplotypecaller"})


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


def check_tumors_are_called(samplesheet: pd.DataFrame, tools: list[str], has_intervals: bool):
    """Fail when none of the selected callers would call some of the tumor samples."""

    samples = samplesheet.drop_duplicates(["patient", "sample"])
    tumors = samples.loc[samples["status"] == "1"]
    paired = tumors["patient"].isin(samples.loc[samples["status"] == "0", "patient"])

    # CNVkit builds its tumor-only reference from the intervals, and calls nothing without them
    tumor_only_callers = _TUMOR_ONLY_CALLERS if has_intervals else _TUMOR_ONLY_CALLERS - {"cnvkit"}
    if not paired.all() and not tumor_only_callers.intersection(tools):
        raise ValueError(
            "These tumor samples have no normal sample under the same patient, and none of the "
            f"selected callers calls a tumor alone: {', '.join(tumors.loc[~paired, 'sample'])}. "
            "Give each tumor the patient of its matched normal, or also select one of: "
            f"{', '.join(sorted(tumor_only_callers))}."
        )
    if paired.any() and not _PAIRED_CALLERS.intersection(tools):
        raise ValueError(
            "None of the selected callers calls a tumor against its normal, so these tumor samples "
            f"would not be called: {', '.join(tumors.loc[paired, 'sample'])}. Also select one of: "
            f"{', '.join(sorted(_PAIRED_CALLERS))}."
        )


_MUTECT2_RESOURCE_PARAMS = ("germline_resource", "pon")


def mutect2_resources(params: dict, tools: list[str]) -> dict:
    """The germline resource and panel of normals picked for Mutect2, if it was selected.

    The Germline Resource File reference type stores the VCF's index beside it, so the
    index param points there. sarek indexes a panel of normals itself.
    """
    if "mutect2" not in tools:
        return {}
    resources = {param: params[param] for param in _MUTECT2_RESOURCE_PARAMS if params.get(param)}
    if "germline_resource" in resources:
        resources["germline_resource_tbi"] = f"{resources['germline_resource']}.tbi"
    return resources


def pon_copy_uri(params: dict) -> str | None:
    """Where to copy the panel of normals when its file name is the germline resource's.

    The references library saves every VCF as germline_resource.vcf.gz, and Nextflow
    aborts a process given two inputs of the same name, so the panel of normals is
    copied under its own name into the config/ folder that holds the samplesheet.
    """
    pon, germline_resource = params.get("pon"), params.get("germline_resource")
    if not (pon and germline_resource):
        return None
    pon_name = pon.rsplit("/", 1)[-1]
    if pon_name != germline_resource.rsplit("/", 1)[-1]:
        return None
    return f"{params['input'].rsplit('/', 1)[0]}/pon_{pon_name}"


_DEFAULT_WORKFLOW_VERSION = "3.8.1"

# Params set by Cirro infrastructure that must not be overridden by user-supplied extra JSON.
_PROTECTED_PARAMS = frozenset({
    "input",
    "outdir",
    "vep_cache",
    "snpeff_cache",
    "monochrome_logs",
    "igenomes_ignore",  # hardcoded true — custom genome workflow does not use iGenomes
    "skip_tools",       # hardcoded to "baserecalibrator" by process-input.json
})

# sarek_custom has no compute_multiplier; wes passes through as a native nf-core param
_CIRRO_PASSTHROUGH_PARAMS = frozenset({"optical_duplicate_pixel_distance"})

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

    import boto3
    from cirro.helpers.preprocess_dataset import PreprocessDataset
    from cirro.models.s3_path import S3Path

    ds = PreprocessDataset.from_running()

    ds.logger.info(f"Starting sarek_custom preprocess — workflow_version={ds.params.get('workflow_version', _DEFAULT_WORKFLOW_VERSION)!r}")

    tools = ds.params.get("tools")
    assert tools, "ERROR: You must select at least one variant calling tool."

    samplesheet = read_samplesheet(ds.params.get("input"))
    filled = fill_samplesheet(samplesheet)

    check_one_patient_per_sample(filled)
    check_one_status_per_sample(filled)
    check_unique_lanes(filled)
    check_one_normal_per_tumor_patient(filled)
    check_status_requirements(filled, tools)
    check_tumors_are_called(filled, tools, has_intervals=bool(ds.params.get("intervals")))

    if not filled.equals(samplesheet):
        # Overwrite the file Cirro wrote to the dataset's config/ folder (mapped in process-input.json)
        filled.to_csv(ds.params["input"], index=None)
        ds.logger.info(f"Wrote {filled.shape[0]} row(s) to {ds.params['input']}")

    # Convert the tools parameter from a list to a comma delimited string
    ds.add_param('tools', ",".join(tools), overwrite=True)
    ds.logger.info(f"tools: {ds.params.get('tools')!r}")

    # The Mutect2 pickers keep their values while hidden, so pass only what Mutect2 uses
    resources = mutect2_resources(ds.params, tools)
    for param in _MUTECT2_RESOURCE_PARAMS:
        ds.remove_param(param, force=True)
    for param, uri in resources.items():
        ds.add_param(param, uri, overwrite=True)

    # If an intervals file was not selected, use --no_intervals
    if not ds.params.get("intervals"):
        ds.add_param("no_intervals", True)
        ds.logger.info("No intervals file selected — adding --no_intervals flag")

    # If the user selected to save alignments in BAM format, set the flag
    if ds.params.get("alignment_format", "CRAM") == "BAM":
        ds.add_param("save_output_as_bam", True)
    ds.remove_param("alignment_format", force=True)

    apply_extra_json_params(ds)

    pon_copy = pon_copy_uri(ds.params)
    if pon_copy:
        source, dest = S3Path(ds.params["pon"]), S3Path(pon_copy)
        assert source.valid, f"Cannot copy the panel of normals: {ds.params['pon']} is not an S3 path"
        ds.logger.info(f"Copying the panel of normals to {pon_copy}, apart from the germline resource of the same name")
        boto3.client("s3").copy({"Bucket": source.bucket, "Key": source.key}, dest.bucket, dest.key)
        ds.add_param("pon", pon_copy, overwrite=True)

    filter_params_by_schema(ds)

    ds.logger.info(f"Final params ({len(ds.params)} total): {ds.params}")
