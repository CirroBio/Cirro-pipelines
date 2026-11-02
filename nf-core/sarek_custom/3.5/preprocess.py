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

    from cirro.helpers.preprocess_dataset import PreprocessDataset

    ds = PreprocessDataset.from_running()

    ds.logger.info(f"Starting sarek_custom preprocess — workflow_version={ds.params.get('workflow_version', _DEFAULT_WORKFLOW_VERSION)!r}")

    samplesheet = read_samplesheet(ds.params.get("input"))
    filled = fill_samplesheet(samplesheet)

    for line in filled.to_csv(index=None).split("\n"):
        ds.logger.info(line)

    check_one_patient_per_sample(filled)
    check_unique_lanes(filled)

    if not filled.equals(samplesheet):
        # Overwrite the file Cirro wrote to the dataset's config/ folder (mapped in process-input.json)
        filled.to_csv(ds.params["input"], index=None)
        ds.logger.info(f"Wrote {filled.shape[0]} row(s) to {ds.params['input']}")

    tools = ds.params.get("tools")
    assert tools, "ERROR: You must select at least one variant calling tool."

    # Convert the tools parameter from a list to a comma delimited string
    ds.add_param('tools', ",".join(ds.params.get("tools")), overwrite=True)
    ds.logger.info(f"tools: {ds.params.get('tools')!r}")

    # If an intervals file was not selected, use --no_intervals
    if not ds.params.get("intervals"):
        ds.add_param("no_intervals", True)
        ds.logger.info("No intervals file selected — adding --no_intervals flag")

    # If the user selected to save alignments in BAM format, set the flag
    if ds.params.get("alignment_format", "CRAM") == "BAM":
        ds.add_param("save_output_as_bam", True)
    ds.remove_param("alignment_format", force=True)

    apply_extra_json_params(ds)
    filter_params_by_schema(ds)

    ds.logger.info(f"Final params ({len(ds.params)} total): {ds.params}")
