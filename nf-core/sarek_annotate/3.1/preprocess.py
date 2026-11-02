#!/usr/bin/env python3

from __future__ import annotations

import json
import re
import urllib.error
import urllib.request
from pathlib import PurePosixPath
from typing import TYPE_CHECKING

import pandas as pd

# The SDK is imported in __main__, so test_preprocess.py runs without it
if TYPE_CHECKING:
    from cirro.helpers.preprocess_dataset import PreprocessDataset


def read_samplesheet(path) -> pd.DataFrame:
    """Read the samplesheet filled out in the form, keeping every value as text."""

    # The form field is marked required, but an empty table can still be submitted,
    # in which case Cirro either omits the param or writes an empty file
    if not path:
        raise ValueError("No samplesheet was provided: add at least one row with a VCF file")
    try:
        samplesheet = pd.read_csv(path, dtype=str, keep_default_na=False)
    except pd.errors.EmptyDataError:
        samplesheet = pd.DataFrame()
    if samplesheet.empty:
        raise ValueError("The samplesheet is empty: add at least one row with a VCF file")
    return samplesheet


def vcf_stem(vcf: str) -> str:
    """The VCF's file name without its .vcf or .vcf.gz extension."""
    return re.sub(r"\.vcf(\.gz)?$", "", PurePosixPath(vcf).name)


def distinguish_repeated_samples(samplesheet: pd.DataFrame) -> pd.DataFrame:
    """Give every row a patient and sample of its own.

    sarek stops when two rows share a patient and sample, even when their variant callers
    differ (nf-core/sarek#2201, in every release up to 3.11.0). sarek's variantcalled.csv
    lists one row per VCF, so a sample called by several callers repeats, and so does one
    whose caller writes separate SNV and indel files. Repeated rows get their variant
    caller added to the sample name; rows that still repeat are named after their VCF.
    """
    samplesheet = samplesheet.copy()

    def repeated() -> pd.Series:
        return samplesheet.duplicated(["patient", "sample"], keep=False)

    with_caller = repeated() & (samplesheet["variantcaller"] != "")
    samplesheet.loc[with_caller, "sample"] = (
        samplesheet.loc[with_caller, "sample"] + "_" + samplesheet.loc[with_caller, "variantcaller"]
    )
    still_repeated = repeated()
    samplesheet.loc[still_repeated, "sample"] = samplesheet.loc[still_repeated, "vcf"].map(vcf_stem)

    duplicated = samplesheet.loc[repeated(), ["patient", "sample"]].drop_duplicates()
    if not duplicated.empty:
        listing = ", ".join(f"{patient}/{sample}" for patient, sample in duplicated.itertuples(index=False))
        raise ValueError(
            "Each VCF needs its own patient and sample, but these patient/sample pairs are "
            f"on more than one row even after adding the variant caller and VCF file name: {listing}"
        )
    return samplesheet


def set_annotation_params(ds: PreprocessDataset):

    # The annotate step requires at least one tool (sarek errors otherwise).
    annotation_tool = ds.params.get('annotation_tool') or []
    if not annotation_tool:
        ds.logger.error(
            "No annotation tool selected. The annotate step requires at least one of "
            "VEP or snpEff."
        )

    # Combine the two
    tools = ','.join(map(str, annotation_tool))

    ds.add_param('tools', tools, overwrite=True)
    genome = ds.params.get('genome')
    ds.logger.info(f"Reference genome: {genome!r}")

    # annotation_tool is a Cirro UI staging param — always remove it
    ds.remove_param('annotation_tool', force=True)

    # construct logic for dbNSFP & SpliceAI
    dbnsfp_param = ds.params.get('vep_dbnsfp') # true or null
    spliceai = ds.params.get('vep_spliceai') # true or null
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


_DEFAULT_WORKFLOW_VERSION = "3.8.1"

# Params set by Cirro infrastructure that must not be overridden by user-supplied extra JSON.
_PROTECTED_PARAMS = frozenset({
    "input",
    "outdir",
    "igenomes_base",
    "vep_cache",
    "snpeff_cache",
    "monochrome_logs",
    "step",  # hardcoded to "annotate" by process-input.json
})

# igenomes_base, vep_cache, snpeff_cache, and monochrome_logs are Cirro-pinned S3 paths
# that must survive schema filtering regardless of sarek version.
_CIRRO_PASSTHROUGH_PARAMS = frozenset({
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

    # Get the currently running dataset
    ds = PreprocessDataset.from_running()

    ds.logger.info(f"Starting sarek_annotate preprocess — workflow_version={ds.params.get('workflow_version', _DEFAULT_WORKFLOW_VERSION)!r}")

    samplesheet = read_samplesheet(ds.params.get("input"))
    distinct = distinguish_repeated_samples(samplesheet)
    if not distinct.equals(samplesheet):
        # Overwrite the file Cirro wrote to the dataset's config/ folder (mapped in process-input.json)
        distinct.to_csv(ds.params["input"], index=None)
    ds.logger.info("Samplesheet:")
    ds.logger.info(distinct.to_csv(index=None))

    set_annotation_params(ds)

    apply_extra_json_params(ds)
    filter_params_by_schema(ds)

    ds.logger.info(f"Final params ({len(ds.params)} total): {ds.params}")
