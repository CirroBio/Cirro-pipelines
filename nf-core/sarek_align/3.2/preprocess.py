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


def warn_custom_genome_limitations(ds: PreprocessDataset):
    """Warn when a custom BWA genome is used that base recalibration cannot run.

    Custom genome datasets provide only the FASTA + BWA index — no GATK known-sites
    (dbsnp/known_indels) — so base recalibration will fail. Must be called before
    resolve_reference_genome removes ``genome_source``.
    """
    if ds.params.get("genome_source") != "dataset":
        return
    ds.logger.warning(
        "Custom genome selected: GATK known-sites (dbsnp/known_indels) are not available, "
        "so base recalibration cannot run (it is skipped automatically — see "
        "skip_baserecalibration_without_known_sites)."
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

    For iGenomes the curated ``genome`` key is passed through unchanged, and
    ``aligner`` (bwa-mem/bwa-mem2/parabricks) is a real sarek param that is
    defaulted and left in place for both branches -- unlike ``genome_source``/
    ``bwa_index``/``bwamem2_index``, which are Cirro-form-only fields that must
    never reach Nextflow. For a custom genome the user selects a pre-built BWA
    or BWA-MEM2 index dataset (mutually exclusive in the form, keyed off
    ``aligner``); we point ``--fasta``/``--bwa``-or-``--bwamem2`` at that
    dataset and drop
    ``--genome``/``--igenomes_base``. Both index pipelines publish
    ``genome.fasta`` and their respective flat index files directly into the
    dataset's data directory, so the directory itself serves as the index
    argument (nf-core's bwa/mem and bwamem2/mem modules derive the index
    prefix from the index files present).

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

    # A real sarek param (unlike genome_source/bwa_index/bwamem2_index below), so it
    # is defaulted and written back once here instead of being removed and re-added
    # separately in each branch.
    aligner = ds.params.get("aligner") or "bwa-mem"
    ds.add_param("aligner", aligner, overwrite=True)

    bwa_index = ds.params.get("bwa_index")
    ds.remove_param("bwa_index", force=True)

    bwamem2_index = ds.params.get("bwamem2_index")
    ds.remove_param("bwamem2_index", force=True)

    if genome_source == "dataset":
        use_bwamem2 = aligner == "bwa-mem2"
        custom_index = bwamem2_index if use_bwamem2 else bwa_index
        if not custom_index:
            raise ValueError(
                f"Custom Genome selected with aligner={aligner!r} but no matching "
                "genome index dataset was provided."
            )
        ds.logger.info(f"genome_source=dataset: using custom {aligner} index at {custom_index}")
        ds.add_param("fasta", f"{custom_index}/genome.fasta", overwrite=True)
        ds.add_param("fasta_fai", f"{custom_index}/genome.fasta.fai", overwrite=True)
        if use_bwamem2:
            ds.add_param("bwamem2", custom_index, overwrite=True)
        else:
            ds.add_param("bwa", custom_index, overwrite=True)
        ds.add_param("igenomes_ignore", True, overwrite=True)
        ds.remove_param("genome", force=True)
        ds.remove_param("igenomes_base", force=True)
    else:
        ds.logger.info(
            f"genome_source=igenomes: genome={ds.params.get('genome')!r}, aligner={aligner!r}"
        )


_VCF_PARAM_PAIRS = (
    ("dbsnp", "dbsnp_tbi"),
    ("known_indels", "known_indels_tbi"),
)


def colliding_vcf_copies(params: dict) -> list[tuple[str, str, str]]:
    """The (param, source URI, copy URI) of each known-sites VCF to copy under a name of its own.

    Nextflow stages every input of a process into a single work directory, so dbsnp and
    known_indels picked from the same germline_resource.vcf.gz reference type collide.
    Each offending VCF is copied into the dataset's config/ folder, beside the
    samplesheet, under a name prefixed with its param key, and its index as
    ``<copy>.tbi``.
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
        "(analysis_type.wes, analysis_type.intervals, "
        "analysis_type.genome_selection.igenomes.genome, advanced_options.*, "
        "read_trimming_options.*); a flat parameter dict binds nothing."
    )


_DEFAULT_WORKFLOW_VERSION = "3.8.1"

# Params set by Cirro infrastructure or computed by this script that must not
# be overridden by user-supplied extra JSON.
_PROTECTED_PARAMS = frozenset({
    "input",
    "outdir",
    "igenomes_base",
    "igenomes_ignore",  # set by resolve_reference_genome for custom genomes
    "vep_cache",
    "snpeff_cache",
    "monochrome_logs",
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
        key
        for key in list(ds.params.keys())
        if key not in allowed and key not in _CIRRO_PASSTHROUGH_PARAMS
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

    ds.logger.info(f"Starting sarek_align preprocess — workflow_version={ds.params.get('workflow_version', _DEFAULT_WORKFLOW_VERSION)!r}")

    require_analysis_type_binding(ds)

    samplesheet = read_samplesheet(ds.params.get("input"))
    filled = fill_samplesheet(samplesheet)

    for line in filled.to_csv(index=None).split("\n"):
        ds.logger.info(line)

    check_one_patient_per_sample(filled)
    check_one_status_per_sample(filled)
    check_unique_lanes(filled)
    check_one_normal_per_tumor_patient(filled)

    if not filled.equals(samplesheet):
        # Overwrite the file Cirro wrote to the dataset's config/ folder (mapped in process-input.json)
        filled.to_csv(ds.params["input"], index=None)
        ds.logger.info(f"Wrote {filled.shape[0]} row(s) to {ds.params['input']}")

    # Warn about custom-genome limitations while genome_source is still present.
    warn_custom_genome_limitations(ds)

    # Capture the genome source before resolve_reference_genome removes it.
    is_custom_genome = ds.params.get("genome_source") == "dataset"

    # Resolve the reference genome (iGenomes vs Custom BWA index)
    resolve_reference_genome(ds)

    # `compute_multiplier` == 2 for WGS and 1 for WES; consumed by process-compute.config
    wes = ds.params["wes"]
    compute_multiplier = int(2 - int(wes))
    ds.add_param("compute_multiplier", compute_multiplier)
    ds.logger.info(f"compute_multiplier={compute_multiplier} ({'WES' if wes else 'WGS'})")

    # If an intervals file was not selected, use --no_intervals
    if not ds.params.get("intervals"):
        ds.add_param("no_intervals", True)
        ds.logger.info("No intervals file selected — adding --no_intervals flag")

    apply_extra_json_params(ds)

    copies = colliding_vcf_copies(ds.params)
    if not copies:
        ds.logger.info("VCF inputs: no file name collisions to resolve")
    for param, source_uri, copy_uri in copies:
        source, dest = S3Path(source_uri), S3Path(copy_uri)
        assert source.valid, f"Cannot stage a copy of --{param}: {source_uri} is not an S3 path"
        ds.logger.info(f"VCF inputs: staging {source_uri} as {copy_uri}")
        boto3.client("s3").copy({"Bucket": source.bucket, "Key": source.key}, dest.bucket, dest.key)
        ds.add_param(param, copy_uri, overwrite=True)

    filter_params_by_schema(ds)

    # With all params populated, skip base recalibration for custom genomes that
    # lack known-sites resources (dbsnp/known_indels). iGenomes supplies these.
    skip_baserecalibration_without_known_sites(ds, is_custom_genome)

    ds.logger.info(f"Final params ({len(ds.params)} total): {ds.params}")
