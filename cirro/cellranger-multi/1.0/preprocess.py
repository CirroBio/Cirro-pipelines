#!/usr/bin/env python3

import pandas as pd


def reference_dir(value, references_base: str):
    """A CellRanger reference directory param, made absolute against the references bucket.

    Form values name a path under that bucket rather than a full URI, so the
    configuration is not tied to one deployment. A value that is already absolute
    is left alone, which keeps datasets created before that change re-runnable.
    HealthOmics reads an S3 value with no trailing separator as an object key and
    rejects the run, so the result always ends with one.
    """
    if not value:
        return value
    if not value.startswith("s3://"):
        value = f"{references_base}/{value}"
    return value.rstrip("/") + "/"


def read_optional_table(path) -> pd.DataFrame:
    """Read a table filled out in the form, keeping every value as text.

    Cirro omits the param when the field was never touched, and writes an empty file
    when it has no rows; both read as an empty table.
    """
    if not path:
        return pd.DataFrame()
    try:
        return pd.read_csv(path, dtype=str, keep_default_na=False)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def read_samplesheet(path) -> pd.DataFrame:
    """Read the samplesheet filled out in the form, keeping every value as text."""

    # The form field is marked required, but an empty table can still be submitted
    if not path:
        raise ValueError("No samplesheet was provided: add one row per library")
    samplesheet = read_optional_table(path)
    if samplesheet.empty:
        raise ValueError("The samplesheet is empty: add one row per library")
    return samplesheet


def drop_repeated_rows(samplesheet: pd.DataFrame) -> pd.DataFrame:
    """Rows filled from a dataset's file names repeat each library once per lane and
    read, and each row becomes a library in the CellRanger config."""
    return samplesheet.drop_duplicates(ignore_index=True)


def unset_umi_thresholds(params: dict) -> list:
    """The UMI thresholds left blank, which are dropped so the CellRanger default applies.

    CellRanger added emptydrops-minimum-umis in 7.1 and global-minimum-umis in 8.0,
    so a threshold that is set needs at least that version.
    """
    selected = params["cellranger_version"]
    version = tuple(int(n) for n in selected.split(".")[:2])
    unset = []
    for param, min_version in [("emptydrops_minimum_umis", (7, 1)), ("global_minimum_umis", (8, 0))]:
        if params.get(param) is None:
            unset.append(param)
        elif version < min_version:
            raise ValueError(
                f"{param} requires CellRanger {min_version[0]}.{min_version[1]} or later "
                f"(selected: {selected})"
            )
    return unset


if __name__ == "__main__":

    # Imported here so that test_preprocess.py runs without the SDK
    from cirro.helpers.preprocess_dataset import PreprocessDataset

    ds = PreprocessDataset.from_running()

    for param in ("transcriptome_dir", "vdj_dir"):
        resolved = reference_dir(ds.params.get(param), ds.references_base)
        if resolved != ds.params.get(param):
            ds.add_param(param, resolved, overwrite=True)

    samplesheet = read_samplesheet(ds.params.get("grouping"))
    libraries = drop_repeated_rows(samplesheet)
    if len(libraries) < len(samplesheet):
        # Overwrite the file Cirro wrote to the dataset's config/ folder (mapped in process-input.json)
        libraries.to_csv(ds.params["grouping"], index=None)

    for param in unset_umi_thresholds(ds.params):
        ds.remove_param(param, force=True)

    # Left blank, the workflow's empty feature reference is used
    if ds.params.get("feature_csv") is None:
        ds.remove_param("feature_csv", force=True)

    # Read below, only when the user indicates fixed RNA profiling
    probe_barcodes = pd.DataFrame()

    if ds.params.get("is_frp"):
        ds.logger.info("User indicated that this is fixed RNA profiling")

        # Add the probe set built against the same reference version as the
        # transcriptome -- CellRanger refuses a probe set whose version differs
        is_2024 = "2024-A" in ds.params["transcriptome_dir"]
        if "GRCh38" in ds.params["transcriptome_dir"]:
            ds.logger.info("Adding human reference probe set")
            probe_set = (
                "Chromium_Human_Transcriptome_Probe_Set_v1.1.0_GRCh38-2024-A.csv" if is_2024
                else "Chromium_Human_Transcriptome_Probe_Set_v1.0.1_GRCh38-2020-A.csv"
            )
        else:
            ds.logger.info("Adding mouse reference probe set")
            probe_set = (
                "Chromium_Mouse_Transcriptome_Probe_Set_v1.1.0_GRCm39-2024-A.csv" if is_2024
                else "Chromium_Mouse_Transcriptome_Probe_Set_v1.0.1_mm10-2020-A.csv"
            )
        ds.add_param("probes_csv", f"{ds.references_base}/cellranger/flex/{probe_set}")

        probe_barcodes = read_optional_table(ds.params.get("probe_barcodes"))
        ds.logger.info(f"Detected {probe_barcodes.shape[0]:,} samples in the probe barcode table")

    # Without rows the workflow's empty probe barcode table is used, as it is for
    # singleplex samples
    if probe_barcodes.empty:
        ds.remove_param("probe_barcodes", force=True)

    # Log the parameters present
    for k, v in ds.params.items():
        ds.logger.info(f"{k}: {v}")

    # Consumed above; not a parameter the workflow declares.
    ds.remove_param("is_frp", force=True)
