#!/usr/bin/env python3

import pandas as pd


def read_table(path, name: str, rows: str) -> pd.DataFrame:
    """Read a table filled out in the form, keeping every value as text."""

    # The form field is marked required, but an empty table can still be submitted,
    # in which case Cirro either omits the param or writes an empty file
    if not path:
        raise ValueError(f"No {name} was provided: add {rows}")
    try:
        table = pd.read_csv(path, dtype=str, keep_default_na=False)
    except pd.errors.EmptyDataError:
        table = pd.DataFrame()
    if table.empty:
        raise ValueError(f"The {name} is empty: add {rows}")
    return table


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


def manual_probe_barcodes(params: dict) -> pd.DataFrame:
    """The sample entered for each probe barcode field, bc001 to bc016.

    The form groups these under samples, but process-input.json maps each to a flat
    bcNNN key. A barcode the user left blank names no sample.
    """
    return pd.DataFrame(
        [
            dict(sample_id=params[name], probe_barcode_ids=name.upper())
            for name in sorted(params)
            if name.startswith("bc0") and params[name]
        ],
        columns=["sample_id", "probe_barcode_ids"],
    )


def select_probe_barcodes(table: pd.DataFrame) -> pd.DataFrame:
    """The sample_id and probe_barcode_ids columns of a barcode table, without the rows
    that leave either one blank."""

    # CellRanger 10.1 also names this column sample_barcode_ids; 10.0 accepts only
    # probe_barcode_ids, and so does the workflow's config builder
    table = table.rename(columns={"sample_barcode_ids": "probe_barcode_ids"})
    for cname in ["sample_id", "probe_barcode_ids"]:
        if cname not in table.columns:
            raise ValueError(
                f"Barcode table must have a column named {cname} (found: {', '.join(table.columns)})"
            )
    table = table[["sample_id", "probe_barcode_ids"]].fillna("")
    return table[(table != "").all(axis=1)]


def join_probe_barcodes(probe_barcodes: pd.DataFrame) -> pd.DataFrame:
    """One row per sample, joining the barcodes of a sample listed on several rows with
    a pipe (e.g. BC001|BC002), as CellRanger takes them."""

    if probe_barcodes.empty:
        raise ValueError("No sample was assigned a probe barcode")
    return probe_barcodes.groupby("sample_id")["probe_barcode_ids"].agg("|".join).reset_index()


if __name__ == "__main__":

    # Imported here so that test_preprocess.py runs without the SDK
    from cirro.helpers.preprocess_dataset import PreprocessDataset

    ds = PreprocessDataset.from_running()

    # Datasets launched before the form offered a choice carry no flex_version
    flex_version = ds.params.get("flex_version") or "v1"
    ds.logger.info(f"Flex chemistry: {flex_version}")

    # Flex v2 chemistries were added in CellRanger 10.0
    if flex_version == "v2":
        msg = f"Flex v2 requires CellRanger 10.0 or later (selected: {ds.params['cellranger_version']})"
        assert int(ds.params["cellranger_version"].split(".")[0]) >= 10, msg

    for param in unset_umi_thresholds(ds.params):
        ds.remove_param(param, force=True)

    # The GEX reference and the 10X probe sets built against it, by Flex version.
    # CellRanger refuses a probe set whose reference version differs from the GEX
    # reference, and the v2 probe sets exist for the 2024-A references only.
    transcriptome, probe_sets = {
        "Homo sapiens (GRCh38-2024)": (
            "refdata-gex-GRCh38-2024-A",
            {
                "v1": "Chromium_Human_Transcriptome_Probe_Set_v1.1.0_GRCh38-2024-A.csv",
                "v2": "Chromium_Human_Transcriptome_Probe_Set_v2.0.0_GRCh38-2024-A.csv"
            }
        ),
        "Homo sapiens (GRCh38-2020)": (
            "refdata-gex-GRCh38-2020-A",
            {"v1": "Chromium_Human_Transcriptome_Probe_Set_v1.0.1_GRCh38-2020-A.csv"}
        ),
        "Mus musculus (GRCm39-2024)": (
            "refdata-gex-GRCm39-2024-A",
            {
                "v1": "Chromium_Mouse_Transcriptome_Probe_Set_v1.1.0_GRCm39-2024-A.csv",
                "v2": "Chromium_Mouse_Transcriptome_Probe_Set_v2.0.0_GRCm39-2024-A.csv"
            }
        ),
        "Mus musculus (mm10-2020)": (
            "refdata-gex-mm10-2020-A",
            {"v1": "Chromium_Mouse_Transcriptome_Probe_Set_v1.0.1_mm10-2020-A.csv"}
        )
    }[ds.params["reference"]]
    msg = f"No Flex {flex_version} probe set is available for {ds.params['reference']}"
    assert flex_version in probe_sets, msg
    default_probe_set = probe_sets[flex_version]

    ds.add_param(
        "transcriptome_dir",
        # a reference directory, which HealthOmics reads as an object without it
        f"{ds.references_base}/cellranger/{transcriptome}/"
    )

    # If the user did not select a custom probe set
    if ds.params.get("probe_set") is None or ds.params.get("probe_set") == "":
        ds.add_param(
            "probe_set",
            f"{ds.references_base}/cellranger/flex/{default_probe_set}",
            overwrite=True
        )

    # Datasets launched before the form offered a choice carry no barcode_source
    barcode_source = ds.params.get("barcode_source") or "manual"
    ds.logger.info(f"Probe barcode source: {barcode_source}")

    if barcode_source == "manual":
        probe_barcodes = manual_probe_barcodes(ds.params)
    elif barcode_source == "table":
        probe_barcodes = select_probe_barcodes(read_table(
            ds.params.get("barcode_table"), "probe barcode table", "one row per sample"
        ))
    else:
        # A Spreadsheet (CSV) reference
        probe_barcodes = select_probe_barcodes(pd.read_csv(
            ds.params["barcode_reference"], dtype=str, keep_default_na=False, skipinitialspace=True
        ))
    probe_barcodes = join_probe_barcodes(probe_barcodes)

    # Save the sample barcode spreadsheet to the dataset's config/ folder
    # (mapped in process-input.json)
    ds.logger.info("Sample probe barcodes specified:")
    ds.logger.info(probe_barcodes.to_csv(index=None))
    probe_barcodes.to_csv(ds.params["probe_barcodes"], index=None)

    samplesheet = read_table(ds.params.get("samples"), "samplesheet", "one row per library")
    libraries = drop_repeated_rows(samplesheet)
    ds.logger.info("Samplesheet:")
    ds.logger.info(libraries.to_csv(index=None))
    if len(libraries) < len(samplesheet):
        # Overwrite the file Cirro wrote to the dataset's config/ folder (mapped in process-input.json)
        libraries.to_csv(ds.params["samples"], index=None)

    # Log the parameters present
    for k, v in ds.params.items():
        ds.logger.info(f"{k}: {v}")

    # Consumed above; not parameters the workflow declares.
    barcode_params = [param for param in ds.params if param.startswith("bc0")]
    for param in ["flex_version", "reference", "barcode_source", "barcode_table", "barcode_reference"] + barcode_params:
        ds.remove_param(param, force=True)
