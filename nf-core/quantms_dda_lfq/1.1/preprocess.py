#!/usr/bin/env python3

import pandas as pd

FACTOR = "factor value[condition]"

# Required by the SDRF schema that quantms validates. quantms reads them only as the
# conditions to compare when there is no factor value column.
NOT_APPLICABLE_COLUMNS = [
    "characteristics[organism part]",
    "characteristics[disease]",
    "characteristics[cell type]",
    "characteristics[biological replicate]",
]


def read_samplesheet(path) -> pd.DataFrame:
    """Read the tab-separated SDRF table filled out in the form, keeping every value as text."""

    # The form field is marked required, but an empty table can still be submitted,
    # in which case Cirro either omits the param or writes an empty file
    if not path:
        raise ValueError("No SDRF was provided: add one row per data file")
    try:
        samplesheet = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    except pd.errors.EmptyDataError:
        samplesheet = pd.DataFrame()
    if samplesheet.empty:
        raise ValueError("The SDRF is empty: add one row per data file")
    return samplesheet


def make_sdrf(
    samplesheet: pd.DataFrame,
    organism: str,
    enzyme: str,
    fragment_tolerance: str,
    precursor_tolerance: str
) -> pd.DataFrame:
    """The SDRF quantms reads: the table with its blanks filled, plus the per-run columns."""

    # sdrf-pipelines rejects values with leading or trailing whitespace
    sdrf = samplesheet.apply(lambda values: values.str.strip())

    for column in ["source name", "comment[file uri]"]:
        blank = sdrf.index[sdrf[column] == ""]
        if len(blank) > 0:
            rows = ", ".join(str(i + 1) for i in blank)
            raise ValueError(f"Every row needs a {column}, but it is blank in row(s) {rows}")

    sdrf["comment[data file]"] = sdrf["comment[file uri]"].str.rsplit("/", n=1).str[-1]

    # quantms names each run after its file, and stops on a label-free file listed twice
    repeated = sdrf.loc[sdrf["comment[data file]"].duplicated(), "comment[data file]"].unique()
    if len(repeated) > 0:
        raise ValueError(f"Each data file name can be listed only once, but these are repeated: {', '.join(repeated)}")

    defaults = {
        "assay name": sdrf["source name"],
        "characteristics[organism]": organism,
        "comment[fraction identifier]": "1",
        **{column: "not applicable" for column in NOT_APPLICABLE_COLUMNS},
    }
    for column, default in defaults.items():
        sdrf[column] = sdrf[column].mask(sdrf[column] == "", default)

    # A sample's fractions belong to the same technical replicate
    run_number = sdrf.groupby(["source name", "comment[fraction identifier]"]).cumcount() + 1
    sdrf["comment[technical replicate]"] = sdrf["comment[technical replicate]"].mask(
        sdrf["comment[technical replicate]"] == "", run_number.astype(str)
    )

    # sdrf-pipelines reads both as integers
    for column in ["comment[fraction identifier]", "comment[technical replicate]"]:
        invalid = sdrf.loc[~sdrf[column].str.fullmatch(r"[1-9][0-9]*"), column].unique()
        if len(invalid) > 0:
            raise ValueError(f"{column} must be a whole number from 1, not: {', '.join(invalid)}")

    # sdrf-pipelines stops on any empty cell, and takes conditions from the
    # characteristics columns when there is no factor value column
    with_factor = sdrf[FACTOR] != ""
    if not with_factor.any():
        sdrf = sdrf.drop(columns=[FACTOR])
    elif not with_factor.all():
        rows = ", ".join(str(i + 1) for i in sdrf.index[~with_factor])
        raise ValueError(f"Fill in {FACTOR} in every row or none, but it is blank in row(s) {rows}")

    sdrf = sdrf.assign(**{
        "comment[label]": "NT=label free sample",
        "comment[instrument]": "NT=LTQ Orbitrap XL",
        "comment[cleavage agent details]": f"NT={enzyme}",
        "comment[fragment mass tolerance]": fragment_tolerance,
        "comment[precursor mass tolerance]": precursor_tolerance,
    })
    return sdrf.reindex(columns=sort_columns(sdrf.columns))


def sort_columns(columns) -> list:
    """The column order the SDRF schema requires, with the factor values last."""

    ordered = (
        ["source name"] +
        [cname for cname in columns if cname.startswith("characteristics")] +
        ["assay name"] +
        [cname for cname in columns if cname.startswith("comment")]
    )
    return ordered + [cname for cname in columns if cname not in ordered]


if __name__ == "__main__":

    # Imported here so that test_preprocess.py runs without the SDK
    from cirro.helpers.preprocess_dataset import PreprocessDataset

    ds = PreprocessDataset.from_running()

    sdrf = make_sdrf(
        read_samplesheet(ds.params.get("samplesheet")),
        organism=ds.params["organism"],
        enzyme=ds.params["enzyme"],
        fragment_tolerance=f"{ds.params['fragment_mass_tolerance']} {ds.params['fragment_mass_tolerance_unit']}",
        precursor_tolerance=f"{ds.params['precursor_mass_tolerance']} {ds.params['precursor_mass_tolerance_unit']}",
    )

    # Written apart from the table (mapped in process-input.json): quantms treats its input
    # as an SDRF only when the file name contains "sdrf", and names its outputs after it,
    # which hutch-quantms-qc reads as proteomicslfq/input.sdrf_openms_design_openms.mzTab
    sdrf.to_csv(ds.params["input"], sep="\t", index=None)

    # Used only to build the SDRF, and not quantms params
    ds.remove_param("samplesheet")
    ds.remove_param("organism")
