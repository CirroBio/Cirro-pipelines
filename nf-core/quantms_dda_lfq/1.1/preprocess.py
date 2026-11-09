#!/usr/bin/env python3

import re

import pandas as pd

FACTOR = "factor value[condition]"
MODIFICATION = "comment[modification parameters]"

# The Unimod accession of each modification the form offers. sdrf-pipelines takes the
# modification's name for the search from its accession, and validates the NT= name
# against the Unimod ontology.
UNIMOD_ACCESSIONS = {
    "Acetyl": 1,
    "Carbamidomethyl": 4,
    "Carbamyl": 5,
    "Deamidated": 7,
    "Phospho": 21,
    "Propionamide": 24,
    "Pyro-carbamidomethyl": 26,
    "Glu->pyro-Glu": 27,
    "Gln->pyro-Glu": 28,
    "Methyl": 34,
    "Oxidation": 35,
    "Dimethyl": 36,
    "Trimethyl": 37,
    "Methylthio": 39,
    "GG": 121,
}

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


def split_modifications(value: str) -> list:
    """The modifications picked in the form, which Cirro joins with commas."""

    return [modification.strip() for modification in value.split(",") if modification.strip()]


def sdrf_modification(modification: str, mod_type: str) -> str:
    """
    The SDRF modification parameters for a modification in OpenMS notation, such as
    'Oxidation (M)' or 'Gln->pyro-Glu (N-term Q)'. sdrf-pipelines turns them back into
    that notation for the search engines.
    """

    match = re.fullmatch(r"(.+) \((.+)\)", modification)
    if match is None or match.group(1) not in UNIMOD_ACCESSIONS:
        raise ValueError(f"Unsupported modification: {modification}")
    name, site = match.groups()

    if site in ("Protein N-term", "Protein C-term"):
        position, target = site, site.split()[-1]
    elif site.startswith(("N-term", "C-term")):
        terminus, _, residue = site.partition(" ")
        position, target = f"Any {terminus}", residue or terminus
    else:
        position, target = "Anywhere", site

    return f"NT={name};AC=UNIMOD:{UNIMOD_ACCESSIONS[name]};PP={position};TA={target};MT={mod_type}"


def make_sdrf(
    samplesheet: pd.DataFrame,
    organism: str,
    enzyme: str,
    fragment_tolerance: str,
    precursor_tolerance: str,
    fixed_mods: list,
    variable_mods: list
) -> pd.DataFrame:
    """The SDRF quantms reads: the table with its blanks filled, plus the per-run columns."""

    both = [modification for modification in fixed_mods if modification in variable_mods]
    if both:
        raise ValueError(f"A modification can be fixed or variable, but not both: {', '.join(both)}")

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
    sdrf = sdrf.reindex(columns=sort_columns(sdrf.columns))

    # One column per modification, all with the same name, ahead of the factor values.
    # Without them the search uses only Carbamidomethyl (C), as a fixed modification,
    # whatever quantms's params say.
    modifications = (
        [sdrf_modification(modification, "Fixed") for modification in fixed_mods] +
        [sdrf_modification(modification, "Variable") for modification in variable_mods]
    )
    position = len([cname for cname in sdrf.columns if not cname.startswith("factor value")])
    for i, parameters in enumerate(modifications):
        sdrf.insert(position + i, MODIFICATION, parameters, allow_duplicates=True)
    return sdrf


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
        # Absent or blank when no box is ticked
        fixed_mods=split_modifications(ds.params.get("fixed_mods", "")),
        variable_mods=split_modifications(ds.params.get("variable_mods", "")),
    )

    # Written apart from the table (mapped in process-input.json): quantms treats its input
    # as an SDRF only when the file name contains "sdrf", and names its outputs after it,
    # which hutch-quantms-qc reads as proteomicslfq/input.sdrf_openms_design_openms.mzTab
    sdrf.to_csv(ds.params["input"], sep="\t", index=None)

    # Used only to build the SDRF, and not quantms params
    ds.remove_param("samplesheet")
    ds.remove_param("organism")
