import pandas as pd


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


def drop_blank_fastq_2(samplesheet: pd.DataFrame) -> pd.DataFrame:
    """Remove the R2 column of a single-end samplesheet.

    gig-map calls file() on every fastq_2 cell, which fails on a blank one, and skips
    R2 only when the column is absent. So R2 must be given for every row or none.
    """
    blank = samplesheet["fastq_2"] == ""
    if blank.all():
        return samplesheet.drop(columns="fastq_2")
    if blank.any():
        samples = ", ".join(samplesheet.loc[blank, "sample"].unique())
        raise ValueError(
            f"Some rows have no FASTQ (R2) file ({samples}): "
            "give every row an R2 file, or leave R2 blank in all of them"
        )
    return samplesheet


def pangenome_paths(pangenome: str) -> dict:
    """The gig-map inputs inside the picked pangenome dataset.

    Built here rather than with a `|/...` suffix in process-input.json: a whole dataset
    picked in the folder picker arrives as its data/ folder with a trailing slash.
    """
    pangenome = pangenome.rstrip("/")
    return {
        "genes": f"{pangenome}/gene_catalog/centroids.faa.gz",
        "gene_bins": f"{pangenome}/bin_pangenome/gene_bins.csv",
        "group_profile": f"{pangenome}/bin_pangenome/group_profile.csv",
        "genome_groups": f"{pangenome}/bin_pangenome/genome_groups.csv",
    }


def check_metadata(metadata: pd.DataFrame, formula):
    """The metadata is only read after every sample is aligned, so check its columns first."""

    missing = [
        cname for cname in ["sample"] + (formula or "").replace("+", " ").split()
        if cname not in metadata.columns
    ]
    if missing:
        raise ValueError(
            f"The sample metadata has no column {', '.join(missing)} "
            f"(found: {', '.join(metadata.columns)})"
        )


if __name__ == "__main__":

    # Imported here so that test_preprocess.py runs without the SDK
    from cirro.helpers.preprocess_dataset import PreprocessDataset

    ds = PreprocessDataset.from_running()

    samplesheet = drop_blank_fastq_2(read_samplesheet(ds.params.get("samplesheet")))
    samplesheet.to_csv(ds.params["samplesheet"], index=None)

    for param, value in pangenome_paths(ds.params["pangenome"]).items():
        ds.add_param(param, value)
    ds.remove_param("pangenome")

    # Without a metadata file gig-map falls back to its own placeholder table
    if ds.params.get("metadata"):
        check_metadata(
            pd.read_csv(ds.params["metadata"], dtype=str, keep_default_na=False),
            ds.params.get("formula")
        )

    # Force params.json to be written: the HealthOmics pre-process Lambda fails the run
    # when the file is absent, and the SDK writes it only when a parameter changes.
    ds.keep_params(list(ds.params.keys()))
