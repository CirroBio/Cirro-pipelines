import pandas as pd


def absolute_reference(value, references_base: str):
    """A reference param's value, made absolute against the references bucket.

    Form values name a path under that bucket rather than a full URI, so the
    configuration is not tied to one deployment. A value that is already absolute
    is returned as is, which keeps datasets created before that change re-runnable.
    """
    if not isinstance(value, str) or not value or value.startswith("s3://"):
        return value
    return f"{references_base}/{value}"


def read_samplesheet(path) -> pd.DataFrame:
    """Read the samplesheet filled out in the form, keeping every value as text."""

    # The form field is marked required, but an empty table can still be submitted, in
    # which case Cirro either omits the param or writes an empty file. The workflow would
    # then finish successfully without running Cell Ranger at all.
    if not path:
        raise ValueError("No samplesheet was provided: add at least one row of FASTQ files")
    try:
        samplesheet = pd.read_csv(path, dtype=str, keep_default_na=False)
    except pd.errors.EmptyDataError:
        samplesheet = pd.DataFrame()
    if samplesheet.empty:
        raise ValueError("The samplesheet is empty: add at least one row of FASTQ files")
    return samplesheet


if __name__ == "__main__":

    # Imported here so that test_preprocess.py runs without the SDK
    from cirro.helpers.preprocess_dataset import PreprocessDataset

    ds = PreprocessDataset.from_running()

    transcriptome_dir = ds.params.get("transcriptome_dir")
    resolved = absolute_reference(transcriptome_dir, ds.references_base)
    if resolved != transcriptome_dir:
        ds.add_param("transcriptome_dir", resolved, overwrite=True)

    # Read only to fail clearly when the table is missing or empty
    read_samplesheet(ds.params.get("samplesheet"))

    # Force params.json to be written: the HealthOmics pre-process Lambda fails the run
    # when the file is absent, and the SDK writes it only when a parameter changes.
    ds.keep_params(list(ds.params.keys()))
