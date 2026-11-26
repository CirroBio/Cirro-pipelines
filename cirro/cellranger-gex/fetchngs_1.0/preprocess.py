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


def order_reads(samplesheet: pd.DataFrame, sizes: dict) -> tuple[pd.DataFrame, bool]:
    """Swap the FASTQ columns when fastq_1 holds the larger files, and say whether it did.

    SRA numbers a run's reads in upload order, so `_1` is not always the barcode read.
    R2 (cDNA) is the longer read, so its files are larger. This compares the median size
    of each column, following
    https://kb.10xgenomics.com/hc/en-us/articles/115003802691-How-do-I-prepare-Sequence-Read-Archive-SRA-data-from-NCBI-for-Cell-Ranger-
    """
    median_1 = samplesheet["fastq_1"].map(sizes).median()
    median_2 = samplesheet["fastq_2"].map(sizes).median()
    if median_1 <= median_2:
        return samplesheet, False
    return samplesheet.rename(columns={"fastq_1": "fastq_2", "fastq_2": "fastq_1"})[samplesheet.columns], True


def file_size(s3, uri: str) -> int:
    bucket, key = uri[len("s3://"):].split("/", 1)
    return s3.head_object(Bucket=bucket, Key=key)["ContentLength"]


if __name__ == "__main__":

    # Imported here so that test_preprocess.py runs without the SDK
    import boto3
    from cirro.helpers.preprocess_dataset import PreprocessDataset

    ds = PreprocessDataset.from_running()

    transcriptome_dir = ds.params.get("transcriptome_dir")
    resolved = absolute_reference(transcriptome_dir, ds.references_base)
    if resolved != transcriptome_dir:
        ds.add_param("transcriptome_dir", resolved, overwrite=True)

    samplesheet = read_samplesheet(ds.params.get("samplesheet"))

    s3 = boto3.client("s3")
    uris = pd.concat([samplesheet["fastq_1"], samplesheet["fastq_2"]]).unique()
    samplesheet, swapped = order_reads(samplesheet, {uri: file_size(s3, uri) for uri in uris})
    if swapped:
        ds.logger.info("fastq_1 held the larger files, so it was swapped with fastq_2 to make R1 the barcode read")
        samplesheet.to_csv(ds.params["samplesheet"], index=None)

    # Force params.json to be written: the HealthOmics pre-process Lambda fails the run
    # when the file is absent, and the SDK writes it only when a parameter changes.
    ds.keep_params(list(ds.params.keys()))
