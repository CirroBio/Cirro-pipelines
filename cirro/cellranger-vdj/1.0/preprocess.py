#!/usr/bin/env python3


def reference_dir(value: str, references_base: str) -> str:
    """The V(D)J reference param's value, made absolute and marked as a folder.

    Form values name a path under the references bucket rather than a full URI, so the
    configuration is not tied to one deployment. A value that is already absolute is
    kept, which keeps datasets created before that change re-runnable. HealthOmics reads
    an S3 value with no trailing separator as an object key and rejects the run.
    """
    if not value.startswith("s3://"):
        value = f"{references_base}/{value}"
    return value.rstrip("/") + "/"


def fastq_folder(value) -> str:
    """The picked FASTQ folder, without the trailing slash a whole-dataset pick has.

    The workflow joins `/*` and `/**.fastq.gz` onto it, and would otherwise glob a path
    containing `//`.
    """
    if not value:
        raise ValueError("No FASTQ folder was selected")
    return value.rstrip("/")


if __name__ == "__main__":

    # Imported here so that test_preprocess.py runs without the SDK
    from cirro.helpers.preprocess_dataset import PreprocessDataset

    ds = PreprocessDataset.from_running()

    ds.add_param("vdj_dir", reference_dir(ds.params["vdj_dir"], ds.references_base), overwrite=True)
    ds.add_param("fastq_dir", fastq_folder(ds.params.get("fastq_dir")), overwrite=True)

    # Log the parameters present
    for k, v in ds.params.items():
        ds.logger.info(f"{k}: {v}")
