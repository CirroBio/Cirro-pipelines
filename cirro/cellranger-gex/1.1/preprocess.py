from cirro.helpers.preprocess_dataset import PreprocessDataset
import pandas as pd


def resolve_references(ds: PreprocessDataset, *params: str):
    """Make reference params absolute against the references bucket.

    Form values name a path under that bucket rather than a full URI, so the
    configuration is not tied to one deployment. A value that is already absolute
    is left alone, which keeps datasets created before that change re-runnable.
    """
    for param in params:
        value = ds.params.get(param)
        if not isinstance(value, str) or not value or value.startswith("s3://"):
            continue
        ds.add_param(param, f"{ds.references_base}/{value}", overwrite=True)

ds = PreprocessDataset.from_running()

resolve_references(ds, "transcriptome_dir")

# The form field is marked required, but an empty table can still be submitted, in
# which case Cirro either omits the param or writes an empty file. The workflow would
# then finish successfully without running Cell Ranger at all.
if not ds.params.get("samplesheet"):
    raise ValueError("No samplesheet was provided: add at least one row of FASTQ files")
try:
    samplesheet = pd.read_csv(ds.params["samplesheet"], dtype=str, keep_default_na=False)
except pd.errors.EmptyDataError:
    raise ValueError("The samplesheet is empty: add at least one row of FASTQ files")

ds.logger.info("Samplesheet:")
ds.logger.info(samplesheet.to_csv(index=None))

# Force params.json to be written: the HealthOmics pre-process Lambda fails the run
# when the file is absent, and the SDK writes it only when a parameter changes.
ds.keep_params(list(ds.params.keys()))
