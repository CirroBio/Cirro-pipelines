#!/usr/bin/env python3

from cirro.helpers.preprocess_dataset import PreprocessDataset


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

resolve_references(ds, "vdj_dir")

# These name CellRanger reference directories. HealthOmics reads an S3 value with no
# trailing separator as an object key and rejects the run, so mark them as folders.
for ref_dir in ("vdj_dir",):
    ref_path = ds.params.get(ref_dir)
    if ref_path:
        ds.add_param(ref_dir, ref_path.rstrip("/") + "/", overwrite=True)

# Build fastq_dir as a comma-delimited list of all input dataset paths
data_paths = [dataset['dataPath'] for dataset in ds.metadata['inputs']]
assert len(data_paths) > 0, "No input datasets found"
ds.add_param("fastq_dir", ",".join(data_paths))
ds.logger.info(f"fastq_dir: {ds.params['fastq_dir']}")

# Log the parameters present
for k, v in ds.params.items():
    ds.logger.info(f"{k}: {v}")
