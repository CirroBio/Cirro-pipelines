#!/usr/bin/env python3

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


def log_lines(title: str, dat: pd.DataFrame, **kwargs):
    ds.logger.info(title)
    for line in dat.to_csv(**kwargs).split("\n"):
        ds.logger.info(line)


resolve_references(ds, "reference_dir")

# Log inputs
log_lines("Files:", ds.files, index=None)
assert ds.files.shape[0] > 0, "No files found"

# Remove index reads, if present
if (
    ("readType" in ds.files.columns.values) and
    (ds.files["readType"] == "I").any()
):
    print("Removing index reads")
    ds.files = ds.files.query("readType != 'I'")
    log_lines("Filtered Files:", ds.files, index=None)

log_lines("Samples:", ds.samplesheet, index=None)
assert ds.samplesheet.shape[0] > 0, "No sample metadata found"

# Merge the files table with the sample metadata table
df = (
    ds
    .wide_samplesheet()
    .merge(
        ds.samplesheet,
        left_on="sample",
        right_on="sample"
    )
)
log_lines("Merged:", df)

# The column 'grouping' is used to pair samples
# The column 'feature_types' is used to indicate analysis modality
for cname in ['grouping', 'feature_types']:
    msg = f"'{cname}' metadata required (see documentation)"
    assert cname in df.columns.values, msg
    msg = f"Missing values for '{cname}' metadata (see documentation)"
    assert df[cname].notnull().all(), msg

allowed_types = ['Gene Expression', 'Chromatin Accessibility']
msg = f"Allowed feature_types are '{', '.join(allowed_types)}'"
assert df['feature_types'].isin(allowed_types).all(), msg

# Reformat the samplesheet
df = (
    df.replace(
        to_replace=dict(
            feature_types=dict(zip(
                allowed_types,
                ['gex', 'atac']
            ))
        )
    )
    .reindex(
        columns=['grouping', 'feature_types', 'fastq_1', 'fastq_2']
    )
    .rename(
        columns=dict(
            grouping='sample',
            feature_types='fastq_type'
        )
    )
)
log_lines("Reformatted:", df, index=None)

# Write to the dataset's config/ folder (mapped in process-input.json)
df.to_csv(ds.params["samplesheet"], index=None)

# Force params.json to be written: the HealthOmics pre-process Lambda fails the run
# when the file is absent, and the SDK writes it only when a parameter changes.
ds.keep_params(list(ds.params.keys()))
