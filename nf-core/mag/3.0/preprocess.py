#!/usr/bin/env python3

from cirro.helpers.preprocess_dataset import PreprocessDataset

ds = PreprocessDataset.from_running()

# Based on run_gtdb, set skip_gtdbtk
ds.add_param(
    "skip_gtdbtk",
    not ds.params.get("run_gtdb", False)
)
ds.remove_param("run_gtdb")

for run_param, db_param in [
    ("run_kraken2", "kraken2_db"),
    ("run_cat", "cat_db"),
    ("run_genomad", "genomad_db")
]:
    # Based on the run param, optionally unset the db param
    if not ds.params.get(run_param, False):
        ds.remove_param(db_param, force=True)
    ds.remove_param(run_param, force=True)

# Set up the samplesheet input
# sample, short_reads_1, short_reads_2
samplesheet = (
    ds.files
    .assign(
        readType=lambda d: d.apply(
            lambda r: r.get("readType", "R"),
            axis=1
        )
    )
    .query("readType == 'R'")
    .pivot(
        index=["sampleIndex", "sample", "dataset"],
        columns="read",
        values="file"
    )
    .rename(columns=lambda i: f"short_reads_{int(i)}")
    .reset_index()
    .reindex(columns=[
        "sample",
        "short_reads_1",
        "short_reads_2",
        "long_reads"
    ])
    .merge(
        ds.samplesheet,
        left_on="sample",
        right_on="sample"
    )
    .assign(
        group=lambda d: d.apply(
            lambda r: r.get("group", 0),
            axis=1
        )
    )
)

ds.logger.info("Formatted samplesheet:")
ds.logger.info(samplesheet.to_csv(index=None))
msg = "No files detected -- there may be an error with data ingest"
assert samplesheet.shape[0] > 0, msg

# Write to the dataset's config/ folder (mapped in process-input.json)
samplesheet.to_csv(ds.params["input"], index=None)

# Auto-detect single-end
if samplesheet["short_reads_2"].isnull().all():
    ds.add_param("single_end", True)
