#!/usr/bin/env python3

from __future__ import annotations

from typing import TYPE_CHECKING

import pandas as pd

# The SDK is imported in __main__, so test_preprocess.py runs without it
if TYPE_CHECKING:
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


def read_samplesheet(path) -> pd.DataFrame:
    """Read the samplesheet filled out in the form, keeping every value as text."""

    # The form field is marked required, but an empty table can still be submitted,
    # in which case Cirro either omits the param or writes an empty file
    if not path:
        raise ValueError("No samplesheet was provided: add at least one row of FASTQ files")
    try:
        # Read every value as text, so that sample names like "01" are kept as-is
        samplesheet = pd.read_csv(path, dtype=str, keep_default_na=False)
    except pd.errors.EmptyDataError:
        samplesheet = pd.DataFrame()
    if samplesheet.empty:
        raise ValueError("The samplesheet is empty: add at least one row of FASTQ files")
    return samplesheet


def fill_group(samplesheet: pd.DataFrame) -> pd.DataFrame:
    """Put every row that leaves group blank in group 0; the pipeline requires a group."""

    return samplesheet.assign(group=samplesheet["group"].replace("", "0"))


def check_unique_runs(samplesheet: pd.DataFrame):
    """Fail if two rows of a sample share a run ID, which the pipeline rejects.

    The pipeline reads a blank run as run 0, so a blank and a 0 also collide.
    """

    keys = samplesheet[["sample", "run"]].replace({"run": {"": "0"}})
    duplicated = keys.loc[keys.duplicated(keep=False), "sample"].unique()
    if len(duplicated) > 0:
        raise ValueError(
            f"Each row of a sample needs a different run ID; check sample(s) {', '.join(duplicated)}"
        )


def is_single_end(samplesheet: pd.DataFrame) -> bool:
    """Whether the short reads are single-end, which the pipeline cannot mix with paired-end."""

    paired = samplesheet["short_reads_2"] != ""
    if paired.all():
        return False
    if not paired.any():
        return True
    raise ValueError(
        "Short reads must be all single-end or all paired-end: "
        "fill in Short reads (R2) for every row, or for none"
    )


if __name__ == "__main__":

    from cirro.helpers.preprocess_dataset import PreprocessDataset

    ds = PreprocessDataset.from_running()

    resolve_references(ds, "kraken2_db", "metaeuk_db", "busco_db")

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

    samplesheet = fill_group(read_samplesheet(ds.params.get("input")))
    ds.logger.info("Samplesheet:")
    ds.logger.info(samplesheet.to_csv(index=None))
    check_unique_runs(samplesheet)

    # Overwrite the file Cirro wrote to the dataset's config/ folder (mapped in process-input.json)
    samplesheet.to_csv(ds.params["input"], index=None)

    if is_single_end(samplesheet):
        ds.add_param("single_end", True)
