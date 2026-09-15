#!/usr/bin/env python3

from cirro.helpers.preprocess_dataset import PreprocessDataset
import pandas as pd


def make_manifest(ds: PreprocessDataset) -> pd.DataFrame:
    """Construct a manifest with the paired FASTQ files in the input"""

    ds.logger.info(f"Number of files in dataset: {ds.files.shape[0]:,}")
    ds.logger.info(f"Number of samples in dataset: {ds.samplesheet.shape[0]:,}")

    # Build the manifest
    manifest = (
        ds.files
        .reindex(columns=["sampleIndex", "libraryIndex", "sample", "dataset", "readType", "read", "file"])
        # A row with no read number cannot be paired, and would pivot
        # into a column the rename below cannot name.
        .loc[lambda d: d["read"].notna()]
        .assign(
            fastq_cname_ix=lambda df: df.apply(
                lambda r: r["read"] + (2 if r.get("readType", "R") == "I" else 0),
                axis=1
            )
        )
        .pivot(
            index=["sampleIndex", "libraryIndex", "sample", "dataset"],
            columns="fastq_cname_ix",
            values="file"
        )
        .rename(columns=lambda i: f"fastq_{int(i)}")
        .sort_index()
        .sort_index(axis=1)
        .reset_index()
        .drop(columns=["sampleIndex", "libraryIndex", "dataset"])
    )

    log_table(ds, "Manifest", manifest)

    # Get the read_structure attribute for each sample (if any exists)
    read_structure = ds.samplesheet.reindex(
        columns=["sample", "read_structure"]
    ).set_index(
        "sample"
    )["read_structure"].dropna()

    # Log the number of samples which have read_structure defined
    if read_structure.shape[0] == 0:
        ds.logger.info("No samples have read_structure defined by the samplesheet")
    else:
        for sample, rs in read_structure.items():
            ds.logger.info(f"Sample: {sample} - Read Structure: {rs}")

    ds.logger.info(f"Any missing read_structures will be populated with {ds.params.get('read_structure', 'Missing')}")

    # Add that information to the manifest, filling in the form value when missing
    manifest = manifest.assign(
        read_structure=manifest["sample"].apply(
            read_structure.get
        ).fillna(
            ds.params.get("read_structure")
        )
    )
    assert not manifest["read_structure"].isnull().any(), "Missing read structure information"
    log_table(ds, "Manifest", manifest)

    return manifest


def log_table(ds: PreprocessDataset, title: str, df: pd.DataFrame):
    ds.logger.info(f"{title}:")
    for line in df.to_csv(index=None).split("\n"):
        ds.logger.info(line)


if __name__ == "__main__":
    ds = PreprocessDataset.from_running()

    manifest = make_manifest(ds)

    # Write to the dataset's config/ folder (mapped in process-input.json)
    manifest.to_csv(ds.params["input"], index=None)
    ds.logger.info(f"Wrote out {manifest.shape[0]:,} lines to {ds.params['input']}")
