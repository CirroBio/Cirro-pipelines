#!/usr/bin/env python3

from cirro.helpers.preprocess_dataset import PreprocessDataset
import pandas as pd


def main():
    ds = PreprocessDataset.from_running()

    # Get the table of FASTQ inputs
    fastq_df = make_fastq_df(ds)

    # Get the table of images
    img_df = make_img_df(ds)

    # Make sure that the samples line up
    fastq_samples = set(fastq_df["sample"].tolist())
    img_samples = set(img_df["sample"].tolist())

    ds.logger.info(f"Number of samples with FASTQs: {len(fastq_samples):,}")
    ds.logger.info(f"Number of samples with images: {len(img_samples):,}")

    valid = True
    for diff, msg in [
        (img_samples - fastq_samples, "images only"),
        (fastq_samples - img_samples, "FASTQs only"),
    ]:
        if len(diff) > 0:
            valid = False
            for n in list(diff):
                ds.logger.info(f"Sample {n} has {msg}")

    assert valid, "All FASTQ data must match up to image files"

    # Write to the dataset's config/ folder (mapped in process-input.json)
    fastq_df.to_csv(ds.params["fastq_manifest"], index=None)
    img_df.to_csv(ds.params["image_manifest"], index=None)

    # make_img_df joins the samplesheet name and each image name onto this prefix,
    # so it is deliberately left without a trailing separator. It is consumed here
    # and never reaches the run, which is why the workflow need not declare it.
    ds.remove_param("images", force=True)

    # Force params.json to be written: the HealthOmics pre-process Lambda fails the run
    # when the file is absent, and the SDK writes it only when a parameter changes.
    ds.keep_params(list(ds.params.keys()))


def make_img_df(ds: PreprocessDataset):
    img_prefix = ds.params["images"]
    ds.logger.info(f"Reading images from {img_prefix}")
    ds.logger.info("Loading samplesheet.csv")
    img_df = pd.read_csv(f"{img_prefix}/samplesheet.csv")
    for kw in ["sample", "file"]:
        assert kw in img_df.columns.values, f"Expected '{kw}' column"
    # Add the full path to the file
    img_df = img_df.assign(file=img_df["file"].apply(lambda fn: f"{img_prefix}/{fn}"))
    ds.logger.info(img_df.to_csv(index=None))
    return img_df


def make_fastq_df(ds: PreprocessDataset) -> pd.DataFrame:
    """Format the FASTQ inputs in wide format."""

    # Format as a wide dataset
    ds.logger.info("Creating paired table of FASTQ inputs")
    reads = (
        ds.files
        .assign(
            readType=ds.files.reindex(columns=["readType"])["readType"].fillna("R")
        )
        .query("readType == 'R'")
        .reindex(columns=["sampleIndex", "sample", "lane", "read", "file"])
    )

    # A run folder holds more than the sample FASTQs -- QC reports, undetermined
    # reads -- and ingest parses no read number from those. readType is filled in
    # above, so they pass the filter with read unset and pivot into a NaN column,
    # which cannot be named.
    unparsed = reads["read"].isna()
    if unparsed.any():
        ds.logger.info(f"Ignoring {unparsed.sum():,} file(s) with no read number:")
        for fn in reads.loc[unparsed, "file"]:
            ds.logger.info(f"  {fn}")
        reads = reads.loc[~unparsed]

    fastq_df = (
        reads
        .pivot(
            index=["sampleIndex", "sample", "lane"],
            columns="read",
            values="file"
        )
        .rename(columns=lambda i: f"fastq_{int(i)}")
        .reset_index()
        .merge(ds.samplesheet, on="sample")
    )
    assert fastq_df.shape[0] > 0, \
        "No paired FASTQ files found -- there may be an error with data ingest"

    ds.logger.info("Creating paired table of inputs - DONE")
    ds.logger.info(fastq_df.to_csv(index=None))

    return fastq_df


if __name__ == "__main__":
    main()
