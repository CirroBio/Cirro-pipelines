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

    fastq_df.to_csv("fastq_manifest.csv", index=None)
    ds.add_param("fastq_manifest", "fastq_manifest.csv")
    img_df.to_csv("image_manifest.csv", index=None)
    ds.add_param("image_manifest", "image_manifest.csv")


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
    fastq_df = (
        ds.files
        .assign(
            readType=ds.files.reindex(columns=["readType"])["readType"].fillna("R")
        )
        .query("readType == 'R'")
        .reindex(columns=["sampleIndex", "sample", "lane", "read", "file"])
        .pivot(
            index=["sampleIndex", "sample", "lane"],
            columns="read",
            values="file"
        )
        .rename(columns=lambda i: f"fastq_{int(i)}")
        .reset_index()
        .merge(ds.samplesheet, on="sample")
    )
    ds.logger.info("Creating paired table of inputs - DONE")
    ds.logger.info(fastq_df.to_csv(index=None))

    return fastq_df


if __name__ == "__main__":
    main()
