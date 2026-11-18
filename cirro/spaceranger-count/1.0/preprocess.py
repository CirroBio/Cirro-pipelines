#!/usr/bin/env python3

import pandas as pd


def read_samplesheet(path) -> pd.DataFrame:
    """Read the samplesheet filled out in the form, keeping every value as text."""

    # The form field is marked required, but an empty table can still be submitted,
    # in which case Cirro either omits the param or writes an empty file
    if not path:
        raise ValueError("No samplesheet was provided: add at least one row of FASTQ files")
    try:
        samplesheet = pd.read_csv(path, dtype=str, keep_default_na=False)
    except pd.errors.EmptyDataError:
        samplesheet = pd.DataFrame()
    if samplesheet.empty:
        raise ValueError("The samplesheet is empty: add at least one row of FASTQ files")
    return samplesheet


def make_image_manifest(images: pd.DataFrame, images_dir: str) -> pd.DataFrame:
    """The samplesheet.csv of the image dataset, with each file made a full path in it."""

    for kw in ["sample", "file"]:
        if kw not in images.columns:
            raise ValueError(
                f"The samplesheet.csv of the image dataset must have a '{kw}' column "
                f"(found: {', '.join(images.columns)})"
            )
    return images.assign(file=images["file"].apply(lambda fn: f"{images_dir}/{fn}"))


def check_samples_match(samplesheet: pd.DataFrame, images: pd.DataFrame):
    """Space Ranger runs each sample's FASTQs with its image, so every sample needs both."""

    fastq_samples = set(samplesheet["sample"])
    img_samples = set(images["sample"])
    unmatched = [
        f"{sample} has {what} only"
        for samples, what in [(img_samples - fastq_samples, "images"), (fastq_samples - img_samples, "FASTQs")]
        for sample in sorted(samples)
    ]
    if unmatched:
        raise ValueError(f"All FASTQ data must match up to image files: {'; '.join(unmatched)}")


def reference_paths(ref: str) -> dict:
    """The Space Ranger reference folder and probe set inside the picked reference dataset.

    Built here rather than with a `|/...` suffix in process-input.json: a whole dataset
    picked in the folder picker arrives as its data/ folder with a trailing slash, and
    Nextflow stages a folder whose path contains `//` as an empty one.
    """
    ref = ref.rstrip("/")
    return {"ref": f"{ref}/custom_SpaceRanger_reference/", "probes": f"{ref}/probes.csv"}


if __name__ == "__main__":

    # Imported here so that test_preprocess.py runs without the SDK
    from cirro.helpers.preprocess_dataset import PreprocessDataset

    ds = PreprocessDataset.from_running()

    samplesheet = read_samplesheet(ds.params.get("fastq_manifest"))

    # The image file names are joined onto this prefix. A whole dataset picked in the
    # folder picker arrives as its data/ folder with a trailing slash, which pandas
    # would keep as an empty path segment. It is consumed here and never reaches the
    # run, which is why the workflow need not declare it.
    images_dir = ds.params["images"].rstrip("/")
    ds.logger.info(f"Reading images from {images_dir}")
    images = make_image_manifest(
        pd.read_csv(f"{images_dir}/samplesheet.csv", dtype=str, keep_default_na=False),
        images_dir
    )

    check_samples_match(samplesheet, images)

    for param, value in reference_paths(ds.params["ref"]).items():
        ds.add_param(param, value, overwrite=True)

    # Write to the dataset's config/ folder (mapped in process-input.json)
    images.to_csv(ds.params["image_manifest"], index=None)
    ds.remove_param("images", force=True)

    # Force params.json to be written: the HealthOmics pre-process Lambda fails the run
    # when the file is absent, and the SDK writes it only when a parameter changes.
    ds.keep_params(list(ds.params.keys()))
