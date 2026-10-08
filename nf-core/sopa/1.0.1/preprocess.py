#!/usr/bin/env python3
"""Build the nf-core/sopa samplesheet for a Cirro launch.

sopa's --input is a CSV with columns sample,data_path where data_path is the
raw machine-output directory for one sample -- except for phenocycler and
ome_tif, where data_path must be the all-channel image file itself, which the
user picks in the form. Cirro selects a dataset rather than a samplesheet, so
the CSV is written here from the selected dataset's data path.

A dataset can hold several samples: a Xenium, MERSCOPE or CosMx run ingested
whole has one output folder per region, found by the file its reader opens
first, and the OME-TIFF picker accepts several images. Each becomes a row.
"""

from fnmatch import fnmatch
from pathlib import PurePosixPath

import boto3
import pandas as pd
from cirro.helpers.preprocess_dataset import PreprocessDataset
from cirro.models.s3_path import S3Path

# Technologies whose data_path is an image file picked in the form
SINGLE_FILE_TECHNOLOGIES = ("phenocycler", "ome_tif")

# Stripped off an image file name to get the sample name. '.ome' is part of the
# extension of an OME-TIFF, so scan.ome.tif is sample 'scan', not 'scan.ome'.
IMAGE_EXTENSIONS = (".qptiff", ".tiff", ".tif", ".ome")

# The file at the top of one output folder, which the sopa reader opens from
# data_path. The CosMx reader insists on exactly one fov_positions file under
# data_path, so a dataset with two slides only reads when split this way.
OUTPUT_MARKERS = {
    "xenium": ("experiment.xenium",),
    "merscope": ("detected_transcripts.csv",),
    "cosmx": ("*_fov_positions_file.csv", "*_fov_positions_file.csv.gz"),
}


def find_output_folders(data_path: str, technology: str) -> list[str]:
    """Return every folder under data_path that holds the technology's marker file."""
    patterns = OUTPUT_MARKERS[technology]
    s3_path = S3Path(data_path)
    assert s3_path.valid, f"Not an S3 URI: {data_path}"
    prefix = s3_path.key.rstrip("/") + "/"
    paginator = boto3.client("s3").get_paginator("list_objects_v2")

    found = {
        f"s3://{s3_path.bucket}/{PurePosixPath(obj['Key']).parent}"
        for page in paginator.paginate(Bucket=s3_path.bucket, Prefix=prefix)
        for obj in page.get("Contents", [])
        # Skip hidden files, e.g. the ._ copies macOS leaves when copying folders
        if not (name := PurePosixPath(obj["Key"]).name).startswith(".")
        and any(fnmatch(name, pattern) for pattern in patterns)
    }
    if not found:
        raise ValueError(
            f"No {' or '.join(patterns)} found under {data_path} - the dataset "
            f"does not contain {technology} output"
        )
    return sorted(found)


def folder_sample_name(folder: str) -> str:
    sample = PurePosixPath(folder).name or "sample"
    if sample == "data":
        # Cirro data folders are all named 'data'; use the parent (dataset id)
        sample = PurePosixPath(folder).parent.name or "sample"
    return sample


def image_sample_name(image_file: str) -> str:
    sample = PurePosixPath(image_file).name
    while (ext := PurePosixPath(sample).suffix.lower()) in IMAGE_EXTENSIONS:
        sample = sample[:-len(ext)]
    return sample or "sample"


def build_samplesheet(
    data_path: str,
    technology: str,
    image_files: list[str] | None = None,
    output_folders: list[str] | None = None,
) -> pd.DataFrame:
    """Return the samplesheet for the selected dataset.

    data_path: the input dataset's data folder (no trailing slash required)
    technology: the sopa reader in use
    image_files: the images picked in the form, for the single-file technologies
    output_folders: the output folders found in the dataset, for OUTPUT_MARKERS
    """
    data_path = data_path.rstrip("/")

    if technology in SINGLE_FILE_TECHNOLOGIES:
        if not image_files:
            raise ValueError(
                f"technology={technology} needs an image file to be selected in the form"
            )
        samplesheet = pd.DataFrame([
            dict(sample=image_sample_name(image_file), data_path=image_file)
            for image_file in image_files
        ])
    else:
        samplesheet = pd.DataFrame([
            dict(sample=folder_sample_name(folder), data_path=folder)
            for folder in (output_folders if technology in OUTPUT_MARKERS else [data_path])
        ])

    duplicated = samplesheet.loc[samplesheet["sample"].duplicated(keep=False)]
    if not duplicated.empty:
        raise ValueError(
            "Sample names must be unique, as sopa names each sample's outputs after it:\n"
            + duplicated.to_csv(index=None)
        )
    return samplesheet


if __name__ == "__main__":

    ds = PreprocessDataset.from_running()

    technology = ds.params.get("technology", "xenium")
    # A multiple-file picker sends its selection as one comma-joined string
    image_files = [
        image_file.strip()
        for image_file in (ds.params.get("image_file") or "").split(",")
        if image_file.strip()
    ]
    output_folders = (
        find_output_folders(ds.params["spatial_data"], technology)
        if technology in OUTPUT_MARKERS else None
    )
    samplesheet = build_samplesheet(
        ds.params["spatial_data"],
        technology,
        image_files,
        output_folders,
    )
    ds.logger.info("samplesheet:")
    ds.logger.info(samplesheet.to_csv(index=None))

    # Write the samplesheet to the S3 location which process-input.json maps the
    # input param to -- the dataset's own config/ folder, as nf-core/ampliseq does.
    # That records the exact input to the run alongside its results, instead of
    # leaving it in the working directory this script runs in.
    samplesheet.to_csv(ds.params["input"], index=None)

    # These exist only to carry values into this script, and are not sopa params
    ds.remove_param("spatial_data", force=True)
    ds.remove_param("image_file", force=True)
