#!/usr/bin/env python3

import re
from typing import TYPE_CHECKING, List

import pandas as pd

if TYPE_CHECKING:
    from cirro.helpers.preprocess_dataset import PreprocessDataset

# Markers on channels that carry no stain; their images are left out of the run
BLANK_MARKERS = ["Blank", "Empty"]


def read_markers(path) -> pd.DataFrame:
    """Read the images and markers table filled out in the form, keeping every value as text."""

    # The form field is marked required, but an empty table can still be submitted, in
    # which case Cirro either omits the param or writes an empty file.
    if not path:
        raise ValueError("No markers table was provided: add a row for each image")
    try:
        markers = pd.read_csv(path, dtype=str, keep_default_na=False)
    except pd.errors.EmptyDataError:
        markers = pd.DataFrame()
    if markers.empty:
        raise ValueError("The markers table is empty: add a row for each image")
    return markers


def image_cycle(file_name: str) -> int:
    """The cycle number in an image's file name, such as 1 for cyc001_reg001.tif."""
    match = re.match(r"cyc(\d+)_", file_name)
    if match is None:
        raise ValueError(
            f"Could not read a cycle number from {file_name}: expected a name such as cyc001_reg001.tif"
        )
    return int(match.group(1))


def marker_table(markers: pd.DataFrame) -> pd.DataFrame:
    """The images to stage and their markers.csv rows, in the order MCMICRO reads them.

    MCMICRO registers the images in raw/ sorted by file name and matches markers.csv to
    the resulting channels by row, so the rows follow file name order. Every image is
    copied into the one raw/ folder, so file names must be unique across the table.
    """
    file_names = markers["image"].str.rsplit("/", n=1).str[-1]
    duplicated = sorted(set(file_names[file_names.duplicated()]))
    if duplicated:
        raise ValueError(f"Image file names must be unique, but these repeat: {', '.join(duplicated)}")

    table = (
        markers
        .assign(file_name=file_names, cycle=[image_cycle(name) for name in file_names])
        .sort_values("file_name")
    )
    table = table[~table["marker_name"].isin(BLANK_MARKERS)]
    if table.empty:
        raise ValueError(f"Every image in the markers table is marked {' or '.join(BLANK_MARKERS)}")
    return table.reindex(columns=["image", "cycle", "marker_name"]).reset_index(drop=True)


def stage_images(images: List[str], raw_folder: str):
    """Copy the images into raw/ of the MCMICRO project folder, which is also the output folder."""

    # Imported here so that test_preprocess.py runs without the SDK or boto3
    import boto3
    from cirro.models.s3_path import S3Path

    s3 = boto3.client("s3")
    for image in images:
        source = S3Path(image)
        dest = S3Path(f"{raw_folder}/{image.rsplit('/', 1)[-1]}")
        s3.copy({"Bucket": source.bucket, "Key": source.key}, dest.bucket, dest.key)


def setup_params(ds: "PreprocessDataset"):
    """Set up the parameters used for execution."""

    # Set up the workflow object
    workflow = dict()

    # Tissue Microarray vs. Whole-Slide Image
    ds.logger.info(f"Image Type: {ds.params['image_type']}")
    if ds.params["image_type"] == "Tissue Microarray":
        workflow["tma"] = True
    else:
        assert ds.params["image_type"] == "Whole-Slide Image"
        workflow["tma"] = False
    ds.remove_param("image_type")

    # Segmentation options
    setup_param_list(
        ds,
        workflow,
        "segmentation",
        ["unmicst", "ilastik", "cypository"]
    )

    # Downstream modules options
    setup_param_list(
        ds,
        workflow,
        "downstream",
        ["naivestates", "scimap", "fastpg", "scanpy", "flowsom"]
    )
    if len(workflow["downstream"]) > 0:
        ds.add_param("stop-at", "downstream")

    ds.add_param("workflow", workflow)


def setup_param_list(
    ds: "PreprocessDataset",
    workflow: dict,
    list_key: str,
    list_items: List[str]
):

    workflow[list_key] = [
        i for i in list_items if ds.params.get(i, False)
    ]

    for i in list_items:
        ds.remove_param(i, force=True)


if __name__ == "__main__":

    # Imported here so that test_preprocess.py runs without the SDK
    from cirro.helpers.preprocess_dataset import PreprocessDataset

    ds = PreprocessDataset.from_running()

    markers = marker_table(read_markers(ds.params.get("markers")))
    ds.remove_param("markers")
    ds.logger.info(f"Staging {len(markers):,} images")

    # MCMICRO reads markers.csv and raw/ from its project folder (--in) and writes its
    # outputs beside them, so the project folder is the output folder
    project_folder = ds.params["in"].rstrip("/")
    stage_images(markers["image"].tolist(), f"{project_folder}/raw")
    markers.drop(columns=["image"]).to_csv(f"{project_folder}/markers.csv", index=None)

    setup_params(ds)

    ds.logger.info(ds.params)
