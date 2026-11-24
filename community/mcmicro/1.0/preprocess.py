#!/usr/bin/env python3

from typing import TYPE_CHECKING, List

import pandas as pd

if TYPE_CHECKING:
    from cirro.helpers.preprocess_dataset import PreprocessDataset


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


def marker_table(markers: pd.DataFrame) -> pd.DataFrame:
    """The images to stage and their markers.csv rows, in the order MCMICRO reads them.

    MCMICRO registers the images in raw/ sorted by file name and matches markers.csv to
    the resulting channels by row, so the rows follow file name order, and the cycles
    must not go down in that order. Every image is copied into the one raw/ folder, so
    file names must be unique across the table.
    """
    file_names = markers["image"].str.rsplit("/", n=1).str[-1]
    duplicated = sorted(set(file_names[file_names.duplicated()]))
    if duplicated:
        raise ValueError(f"Image file names must be unique, but these repeat: {', '.join(duplicated)}")

    # The form checks the column's type only as advice, so it may still hold any text
    not_whole = markers.loc[~markers["cycle"].str.fullmatch(r"\d+"), "cycle"]
    if not not_whole.empty:
        raise ValueError(f"Cycles must be whole numbers, not: {', '.join(not_whole.unique())}")

    table = (
        markers
        .assign(file_name=file_names, cycle=markers["cycle"].astype(int))
        .sort_values("file_name")
    )
    if not table["cycle"].is_monotonic_increasing:
        order = ", ".join(f"{row.file_name} ({row.cycle})" for row in table.itertuples())
        raise ValueError(
            f"MCMICRO reads the images in file name order, so their cycles must not go down in that order: {order}"
        )
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
