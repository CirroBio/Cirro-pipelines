#!/usr/bin/env python3

from typing import List
from cirro.helpers.preprocess_dataset import PreprocessDataset, read_json
from cirro.models.s3_path import S3Path
import boto3
import pandas as pd


def parse_img_files(ds: PreprocessDataset):
    """Parse the list of image files available in the workflow."""

    parent_manifest = read_json(
        ds.params["parent_dataset"] + "/web/manifest.json"
    )
    ds.logger.info("Parent dataset manifest:")
    ds.logger.info(parent_manifest)

    file_list = [
        ds.params["parent_dataset"] + "/" + file["path"]
        for file in parent_manifest["files"]
    ]
    ds.remove_param("parent_dataset")

    ds.logger.info(f"Total number of input files: {len(file_list):,}")

    # Make a list of the input files which should be copied
    img_files = [
        file_uri
        for file_uri in file_list
        if file_uri_is_image(file_uri)
    ]

    # Sort the list
    img_files.sort()

    # If there are no such files, raise an error
    msg = "No .tif files found whose names start with cyc*"
    assert len(img_files) > 0, msg

    ds.logger.info(f"Found {len(img_files):,} image files to process")
    for fn in img_files:
        ds.logger.info(fn)
    return img_files


def file_uri_is_image(file_uri: str) -> bool:
    file_name = file_uri.split("/")[-1]
    if not file_name.endswith(".tif"):
        return False
    elif not file_name.startswith('cyc'):
        return False
    else:
        return True


def stage_inputs(ds: PreprocessDataset, img_files: List[str]):
    """
    Copy the .tif files from the source dataset to params.in/raw/.
    This will mimic the intended behavior of the workflow, which will populate
    the rest of the outputs in the same "input" folder
    (even though we would consider it to be the output folder).
    """

    # Copy the files
    s3 = boto3.resource('s3')
    for file_uri in img_files:
        ds.logger.info(f"Copying {file_uri} to {ds.params['in']}")
        copy_file_s3(
            s3,
            file_uri,
            ds.params['in'] + "/raw/"
        )

    ds.logger.info("Finished copying all files")


def copy_file_s3(s3, source_file, dest_folder):
    if not dest_folder.endswith("/"):
        dest_folder = dest_folder + "/"
    dest_file = dest_folder + source_file.split("/")[-1]

    source = S3Path(source_file)
    dest = S3Path(dest_file)

    s3.meta.client.copy(
        {
            "Bucket": source.bucket,
            "Key": source.key,
        },
        dest.bucket,
        dest.key
    )


def stage_markers(ds: PreprocessDataset, img_files: List[str]):
    """Set up the markers.csv file expected by the workflow."""

    # Read the provided file
    markers_uri = ds.params["channel_names"]
    try:
        markers = pd.read_csv(
            markers_uri,
            header=None,
            names=["marker_name"]
        )
    except Exception as e:
        ds.logger.info(f"Could not read marker information from {markers_uri}")
        ds.logger.info(str(e))
        raise e

    # The number of listed markers should match the number of images
    msg = f"Number of markers ({markers.shape[0]:,}) != number of files ({len(img_files):,})"
    assert len(img_files) == markers.shape[0], msg

    # Parse the cycle number for each image
    cycle_numbers = []
    for fp in img_files:
        fn = fp.split("/")[-1]
        msg = f"Expected file names to start with cyc, not {fn}"
        assert fn.startswith("cyc"), msg

        try:
            cycle = int(fn[len("cyc"):].split("_")[0].lstrip("0"))
        except Exception as e:
            ds.logger.info(f"Could not parse cycle number for file: {fp}")
            raise e

        cycle_numbers.append(cycle)

    ds.logger.info("Found cycles:")
    ds.logger.info(cycle_numbers)

    markers = markers.assign(
        cycle=cycle_numbers,
        uri=img_files
    ).reindex(
        columns=["cycle", "marker_name", "uri"]
    )
    ds.logger.info("Formated marker table")
    ds.logger.info(markers.to_csv(index=None))

    ds.logger.info("Removing Blanks")
    markers = (
        markers
        .query("marker_name != 'Blank'")
        .query("marker_name != 'Empty'")
    )
    ds.logger.info(markers.to_csv(index=None))

    # Get the list of URIs which omit those blanks
    img_files = markers['uri'].tolist()

    # Save the marker table
    markers_uri = f"{ds.params['in']}/markers.csv"
    ds.logger.info(f"Saving to {markers_uri}")
    try:
        (
            markers
            .drop(columns=['uri'])
            .to_csv(
                markers_uri,
                index=None
            )
        )
    except Exception as e:
        ds.logger.info("Error encountered while trying to save the marker table")
        ds.logger.info(str(e))
        raise e

    # Delete the parameter used for the input file
    ds.remove_param("channel_names")

    return img_files


def setup_params(ds: PreprocessDataset):
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
    ds: PreprocessDataset,
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

    ds = PreprocessDataset.from_running()

    # Parse the list of image files which are available
    img_files = parse_img_files(ds)

    # Copy the markers file from the inputs
    # This filters the list of image files to remove any Blanks
    img_files = stage_markers(ds, img_files)

    # Copy the .tif files from the input to the $output/raw/ directory
    stage_inputs(ds, img_files)

    # Set up parameters
    setup_params(ds)

    # log
    ds.logger.info(ds.params)
