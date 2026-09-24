#!/usr/bin/env python3
"""Build the nf-core/sopa samplesheet for a Cirro launch.

sopa's --input is a CSV with columns sample,data_path where data_path is the
raw machine-output directory for one sample -- except for phenocycler and
ome_tif, where data_path must be the all-channel image file itself, which the
user picks in the form. Cirro selects a dataset rather than a samplesheet, so
the CSV is written here from the selected dataset's data path.
"""

from pathlib import PurePosixPath

import pandas as pd
from cirro.helpers.preprocess_dataset import PreprocessDataset

# Technologies whose data_path is one image file picked in the form
SINGLE_FILE_TECHNOLOGIES = ("phenocycler", "ome_tif")

# Stripped off an image file name to get the sample name. '.ome' is part of the
# extension of an OME-TIFF, so scan.ome.tif is sample 'scan', not 'scan.ome'.
IMAGE_EXTENSIONS = (".qptiff", ".tiff", ".tif", ".ome")


def build_samplesheet(
    data_path: str,
    technology: str,
    image_file: str | None = None,
) -> pd.DataFrame:
    """Return the one-row samplesheet for the selected dataset.

    data_path: the input dataset's data folder (no trailing slash required)
    technology: the sopa reader in use
    image_file: the image file picked in the form, for the single-file technologies
    """
    data_path = data_path.rstrip("/")

    if technology in SINGLE_FILE_TECHNOLOGIES:
        if not image_file:
            raise ValueError(
                f"technology={technology} needs an image file to be selected in the form"
            )
        target = image_file
        sample = PurePosixPath(image_file).name
        while (ext := PurePosixPath(sample).suffix.lower()) in IMAGE_EXTENSIONS:
            sample = sample[:-len(ext)]
        sample = sample or "sample"
    else:
        target = data_path
        sample = PurePosixPath(data_path).name or "sample"
        if sample == "data":
            # Cirro data folders are all named 'data'; use the parent (dataset id)
            sample = PurePosixPath(data_path).parent.name or "sample"

    return pd.DataFrame([dict(sample=sample, data_path=target)])


if __name__ == "__main__":

    ds = PreprocessDataset.from_running()

    samplesheet = build_samplesheet(
        ds.params["spatial_data"],
        ds.params.get("technology", "xenium"),
        ds.params.get("image_file"),
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
