#!/usr/bin/env python3
"""Build the nf-core/sopa samplesheet for a Cirro launch.

sopa's --input is a CSV with columns sample,data_path where data_path is the
raw machine-output directory for one sample -- except for phenocycler and
ome_tif, where data_path must be the single all-channel image file itself.
Cirro selects a dataset rather than a samplesheet, so the CSV is written here
from the selected dataset's data path.
"""

from pathlib import PurePosixPath

import pandas as pd
from cirro.helpers.preprocess_dataset import PreprocessDataset

# Technologies whose data_path is one image file rather than a directory
SINGLE_FILE_TECHNOLOGIES = {
    "phenocycler": (".qptiff", ".tif", ".tiff"),
    "ome_tif": (".ome.tif", ".ome.tiff"),
}

# Stripped off an image file name to get the sample name. '.ome' is part of the
# extension of an OME-TIFF, so scan.ome.tif is sample 'scan', not 'scan.ome'.
IMAGE_EXTENSIONS = (".qptiff", ".tiff", ".tif", ".ome")


def build_samplesheet(data_path: str, technology: str, files: list[str]) -> pd.DataFrame:
    """Return the one-row samplesheet for the selected dataset.

    data_path: the input dataset's data folder (no trailing slash required)
    technology: the sopa reader in use
    files: paths of the files in the input dataset
    """
    data_path = data_path.rstrip("/")

    if technology in SINGLE_FILE_TECHNOLOGIES:
        suffixes = SINGLE_FILE_TECHNOLOGIES[technology]
        candidates = [f for f in files if f.lower().endswith(suffixes)]
        if not candidates:
            raise ValueError(
                f"technology={technology} needs a {' or '.join(suffixes)} file "
                f"in the dataset, but none was found among {len(files)} files"
            )
        if len(candidates) > 1:
            raise ValueError(
                f"technology={technology} needs exactly one image file as its "
                f"data_path, but the dataset has {len(candidates)}: "
                + ", ".join(sorted(candidates)[:5])
            )
        target = candidates[0]
        sample = PurePosixPath(target).name
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
        list(ds.files["file"]) if ds.files is not None else []
    )
    ds.logger.info("samplesheet:")
    ds.logger.info(samplesheet.to_csv(index=None))

    # Write the samplesheet to the S3 location which process-input.json maps the
    # input param to -- the dataset's own config/ folder, as nf-core/ampliseq does.
    # That records the exact input to the run alongside its results, instead of
    # leaving it in the working directory this script runs in.
    samplesheet_path = ds.params["input"]
    samplesheet.to_csv(samplesheet_path, index=None)
    ds.add_param("input", samplesheet_path, overwrite=True)

    # spatial_data exists only to carry the dataset path into this script
    ds.remove_param("spatial_data", force=True)
