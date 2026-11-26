#!/usr/bin/env python3

BUNDLE_FILE = "/analysis.zarr.zip"


def bundle_folder(input_uri: str) -> str:
    """The Xenium output bundle that holds the picked analysis.zarr.zip file."""
    if not input_uri:
        raise ValueError("No input was provided: select the analysis.zarr.zip file of a Xenium output bundle")
    if not input_uri.endswith(BUNDLE_FILE):
        raise ValueError(f"Expected the input to end with '{BUNDLE_FILE}', not {input_uri}")
    return input_uri[:-len(BUNDLE_FILE)]


if __name__ == "__main__":

    # Imported here so that test_preprocess.py runs without the SDK
    from cirro.helpers.preprocess_dataset import PreprocessDataset

    ds = PreprocessDataset.from_running()

    ds.add_param("input", bundle_folder(ds.params.get("input")), overwrite=True)
