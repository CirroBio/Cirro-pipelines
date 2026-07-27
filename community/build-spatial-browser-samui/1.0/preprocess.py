#!/usr/bin/env python3
"""Resolve the file selected in the form into the input the workflow expects.

Every platform reader in samui is handed a folder (--folder), but the Cirro form can
only select files. The form therefore asks for one landmark file that sits at a known
depth inside the platform output, and this script turns that path back into the
folder. `spatialdata` needs a store path rather than a folder, and `files` mode passes
its four files straight through.
"""

from cirro.helpers.preprocess_dataset import PreprocessDataset

# Formats read from a folder. `auto` sniffs the folder to pick one of the others.
FOLDER_FORMATS = (
    "xenium",
    "visium",
    "visium_hd",
    "auto",
)

# Path component that sits directly inside the platform output folder, for the
# formats whose landmark file is nested below that folder. The folder is everything
# above the last occurrence of this component; for every other format it is simply
# the directory holding the selected file.
ANCHOR_DIR = {
    "visium": "spatial",         # <outs>/spatial/scalefactors_json.json
    "visium_hd": "binned_outputs",  # <outs>/binned_outputs/<bin>/spatial/scalefactors_json.json
}

# A .zarr store can be pointed at by a file at its root instead of by the store dir.
ZARR_ROOT_FILES = ("zarr.json", ".zgroup", ".zattrs", ".zmetadata")


def platform_folder(selected: str, fmt: str) -> str:
    parts = selected.rstrip("/").split("/")
    anchor = ANCHOR_DIR.get(fmt)
    if anchor is None:
        return "/".join(parts[:-1])
    if anchor not in parts:
        raise ValueError(
            f"{fmt} expects a '{anchor}/' directory in the selected path, got {selected}"
        )
    # Last occurrence: a bucket or dataset prefix may reuse the same name.
    return "/".join(parts[:len(parts) - 1 - parts[::-1].index(anchor)])


def zarr_store(selected: str) -> str:
    parts = selected.rstrip("/").split("/")
    if parts[-1] in ZARR_ROOT_FILES:
        return "/".join(parts[:-1])
    return selected


def main() -> None:
    ds = PreprocessDataset.from_running()

    fmt = ds.params["format"]
    selected = ds.params.get("input_file")
    ds.logger.info(f"format={fmt} input_file={selected}")

    if fmt == "spatialdata":
        if not selected:
            raise ValueError("Select a SpatialData store (.zarr.zip), or enter its S3 path")
        store = zarr_store(selected)
        ds.logger.info(f"SpatialData store: {store}")
        ds.add_param("zarr", store)

    elif fmt == "files":
        missing = [name for name in ("cells", "matrix") if not ds.params.get(name)]
        if missing:
            raise ValueError(f"Individual-file input requires: {', '.join(missing)}")
        for name in ("cells", "matrix", "features", "image"):
            ds.logger.info(f"{name}: {ds.params.get(name)}")

    elif fmt in FOLDER_FORMATS:
        if selected:
            folder = platform_folder(selected, fmt)
        elif fmt == "auto":
            folder = ds.params["dataset_folder"]
        else:
            raise ValueError(f"{fmt} requires a file to be selected to locate its folder")
        ds.logger.info(f"Reading {fmt} from folder: {folder}")
        ds.add_param("folder", folder)

    else:
        raise ValueError(f"Unrecognized data type: {fmt}")

    # Helper params used only to work out the input; the workflow does not take them.
    for name in ("input_file", "dataset_folder"):
        ds.remove_param(name, force=True)


if __name__ == "__main__":
    main()
