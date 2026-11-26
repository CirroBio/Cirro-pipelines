#!/usr/bin/env python3
"""Prepare the launch parameters for cirrobio/spatial-data-studio.

Builds the --input map from the picked folders, and collapses the data-type
multi-select into the comma-separated string the workflow takes.
"""

import json
from collections import Counter


def dataset_names(inputs: list) -> dict:
    """Each input dataset's data folder, without its trailing slash, mapped to its output name.

    Datasets sharing a name would silently collapse into one output folder. Suffix every
    instance of a repeated name with its dataset id, rather than only the ones after
    the first, so which dataset keeps the plain name does not come down to the order
    the inputs happen to arrive in.
    """
    entries = []
    for entry in inputs:
        dataset_id = entry.get("id")
        name = str(entry.get("name") or "").strip().strip("/")
        data_path = entry.get("dataPath")
        if not name:
            raise ValueError(f"input dataset {dataset_id!r} has no name")
        if not data_path:
            raise ValueError(f"input dataset {name!r} has no dataPath")
        entries.append((name, dataset_id, data_path.rstrip("/")))

    repeated = {
        name
        for name, count in Counter(name for name, _, _ in entries).items()
        if count > 1
    }
    return {
        data_path: f"{name}-{str(dataset_id)[:8]}" if name in repeated else name
        for name, dataset_id, data_path in entries
    }


def make_input_map(input_folders: str, inputs: list) -> dict:
    """Map an output prefix to each picked folder.

    The workflow's --input takes either a data folder or a .json/.yaml file mapping an
    output prefix to each root folder. Use the map form, keyed by dataset name (and the
    folder's path inside the dataset when a subfolder was picked), so a run over several
    folders publishes one organised output tree with each one's results under its own name.
    """
    if not input_folders:
        raise ValueError("No input folders were provided: pick at least one dataset or folder")

    names = dataset_names(inputs)
    input_map = {}
    for folder in dict.fromkeys(folder.rstrip("/") for folder in input_folders.split(",")):
        data_path = next(
            (path for path in names if folder == path or folder.startswith(path + "/")),
            None
        )
        if data_path is None:
            raise ValueError(f"{folder} is not inside any of the run's input datasets")
        prefix = "/".join(filter(None, [names[data_path], folder[len(data_path):].strip("/")]))
        if prefix in input_map:
            raise ValueError(f"cannot disambiguate duplicate output prefix {prefix!r}")
        input_map[prefix] = folder
    return input_map


def main():
    # Imported here so that test_preprocess.py runs without the SDK or boto3
    import boto3
    from cirro.helpers.preprocess_dataset import PreprocessDataset
    from cirro.models.s3_path import S3Path

    ds = PreprocessDataset.from_running()
    input_map = make_input_map(ds.params.get("input_folders"), ds.metadata["inputs"])
    ds.remove_param("input_folders", force=True)
    ds.logger.info("Input map:\n" + json.dumps(input_map, indent=4))

    # Written to the dataset's config/ folder (mapped in process-input.json)
    input_map_path = S3Path(ds.params["input"])
    boto3.client("s3").put_object(
        Bucket=input_map_path.bucket,
        Key=input_map_path.key,
        Body=json.dumps(input_map, indent=2).encode()
    )

    # The form offers the catalogued data types as a multi-select; the workflow takes
    # them as one comma-separated string, and reads an empty one as every type.
    data_types = ",".join(ds.params.get("data_types") or [])
    ds.logger.info(f"Data types: {data_types or 'all'}")
    ds.add_param("data_types", data_types, overwrite=True)


if __name__ == "__main__":
    main()
