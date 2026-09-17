#!/usr/bin/env python3
"""Prepare the launch parameters for cirrobio/spatial-data-studio.

Builds the --input map, and collapses the data-type multi-select into the
comma-separated string the workflow takes.
"""

import json
from collections import Counter

from cirro.helpers.preprocess_dataset import PreprocessDataset


def make_input_map(inputs: list) -> dict:
    """Map an output prefix to each source dataset's data folder.

    The workflow's --input takes either a data folder or a .json/.yaml file mapping an
    output prefix to each root folder. Use the map form, keyed by dataset name, so a run
    over several datasets publishes one organised output tree with each dataset's
    results under its own name.
    """
    if not inputs:
        raise ValueError("no input datasets attached to this run")

    entries = []
    for entry in inputs:
        dataset_id = entry.get("id")
        name = str(entry.get("name") or "").strip().strip("/")
        data_path = entry.get("dataPath")
        if not name:
            raise ValueError(f"input dataset {dataset_id!r} has no name")
        if not data_path:
            raise ValueError(f"input dataset {name!r} has no dataPath")
        entries.append((name, dataset_id, data_path))

    # Datasets sharing a name would silently collapse into one JSON key. Suffix every
    # instance of a repeated name with its dataset id, rather than only the ones after
    # the first, so which dataset keeps the plain name does not come down to the order
    # the inputs happen to arrive in.
    repeated = {
        name
        for name, count in Counter(name for name, _, _ in entries).items()
        if count > 1
    }

    input_map = {}
    for name, dataset_id, data_path in entries:
        prefix = f"{name}-{str(dataset_id)[:8]}" if name in repeated else name
        if prefix in input_map:
            raise ValueError(f"cannot disambiguate duplicate dataset name {prefix!r}")
        input_map[prefix] = data_path
    return input_map


def main():
    ds = PreprocessDataset.from_running()
    input_map = make_input_map(ds.metadata["inputs"])
    ds.logger.info("Input map:\n" + json.dumps(input_map, indent=4))
    with open("input-map.json", "w") as handle:
        json.dump(input_map, handle, indent=2)
    ds.add_param("input", "input-map.json", overwrite=True)

    # The form offers the catalogued data types as a multi-select; the workflow takes
    # them as one comma-separated string, and reads an empty one as every type.
    data_types = ",".join(ds.params.get("data_types") or [])
    ds.logger.info(f"Data types: {data_types or 'all'}")
    ds.add_param("data_types", data_types, overwrite=True)


if __name__ == "__main__":
    main()
