#!/usr/bin/env python3

from cirro.helpers.preprocess_dataset import PreprocessDataset

ds = PreprocessDataset.from_running()

input_paths = [dataset["dataPath"] + "/" for dataset in ds.metadata["inputs"]]
assert len(input_paths) > 0, "No input datasets provided"

ds.add_param("input", ",".join(input_paths))
