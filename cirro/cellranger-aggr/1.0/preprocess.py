#!/usr/bin/env python3

from cirro.helpers.preprocess_dataset import PreprocessDataset

ds = PreprocessDataset.from_running()

# Build input as a comma-delimited list of all input dataset paths
data_paths = [dataset['dataPath'] for dataset in ds.metadata['inputs']]
assert len(data_paths) > 0, "No input datasets found"
ds.add_param("input", ",".join(data_paths))
ds.logger.info(f"input: {ds.params['input']}")

# Log the parameters present
for k, v in ds.params.items():
    ds.logger.info(f"{k}: {v}")
