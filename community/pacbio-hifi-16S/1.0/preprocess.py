#!/usr/bin/env python3

import json
from cirro.helpers.preprocess_dataset import PreprocessDataset

ds = PreprocessDataset.from_running()

kwargs = dict(index=None, sep="\t")

# Write out the list of files
input = (
    ds
    .files
    .rename(columns=dict(
        sample="sample-id",
        file="absolute-filepath"
    ))
)

ds.logger.info(input.to_csv(**kwargs))
# Write to the dataset's config/ folder (mapped in process-input.json)
input.to_csv(ds.params["input"], **kwargs)

metadata = (
    ds
    .samplesheet
    .rename(columns=dict(sample="sample_name"))
)

if "condition" not in metadata.columns.values:
    ds.logger.info("Adding dummy 'condition' metadata column")
    metadata = metadata.assign(
        condition=[
            f"cat_{i % 2}" for i in range(metadata.shape[0])
        ]
    )

ds.logger.info(metadata.to_csv(**kwargs))
metadata.to_csv(ds.params["metadata"], **kwargs)

ds.logger.info(json.dumps(ds.params, indent=4))

# Force params.json to be written: the HealthOmics pre-process Lambda fails the run
# when the file is absent, and the SDK writes it only when a parameter changes.
ds.keep_params(list(ds.params.keys()))
