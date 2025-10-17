#!/usr/bin/env python3

from cirro.helpers.preprocess_dataset import PreprocessDataset

ds = PreprocessDataset.from_running()

# Get the amount of memory selected by the user
mem = ds.params["memory"]
ds.logger.info(f"Memory allocation selected: {mem}")
assert mem.endswith(".GB"), mem

# Set the number of CPUs based on the amount of memory
cpus = int(
    max(
        1,
        int(mem[:-3]) / 8
    )
)
ds.add_param(
    "cpus",
    cpus
)
