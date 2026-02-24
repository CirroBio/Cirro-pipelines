#!/usr/bin/env python3

from cirro.helpers.preprocess_dataset import PreprocessDataset

ds = PreprocessDataset.from_running()

ds.logger.info("Files annotated in the dataset:")
ds.logger.info(ds.files.to_csv(index=None))

# Make a simple table with each of the files listed on single row
samplesheet = ds.files.rename(
    columns=dict(
        file="fcs"
    )
).reindex(
    columns=["sample", "fcs"]
)

# If a set of input files was selected
if 'fcs_input' in ds.params and isinstance(ds.params['fcs_input'], str) and len(ds.params['fcs_input']) > 0:

    # fcs_input must be a list
    fcs_input = ds.params['fcs_input'].split(",")

    # Make sure that each FCS file listed in `fcs_input` is valid
    assert all([fcs for fcs in fcs_input])

    # Subset the list
    ds.logger.info(f"Subsetting down to {len(fcs_input):,} FCS files")
    samplesheet = samplesheet.loc[
        samplesheet['fcs'].apply(
            lambda fcs: fcs in fcs_input
        )
    ]

    assert samplesheet.shape[0] == len(fcs_input)

samplesheet.to_csv(
    "samplesheet.csv",
    index=False
)

ds.logger.info("Samplesheet:")
ds.logger.info(samplesheet.to_csv(index=None))

# log
ds.logger.info(ds.params)