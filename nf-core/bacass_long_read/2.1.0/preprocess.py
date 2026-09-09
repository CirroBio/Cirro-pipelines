#!/usr/bin/env python3

from cirro.helpers.preprocess_dataset import PreprocessDataset
import pandas as pd


def make_manifest(ds: PreprocessDataset) -> pd.DataFrame:

    ds.logger.info("Files annotated in the dataset:")
    ds.logger.info(ds.files.to_csv(index=None))

    # Make a wide samplesheet with the columns
    # ID,R1,R2,LongFastQ,Fast5,GenomeSize
    samplesheet = ds.files.rename(
        columns=dict(
            sample="ID",
            file="LongFastQ"
        )
    ).reindex(
        columns=[
            "ID",
            "R1",
            "R2",
            "LongFastQ",
            "Fast5",
            "GenomeSize"
        ]
    ).fillna("NA")

    ds.logger.info("Formatted samplesheet:")
    ds.logger.info(samplesheet.to_csv(index=None))
    assert samplesheet.shape[0] > 0, "No files detected -- there may be an error with data ingest"

    return samplesheet


if __name__ == "__main__":

    # Load the information for this dataset
    ds = PreprocessDataset.from_running()

    # Make the samplesheet
    samplesheet = make_manifest(ds)

    # Write to the dataset's config/ folder (mapped in process-input.json)
    samplesheet.to_csv(ds.params["input"], index=None, sep="\t")

    # Format the unicycler_args based on the mode
    if ds.params.get("unicycler_mode", None) is not None:
        ds.add_param("unicycler_args", f"--mode {ds.params['unicycler_mode']}"),
        ds.remove_param("unicycler_mode")

    # log
    ds.logger.info(ds.params)
