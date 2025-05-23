#!/usr/bin/env python3

from cirro.helpers.preprocess_dataset import PreprocessDataset
import pandas as pd


def make_manifest(ds: PreprocessDataset) -> pd.DataFrame:

    # Filter out any index files that may have been uploaded
    ds.files = ds.files.loc[
        ds.files.apply(
            lambda r: r.get('readType', 'R') == 'R',
            axis=1
        )
    ]

    ds.logger.info("Files annotated in the dataset:")
    ds.logger.info(ds.files.to_csv(index=None))

    # Make a wide samplesheet with the columns
    # ID,R1,R2,LongFastQ,Fast5,GenomeSize
    samplesheet = ds.files.pivot(
        index=["sampleIndex", "sample"],
        columns="read",
        values="file"
    ).rename(
        columns=lambda i: f"R{int(i)}"
    ).reset_index(
    ).set_index(
        "sampleIndex"
    ).rename(
        columns=dict(
            sample="ID"
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
    ds.logger.info(samplesheet.to_csv(index=None, sep="\t"))
    assert samplesheet.shape[0] > 0, "No files detected -- there may be an error with data ingest"

    return samplesheet


if __name__ == "__main__":

    # Load the information for this dataset
    ds = PreprocessDataset.from_running()

    # Make the samplesheet
    samplesheet = make_manifest(ds)

    # Write out to a file
    samplesheet.to_csv("samplesheet.tsv", index=None, sep="\t")

    # Add the param for the samplesheet
    ds.add_param("input", "samplesheet.tsv")

    # Format the unicycler_args based on the mode
    if ds.params.get("unicycler_mode", None) is not None:
        ds.add_param("unicycler_args", f"--mode {ds.params['unicycler_mode']}"),
        ds.remove_param("unicycler_mode")

    # log
    ds.logger.info(ds.params)
