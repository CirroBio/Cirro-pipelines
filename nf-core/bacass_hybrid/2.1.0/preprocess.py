#!/usr/bin/env python3

from cirro.helpers.preprocess_dataset import PreprocessDataset
import pandas as pd


def make_manifest(ds: PreprocessDataset) -> pd.DataFrame:

    # Filter out any index files that may have been uploaded
    ds.files = ds.files.loc[
        ds.files.apply(
            lambda r: r.get('readType', 'R') != 'I',
            axis=1
        )
    ]

    ds.logger.info("Files annotated in the input datasets:")
    ds.logger.info(ds.files.to_csv(index=None))

    ds.files = (
        ds.files
        .merge(ds.samplesheet, on="sample", how="left")
    )

    ds.logger.info("Files annotated with metadata:")
    ds.logger.info(ds.files.to_csv(index=None))

    if ds.params["merge_on"] != "sample":
        ds.logger.info(f"Merging on {ds.params['merge_on']}")
        ds.files = ds.files.rename(columns={ds.params["merge_on"]: "sample"})

    # Files that don't have a readType are assumed to be long reads
    ds.files = ds.files.assign(read=ds.files.apply(
        lambda r: "LongFastQ" if r["process"] == "long_read_fastq" else r["read"],
        axis=1
    ))

    # Make sure that there are both short and long reads present
    if "LongFastQ" not in ds.files["read"].values:
        ds.logger.info("######################################################")
        ds.logger.info("# Must provide long reads as inputs for the assembly #")
        ds.logger.info("######################################################")
        raise ValueError("No long reads detected in the input data")
    if 1 not in ds.files["read"].values:
        ds.logger.info("#######################################################")
        ds.logger.info("# Must provide short reads as inputs for the assembly #")
        ds.logger.info("#######################################################")
        raise ValueError("No short reads detected in the input data")

    # If there is more than one sample present
    if ds.files["sample"].nunique() > 1:

        # Get the number of samples for the short and long reads
        short_read_samples = ds.files.query("read == 1")["sample"].unique()
        long_read_samples = ds.files.query("read == 'LongFastQ'")["sample"].unique()

        # If there is only a single sample for both types of data
        if short_read_samples.shape[0] == 1 and long_read_samples.shape[0] == 1:

            # Just use the short read sample name for both
            ds.logger.info(f"Single sample detected for both short and long reads: {short_read_samples[0]}")
            ds.files = ds.files.replace({"sample": {long_read_samples[0]: short_read_samples[0]}})

    # Make a wide samplesheet with the columns
    # ID,R1,R2,LongFastQ,Fast5,GenomeSize
    samplesheet = (
        ds.files
        .replace({"read": {1: "R1", 2: "R2"}})
        .pivot(
            index=["sample"],
            columns="read",
            values="file"
        )
        .reset_index()
        .rename(columns=dict(sample="ID"))
        .reindex(
            columns=[
                "ID",
                "R1",
                "R2",
                "LongFastQ",
                "Fast5",
                "GenomeSize"
            ]
        )
        .fillna("NA")
    )

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
        ds.add_param("unicycler_args", f"--mode {ds.params['unicycler_mode']}")
        ds.remove_param("unicycler_mode")

    # log
    ds.logger.info(ds.params)

    # Force params.json to be written: the HealthOmics pre-process Lambda fails the run
    # when the file is absent, and the SDK writes it only when a parameter changes.
    ds.keep_params(list(ds.params.keys()))
