
from cirro.helpers.preprocess_dataset import PreprocessDataset
import pandas as pd
import re


def make_samplesheet(ds: PreprocessDataset):

    ds.logger.info("All Input Files:")
    ds.logger.info(ds.files.to_csv(index=None))

    # Rename the columns
    samplesheet = (
        ds.files
        .rename(
            columns=dict(file="genome", sample="name"))
        .reindex(
            columns=["genome", "name"]))
    ds.logger.info("Renamed columns:")
    ds.logger.info(samplesheet.to_csv(index=None))

    # Merge in any sample metadata
    samplesheet = pd.merge(
        samplesheet,
        ds.samplesheet,
        left_on="name",
        right_on="sample"
    ).drop(columns=["sample"])
    ds.logger.info("Added metadata:")
    ds.logger.info(samplesheet.to_csv(index=None))

    # Replace all disallowed characters
    samplesheet["name"] = samplesheet["name"].apply(
        lambda s: re.sub('[^a-zA-Z0-9_]+', '_', s)
    )

    assert samplesheet.shape[0] > 0, "No files detected"

    return samplesheet


if __name__ == "__main__":

    # Instantiate the Cirro dataset object
    ds = PreprocessDataset.from_running()

    ###############
    # SAMPLESHEET #
    ###############

    # Make the samplesheet
    samplesheet = make_samplesheet(ds)

    # Write to the dataset's config/ folder (mapped in process-input.json)
    ds.logger.info(
        f"Writing out {samplesheet.shape[0]:,} lines to {ds.params['sample_sheet']}"
    )
    samplesheet.to_csv(
        ds.params["sample_sheet"],
        index=None
    )

    # OUTPUT NAME #
    # Remove any spaces from the output name
    ds.add_param(
        "output_name",
        ds.params["output_name"].replace(" ", "_"),
        overwrite=True
    )

    # COMPARISON GROUPS #
    # If the user has specified a 'category_name'
    if len(ds.params.get("category_name", "")) > 0:

        category_name = ds.params["category_name"]

        # Make sure that the category is in the samplesheet
        assert category_name in samplesheet.columns, \
            f"Category '{category_name}' not found in samplesheet"

        # Make sure that there are at least two unique values
        assert samplesheet[category_name].nunique() > 1, \
            f"Category '{category_name}' must have at least two unique values"
