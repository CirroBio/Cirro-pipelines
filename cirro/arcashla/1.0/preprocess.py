
from cirro.helpers.preprocess_dataset import PreprocessDataset


def make_samplesheet(ds: PreprocessDataset):

    ds.logger.info("All Input Files:")
    ds.logger.info(ds.files.to_csv(index=None))

    # Make a wide samplesheet
    samplesheet = (
        ds.files
        .assign(ext=ds.files["file"].apply(lambda s: s.split(".")[-1]))
        .pivot(
            index="sample",
            columns="ext",
            values="file"
        )
        .reindex(columns=["bam", "bai"])
        .reset_index()
    )
    ds.logger.info("Wide Samplesheet")
    ds.logger.info(samplesheet.to_csv(index=None))
    samplesheet = samplesheet.dropna()
    ds.logger.info(f"Samples with both .bam and .bam.bai: {samplesheet.shape[0]:,}")
    ds.logger.info(samplesheet.to_csv(index=None))

    assert samplesheet.shape[0] > 0, "No files detected"

    ds.logger.info("Samplesheet:")
    ds.logger.info(samplesheet.to_csv(index=None))

    return samplesheet


if __name__ == "__main__":

    # Instantiate the Cirro dataset object
    ds = PreprocessDataset.from_running()

    ###############
    # SAMPLESHEET #
    ###############

    # Make the samplesheet
    samplesheet = make_samplesheet(ds)

    # Write out the sample sheet
    ds.logger.info(
        f"Writing out {samplesheet.shape[0]:,} lines to samplesheet.csv"
    )
    samplesheet.to_csv(
        "samplesheet.csv",
        index=None
    )

    # Add it to the params
    ds.add_param(
        "samplesheet",
        "samplesheet.csv"
    )

    #########
    # GENES #
    #########
    # If 'all' was selected
    if 'all' in ds.params["genes"]:
        # Just use that
        ds.add_param("genes", "all", overwrite=True)

    # Otherwise, make a comma-separated list
    else:
        assert len(ds.params["genes"]) > 0, "Must specify at least 1 gene"
        ds.add_param("genes", ",".join(ds.params["genes"]), overwrite=True)
