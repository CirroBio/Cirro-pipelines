
from cirro.helpers.preprocess_dataset import PreprocessDataset


def make_samplesheet(ds: PreprocessDataset):

    ds.logger.info("All Input Files:")
    ds.logger.info(ds.files.to_csv(index=None))

    # Make a wide samplesheet
    samplesheet = (
        ds.files
        .assign(ext=ds.files["file"].apply(lambda s: s.split(".")[-1]))
        .pipe(lambda d: d.loc[d["ext"].isin(["bam", "bai"])])
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

    # Write to the dataset's config/ folder (mapped in process-input.json)
    ds.logger.info(
        f"Writing out {samplesheet.shape[0]:,} lines to {ds.params['samplesheet']}"
    )
    samplesheet.to_csv(
        ds.params["samplesheet"],
        index=None
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

    # Force params.json to be written: the HealthOmics pre-process Lambda fails the run
    # when the file is absent, and the SDK writes it only when a parameter changes.
    ds.keep_params(list(ds.params.keys()))
