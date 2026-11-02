#!/usr/bin/env python3

import pandas as pd

# bacass 2.1.0 builds its input channel from all six columns, and tells an absent input
# apart only by the literal "NA": a missing column or a blank cell reaches the workflow
# as an empty list, which it treats as a file to process
BACASS_COLUMNS = ["ID", "R1", "R2", "LongFastQ", "Fast5", "GenomeSize"]


def read_samplesheet(path) -> pd.DataFrame:
    """Read the tab-separated samplesheet filled out in the form, keeping every value as text."""

    # The form field is marked required, but an empty table can still be submitted,
    # in which case Cirro either omits the param or writes an empty file
    if not path:
        raise ValueError("No samplesheet was provided: add at least one sample")
    try:
        samplesheet = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    except pd.errors.EmptyDataError:
        samplesheet = pd.DataFrame()
    if samplesheet.empty:
        raise ValueError("The samplesheet is empty: add at least one sample")
    return samplesheet


def check_unique_ids(samplesheet: pd.DataFrame):
    """Fail before the run on a repeated ID, which bacass's own samplesheet check rejects.

    bacass 2.1.0 cannot combine several FASTQ files of one sample.
    """

    repeated = samplesheet.loc[samplesheet["ID"].duplicated(), "ID"].unique()
    if len(repeated) > 0:
        raise ValueError(f"Each sample can have only one row, but these are repeated: {', '.join(repeated)}")


def to_bacass_samplesheet(samplesheet: pd.DataFrame) -> pd.DataFrame:
    """Add the columns the form does not ask for, and mark every empty cell NA."""

    return samplesheet.reindex(columns=BACASS_COLUMNS, fill_value="NA").replace("", "NA")


if __name__ == "__main__":

    # Imported here so that test_preprocess.py runs without the SDK
    from cirro.helpers.preprocess_dataset import PreprocessDataset

    ds = PreprocessDataset.from_running()

    samplesheet = read_samplesheet(ds.params.get("input"))
    check_unique_ids(samplesheet)
    samplesheet = to_bacass_samplesheet(samplesheet)
    ds.logger.info("Samplesheet:")
    ds.logger.info(samplesheet.to_csv(index=None, sep="\t"))

    # Overwrite the file Cirro wrote to the dataset's config/ folder (mapped in process-input.json)
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
