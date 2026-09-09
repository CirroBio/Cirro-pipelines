import pandas as pd
from cirro.helpers.preprocess_dataset import PreprocessDataset
import json


def make_samplesheet(ds: PreprocessDataset):
    samplesheet = ds.samplesheet
    # Drop rows where sample column containing 'samtools'
    samplesheet = samplesheet[~samplesheet["sample"].str.contains("samtools")]

    variable = ds.params["variable"]
    if variable not in samplesheet:
        raise ValueError(f"Column {variable} not found in samplesheet")

    # Write to the dataset's config/ folder (mapped in process-input.json)
    samplesheet.to_csv(ds.params["input"], index=None)
    ds.logger.info(samplesheet.to_csv(index=None))


def make_contrasts(ds: PreprocessDataset):
    # Generate the contrasts.csv
    variable = ds.params["variable"]
    reference = ds.params["reference"]
    target = ds.params["target"]
    ds.remove_param("reference")
    ds.remove_param("target")
    ds.remove_param("variable")

    contrasts = pd.DataFrame([dict(
        id=f"{reference}_vs_{target}",
        variable=variable,
        reference=reference,
        target=target
    )])

    # Write to the dataset's config/ folder (mapped in process-input.json)
    contrasts.to_csv(ds.params["contrasts"], index=None)
    ds.logger.info(contrasts.to_csv(index=None))


def set_genome(ds: PreprocessDataset):
    """
    Use the genome parameter which was selected for the input dataset.
    """

    # Get the metadata set up for this dataset, which
    # includes the params of the input dataset
    input_params = ds.metadata["inputs"][0]["params"]

    ds.add_param("genome", input_params["igenomes"]["genome"])


if __name__ == "__main__":
    ds = PreprocessDataset.from_running()
    make_samplesheet(ds)
    make_contrasts(ds)
    set_genome(ds)
