#!/usr/bin/env python3
import pandas as pd
from cirro.helpers.preprocess_dataset import PreprocessDataset

def make_manifest(ds: PreprocessDataset) -> pd.DataFrame:

    # Filter out any index files that may have been uploaded
    ds.files = ds.files.loc[
        ds.files.apply(
            lambda r: r.get('readType', 'R') == 'R',
            axis=1
        )
    ]

    # Make a wide manifest
    manifest = ds.files.pivot(
        index=["sampleIndex", "sample"],
        columns="read",
        values="file"
    ).rename(
        columns=lambda i: f"fastq_{i}"
    ).reset_index(
    ).set_index(
        "sampleIndex"
    )
    ds.logger.info("Manifest")
    ds.logger.info(manifest)
    assert manifest.shape[0] > 0, "No files detected -- there may be an error with data ingest"

    # Get the sample-level metadata
    samples = ds.samplesheet.set_index("sample")
    ds.logger.info("Sample metadata")
    ds.logger.info(samples)

    # The `grouping` and `replicate` annotation must have been assigned
    for cname in ['grouping', 'replicate']:
        msg = f"The {cname} must be assigned for each sample"
        assert cname in samples.columns.values, msg

    # If control normalization is being used
    ds.logger.info("Parameter - use_control:")
    ds.logger.info(ds.params["use_control"])
    if ds.params["use_control"]:

        # The 'control' annotation must have been assigned
        assert 'control' in samples.columns.values, "The 'control' must be assigned for each sample"

    # Annotate the manifest with the sample-level metadata
    ds.logger.info("Adding group, replicate, and control to manifest")
    manifest = manifest.assign(
        group=manifest["sample"].apply(samples["grouping"].get),
        replicate=manifest["sample"].apply(samples["replicate"].get),
        control=manifest["sample"].apply(samples.reindex(columns=["control"])["control"].fillna("").get)
    )

    # There must be a group value for each sample
    ds.logger.info("Checking for samples lacking 'grouping'")
    unannotated_samples = [
        sample
        for sample, group in manifest["group"].items()
        if pd.isnull(group) or len(group) == 0
    ]
    msg = f"Samples are missing a 'grouping': {', '.join(unannotated_samples)}"
    assert len(unannotated_samples) == 0, msg

    # Each of the control values must either be blank, or correspond to a group
    ds.logger.info("Checking for controls which do not match a grouping")
    unmatched_controls = [
        f"{sample}: {control}"
        for sample, control in manifest["control"].items()
        if len(control) > 0 and control not in manifest["group"].values
    ]
    msg = f"Found 'control' labels which do not correspond to a 'grouping': {', '.join(unmatched_controls)}"
    assert len(unmatched_controls) == 0, msg

    # Arrange the table to have the expected order
    manifest = manifest.reindex(
        columns=[
            "group",
            "replicate",
            "fastq_1",
            "fastq_2",
            "control"
        ]
    )
    ds.logger.info("Reordered manifest")
    ds.logger.info(manifest)

    return manifest


if __name__ == "__main__":

    # Load the information for this dataset
    ds = PreprocessDataset.from_running()
    manifest = make_manifest(ds)

    # Write to the dataset's config/ folder (mapped in process-input.json)
    manifest.to_csv(ds.params["input"], index=None)

    # log
    ds.logger.info(ds.params)
