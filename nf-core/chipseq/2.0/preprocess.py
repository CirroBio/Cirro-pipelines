#!/usr/bin/env python3

from cirro.helpers.preprocess_dataset import PreprocessDataset
import pandas as pd


def make_manifest(ds: PreprocessDataset) -> pd.DataFrame:
    if 'readType' in ds.samplesheet.columns:
        file_filter = 'readType == "R"'
    else:
        file_filter = None

    manifest = ds.pivot_samplesheet(
        pivot_columns=['read'],
        column_prefix='fastq_',
        metadata_columns=['group', 'antibody', 'replicate', 'control'],
        file_filter_predicate=file_filter
    )
    assert manifest.shape[0] > 0, "No files detected -- there may be an error with data ingest"

    manifest = manifest.sort_values(by="sample")

    # The user must have specified the experimental design with group, replicate, and control
    # The antibody sample annotation is optional
    msg = "Samples must be annotated by group, replicate, and control"
    for cname in ['group', 'replicate', 'control']:
        assert cname in ds.samplesheet.columns.values, msg

    # All of the values in the `control` columns must match a value in `group`
    all_groups = ds.samplesheet["group"].dropna().drop_duplicates().tolist()
    all_controls = ds.samplesheet["control"].dropna().drop_duplicates().tolist()
    ds.logger.info(f"Groups: {', '.join(all_groups)}")
    ds.logger.info(f"Controls: {', '.join(all_controls)}")
    for control in all_controls:
        assert control in all_groups, f"Control '{control}' not found in group column"

    # All replicates must be integers
    ds.logger.info("Making sure that all replicate values are integers")
    manifest['replicate'] = manifest['replicate'].apply(lambda v: int(float(v)))
    ds.logger.debug(manifest.to_csv(index=False))

    ds.logger.info("Rearranging the columns")
    # Change group column into sample
    manifest['sample'] = manifest['group']
    # Control replicate should not be defined if it is a control row, and should default to 1 if not defined
    manifest['control_replicate'] = manifest.apply(
        lambda row: None if not row.get('control') else (row.get('control_replicate') or 1),
        axis=1
    )
    # Drop other columns
    manifest = manifest.reindex(
        columns=["sample", "fastq_1", "fastq_2", "replicate", "antibody", "control", "control_replicate"]
    )
    ds.logger.info(manifest.to_csv(index=False))

    return manifest


if __name__ == "__main__":

    ds = PreprocessDataset.from_running()
    manifest = make_manifest(ds)

    # Save the manifest
    manifest.to_csv("design.csv", index=None)

    # Add the param for the manifest
    ds.add_param("input", "design.csv")

    # log
    ds.logger.info(ds.params)
