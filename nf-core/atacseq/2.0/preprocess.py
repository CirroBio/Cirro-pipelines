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

    # Make a wide manifest
    manifest = ds.wide_samplesheet(
        index=["sampleIndex", "sample", "lane"],
        columns="read",
        values="file",
        column_prefix="fastq_"
    ).sort_values(
        by="sample"
    )
    assert manifest.shape[0] > 0, "No files detected -- there may be an error with data ingest"

    # Add an arbitrary `replicate` label for each
    ds.logger.info("Adding replicate column")
    sample_vc = manifest["sample"].value_counts()
    manifest = manifest.assign(
        replicate=[
            i + 1
            for _, n_replicates in sample_vc.items()
            for i in range(n_replicates)
        ]
    )
    ds.logger.info(manifest.to_csv(index=None))

    # Only keep the expected columns, setting sample -> group
    ds.logger.info("Rearranging the columns")
    manifest = manifest.reindex(
        columns=["sample", "fastq_1", "fastq_2", "replicate"]
    )
    ds.logger.info(manifest.to_csv(index=None))

    return manifest


if __name__ == "__main__":

    ds = PreprocessDataset.from_running()
    manifest = make_manifest(ds)

    # Write to the dataset's config/ folder (mapped in process-input.json)
    manifest.to_csv(ds.params["input"], index=None)

    # log
    ds.logger.info(ds.params)
