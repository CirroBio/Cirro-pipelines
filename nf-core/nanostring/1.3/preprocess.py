#!/usr/bin/env python3

from cirro.helpers.preprocess_dataset import PreprocessDataset

ds = PreprocessDataset.from_running()

manifest = (
    ds.files
    .merge(
        ds.samplesheet,
        on="sample",
    )
    .rename(
        columns=dict(
            file="RCC_FILE",
            sample="SAMPLE_ID"
        )
    )
)

# Drop any columns which aren't being used for plotting
manifest = manifest.reindex(
    columns=[
        cname for cname in manifest.columns
        if cname in ["RCC_FILE", "SAMPLE_ID", ds.params["heatmap_id_column"]]
    ]
)

ds.logger.info("Analysis Manifest:")
for line in manifest.to_csv(index=False).split("\n"):
    ds.logger.info(line)

manifest.to_csv("manifest.csv", index=False)
