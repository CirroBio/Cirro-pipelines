#!/usr/bin/env python3
import pandas as pd
from cirro.helpers.preprocess_dataset import PreprocessDataset


def resolve_references(ds: PreprocessDataset, *params: str):
    """Make reference params absolute against the references bucket.

    Form values name a path under that bucket rather than a full URI, so the
    configuration is not tied to one deployment. A value that is already absolute
    is left alone, which keeps datasets created before that change re-runnable.
    """
    for param in params:
        value = ds.params.get(param)
        if not isinstance(value, str) or not value or value.startswith("s3://"):
            continue
        ds.add_param(param, f"{ds.references_base}/{value}", overwrite=True)


def make_samplesheet(ds: PreprocessDataset) -> pd.DataFrame:
    msg = "No files detected -- there may be a data ingest error"
    assert ds.files.shape[0] > 0, msg

    ds.logger.info("File listing for dataset:")
    ds.logger.info(ds.files.to_csv(index=None))

    for cname in ['sample', 'file']:
        assert cname in ds.files.columns.values, f"Missing column {cname}"

    return ds.files.assign(
        sample=lambda d: d['sample'].apply(
            lambda s: s.replace("/", "_")
        )
    ).reindex(
        columns=['sample', 'file']
    ).rename(
        columns=dict(
            file='vcf'
        )
    )


if __name__ == '__main__':
    ds = PreprocessDataset.from_running()

    resolve_references(ds, "data_sources")

    samplesheet = make_samplesheet(ds)
    ds.logger.info("Samplesheet:")
    ds.logger.info(samplesheet.to_csv(index=None))

    # Write to the dataset's config/ folder (mapped in process-input.json)
    samplesheet.to_csv(ds.params["samplesheet"], index=None)

    # Force params.json to be written: the HealthOmics pre-process Lambda fails the run
    # when the file is absent, and the SDK writes it only when a parameter changes.
    ds.keep_params(list(ds.params.keys()))
