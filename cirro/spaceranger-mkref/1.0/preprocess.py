#!/usr/bin/env python3

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


if __name__ == "__main__":

    ds = PreprocessDataset.from_running()

    resolve_references(ds, "genome_dir", "probes")

    # process-input.json derives fasta and genes from the raw form value, which is
    # relative, so both are rebuilt here from the resolved directory. The workflow
    # takes those two paths and not the directory itself.
    genome_dir = ds.params["genome_dir"]
    ds.add_param("fasta", f"{genome_dir}/fasta/genome.fa", overwrite=True)
    ds.add_param("genes", f"{genome_dir}/genes/genes.gtf", overwrite=True)
    ds.remove_param("genome_dir")

    ds.logger.info(ds.params)
