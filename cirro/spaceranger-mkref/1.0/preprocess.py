#!/usr/bin/env python3


def absolute_reference(value: str, references_base: str) -> str:
    """A reference param's value, made absolute against the references bucket.

    Form values name a path under that bucket rather than a full URI, so the
    configuration is not tied to one deployment. A value that is already absolute
    is returned as is, which keeps datasets created before that change re-runnable.
    """
    if value.startswith("s3://"):
        return value
    return f"{references_base}/{value}"


def genome_files(genome_dir: str) -> dict:
    """The FASTA and GTF inside a Cell Ranger reference folder, which the workflow reads."""
    genome_dir = genome_dir.rstrip("/")
    return {"fasta": f"{genome_dir}/fasta/genome.fa", "genes": f"{genome_dir}/genes/genes.gtf"}


def custom_globs(folders) -> str:
    """A glob for every file in each picked folder, comma-separated as the workflow reads `custom`.

    A whole dataset picked in the folder picker arrives as its data/ folder with a
    trailing slash, which is dropped so the glob has no `//`.
    """
    if not folders:
        raise ValueError("No folder of custom sequences was selected")
    return ",".join(f"{folder.rstrip('/')}/**" for folder in folders.split(","))


if __name__ == "__main__":

    # Imported here so that test_preprocess.py runs without the SDK
    from cirro.helpers.preprocess_dataset import PreprocessDataset

    ds = PreprocessDataset.from_running()

    ds.add_param("custom", custom_globs(ds.params.get("custom")), overwrite=True)
    ds.add_param("probes", absolute_reference(ds.params["probes"], ds.references_base), overwrite=True)

    # genome_dir is mapped only so this script can see it
    genome_dir = absolute_reference(ds.params["genome_dir"], ds.references_base)
    for param, value in genome_files(genome_dir).items():
        ds.add_param(param, value, overwrite=True)
    ds.remove_param("genome_dir")

    ds.logger.info(ds.params)
