#!/usr/bin/env python3

from cirro.helpers.preprocess_dataset import PreprocessDataset
import pandas as pd
from cirro.api.models.s3_path import S3Path
import boto3
import json
import urllib.request


def make_manifest(ds: PreprocessDataset) -> pd.DataFrame:
    """Construct a manifest with the paired FASTQ files in the input"""

    ds.logger.info(f"Number of files in dataset: {ds.files.shape[0]:,}")
    ds.logger.info(f"Number of samples in dataset: {ds.samplesheet.shape[0]:,}")

    # Filter out any index files that may have been uploaded
    ds.files = ds.files.loc[
        ds.files.apply(
            lambda r: r.get('readType', 'R') == 'R',
            axis=1
        )
    ]

    # If the files were added with a samplesheet.csv, they will include
    # an index indicating which line of the samplesheet.csv they were in
    if "sampleIndex" in ds.files.columns.values:

        # Reconstruct the manifest
        manifest = ds.files.reindex(
            columns=["sampleIndex", "libraryIndex", "sample", "dataset", "read", "file"]
        ).pivot(
            index=["sampleIndex", "libraryIndex", "sample", "dataset"],
            columns="read",
            values="file"
        ).rename(
            columns=lambda i: f"fastq_{int(i)}"
        ).reindex(
            columns=["fastq_1", "fastq_2"]
        ).sort_index(
        ).sort_index(
            axis=1
        ).reset_index(
        ).drop(
            columns=["sampleIndex", "libraryIndex", "dataset"]
        )

    # If the files weren't added with a samplesheet.csv, then we will
    # use different logic to construct the sample sheet.
    # This is intended to capture the scenario when there are multiple
    # pairs of FASTQs for any samples.
    # Note that this does not work with single-end data, which must
    # use the pattern above.
    else:

        manifest = []

        # Iterate over each file
        for (sample_name, _), sample_files in ds.files.groupby(["sample", "dataset"]):

            # Make a sorted list of the files available for this sample
            file_list = sample_files["file"].sort_values().tolist()

            # Add pairs of files to the manifest
            while len(file_list) > 0:

                assert len(file_list) >= 2, f"Unexpected odd number of files found for sample {sample_name}"

                # Add the first two files to the manifest
                manifest.append(
                    dict(
                        sample=sample_name,
                        fastq_1=file_list[0],
                        fastq_2=file_list[1]
                    )
                )

                # Remove those files from the list
                file_list = file_list[2:]

        manifest = pd.DataFrame(manifest)

    log_table(ds, "Manifest", manifest)

    # Get the strandedness attribute for each sample (if any exists)
    strandedness = ds.samplesheet.reindex(
        columns=["sample", "strandedness"]
    ).set_index(
        "sample"
    )["strandedness"]

    # Add that information to the manifest, filling in "unstranded" when missing
    manifest = manifest.assign(
        strandedness=manifest["sample"].apply(
            strandedness.get
        ).fillna(
            "unstranded"
        )
    )

    return manifest


def log_table(ds: PreprocessDataset, title: str, df: pd.DataFrame):
    ds.logger.info(f"{title}:")
    for line in df.to_csv(index=None).split("\n"):
        ds.logger.info(line)


def read_json(path):

    s3_path = S3Path(path)
    s3 = boto3.client('s3')
    retr = s3.get_object(Bucket=s3_path.bucket, Key=s3_path.key)
    text = retr['Body'].read().decode()
    return json.loads(text)


def filter_params_by_schema(ds: PreprocessDataset):
    """Remove any params not present in the nf-core/rnaseq nextflow_schema.json."""

    version = ds.params.get("workflow_version", "3.23.0")
    url = f"https://raw.githubusercontent.com/nf-core/rnaseq/{version}/nextflow_schema.json"

    ds.logger.info(f"Fetching nextflow_schema.json for nf-core/rnaseq {version}")
    try:
        with urllib.request.urlopen(url) as response:
            schema = json.loads(response.read().decode())
    except Exception as e:
        ds.logger.warning(f"Could not fetch nextflow_schema.json: {e}")
        return

    allowed = set()
    for section in schema.get("$defs", {}).values():
        allowed.update(section.get("properties", {}).keys())

    ds.logger.info(f"Schema defines {len(allowed):,} parameters")

    for key in list(ds.params.keys()):
        if key not in allowed:
            ds.logger.info(f"Removing param not in schema: {key}")
            ds.remove_param(key, force=True)


if __name__ == "__main__":

    ds = PreprocessDataset.from_running()

    manifest = make_manifest(ds)

    # Write out the manifest
    manifest.to_csv("manifest.csv", index=None)
    ds.logger.info(f"Wrote out {manifest.shape[0]:,} lines to manifest.csv")

    # If the user selected a Cirro dataset as the genome source, remove
    # igenomes_base and genome and resolve fasta/gtf from the index dataset
    if ds.params.get("genome_source") == "dataset":
        ds.logger.info("genome_source=dataset: removing igenomes_base and genome params")
        ds.remove_param("igenomes_base", force=True)
        ds.remove_param("genome", force=True)

        aligner = ds.params.get("aligner", "")
        if aligner in ("star_salmon", "star_rsem"):
            index_path = ds.params.get("star_index")
        elif aligner == "hisat2":
            index_path = ds.params.get("hisat2_index")
        elif aligner == "bowtie2_salmon":
            index_path = ds.params.get("bowtie2_index")
        else:
            index_path = None

        if index_path:
            index_path = index_path + "/data"
            ds.add_param("fasta", f"{index_path}/data/genome.fasta.gz")
            ds.add_param("gtf", f"{index_path}/data/genome.gtf.gz")
        else:
            ds.logger.warning(f"genome_source=dataset but no index path found for aligner={aligner!r}")

    # Remove any parameters not defined in the nf-core/rnaseq nextflow_schema.json
    filter_params_by_schema(ds)
