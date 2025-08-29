#!/usr/bin/env python3

from cirro.helpers.preprocess_dataset import PreprocessDataset
import pandas as pd
from cirro.api.models.s3_path import S3Path
import boto3
import json


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
        manifest = ds.files.pivot(
            index=["sampleIndex", "sample", "dataset"],
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
            columns=["sampleIndex", "dataset"]
        )

    # If the files weren't added with a samplesheet.csv, then we will
    # use different logic to construct the sample sheet.
    # This is intended to capture the scenario when there are multiple
    # pairs of FASTQs for any samples.
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

    # If a subset of files have been selected
    if len(ds.params.get("subset", [])) > 0:
        ds.logger.info(f"A subset of {len(ds.params['subset']):,} files have been selected:")
        for fn in ds.params["subset"]:
            ds.logger.info(fn)
        manifest = manifest.loc[
            manifest.apply(
                lambda r: (
                    r['fastq_1'] in ds.params["subset"] or
                    r['fastq_2'] in ds.params["subset"]
                ),
                axis=1
            )
        ]
        log_table(ds, "Manifest (subset)", manifest)
        assert manifest.shape[0] > 0, "No FASTQ pairs selected in subset"
    ds.remove_param("subset", force=True)

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


def evaluate_reference_genome(ds: PreprocessDataset):

    # If the user pointed to an existing dataset
    if "compiled_dataset" in ds.params:

        # Point to the reference indexes generated from that dataset
        # Get the params for the upstream dataset
        params = read_json(
            ds.params["compiled_dataset"] + "/config/params.json"
        )

        # Set the FASTA and GTF which were used by the previous dataset
        for kw in ["fasta", "gtf"]:
            msg = f"Selected dataset missing custom genome: '{kw}' not found"
            assert kw in params, msg
            ds.add_param(kw, params[kw], overwrite=True)

        if "gencode" in params:
            ds.add_param("gencode", params["gencode"], overwrite=True)

        # Get the manifest for the upstream dataset
        ds.logger.info("Reading file list from compiled dataset")

        # List of all files in the upstream dataset
        file_list = [
            file["file"]
            for file in read_json(
                ds.params["compiled_dataset"] + "/web/manifest.json"
            ).get("files", [])
        ]

        for path, kw in [
            ("data/genome/genome.bed", "gene_bed"),
            ("data/genome/genome.transcripts.fa", "transcript_fasta"),
            ("data/genome/index/star", "star_index"),
            ("data/genome/rsem", "rsem_index"),
            ("data/genome/index/salmon", "salmon_index")
        ]:
            if any([fp.startswith(path) for fp in file_list]):
                ds.add_param(
                    kw,
                    ds.params["compiled_dataset"] + "/" + path,
                    overwrite=True
                )

        # Delete the unneeded params
        ds.remove_param("compiled_dataset")

    # If the user provided a custom genome
    elif "fasta" in ds.params or "gtf" in ds.params:

        ds.logger.info("Custom genome was selected")

        # They must have provided both a FASTA and GTF
        for kw in ["fasta", "gtf", "gencode"]:
            msg = f"Missing input for reference genome {kw.upper()}"
            assert kw in ds.params, msg

        # Save the indexed genome
        ds.add_param("save_reference", True)

    # If neither was selected, raise an error
    else:
        raise RuntimeError("Requires custom genome input")


if __name__ == "__main__":

    ds = PreprocessDataset.from_running()

    manifest = make_manifest(ds)

    # Write out the manifest
    manifest.to_csv("manifest.csv", index=None)
    ds.logger.info(f"Wrote out {manifest.shape[0]:,} lines to manifest.csv")

    evaluate_reference_genome(ds)
