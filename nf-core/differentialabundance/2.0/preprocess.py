from typing import Optional

import boto3
import pandas as pd
from cirro.helpers.preprocess_dataset import PreprocessDataset
from cirro.models.s3_path import S3Path
import json

# The genome selection nf-core/rnaseq 3 stores for each of its aligners
RNASEQ_GENOME_SELECTIONS = {
    "star_salmon": "star_salmon_genome_selection",
    "star_rsem": "star_rsem_genome_selection",
    "hisat2": "hisat2_genome_selection",
    "bowtie2_salmon": "bowtie2_salmon_genome_selection",
}

# Where nf-core/rnaseq publishes its merged gene matrices, relative to data/,
# in order of preference: the Salmon pseudo-aligner first, then the
# quantification that comes with each aligner. Each prefix is followed by
# .merged.gene_counts.tsv and .merged.gene_lengths.tsv
RNASEQ_MATRIX_PREFIXES = [
    "salmon/salmon",
    "star_salmon/salmon",
    "bowtie2_salmon/salmon",
    "star_rsem/rsem",
    "kallisto/kallisto",
]


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


def s3_object_exists(s3, uri: str) -> bool:
    path = S3Path(uri)
    listing = s3.list_objects_v2(Bucket=path.bucket, Prefix=path.key, MaxKeys=1)
    return any(obj["Key"] == path.key for obj in listing.get("Contents", []))


def set_matrices(ds: PreprocessDataset):
    """
    Point the pipeline at whichever gene matrices the input dataset contains,
    which depends on the aligner and pseudo-aligner nf-core/rnaseq ran with.
    """

    # process-input.json maps the matrix to data/salmon/salmon.merged.gene_counts.tsv
    # of the input dataset; strip the last two parts to get its data/ folder
    data_dir = ds.params["matrix"].rsplit("/", 2)[0]

    s3 = boto3.client("s3")
    for prefix in RNASEQ_MATRIX_PREFIXES:
        matrix = f"{data_dir}/{prefix}.merged.gene_counts.tsv"
        if not s3_object_exists(s3, matrix):
            continue

        ds.logger.info(f"Using count matrix {matrix}")
        ds.add_param("matrix", matrix, overwrite=True)

        lengths = f"{data_dir}/{prefix}.merged.gene_lengths.tsv"
        if s3_object_exists(s3, lengths):
            ds.add_param("feature_length_matrix", lengths, overwrite=True)
        else:
            ds.logger.warning(f"No gene lengths found at {lengths} - running without length normalization")
            ds.remove_param("feature_length_matrix")
        return

    expected = ", ".join(f"{prefix}.merged.gene_counts.tsv" for prefix in RNASEQ_MATRIX_PREFIXES)
    raise ValueError(
        f"No gene count matrix found in {data_dir}/ of the input dataset (looked for {expected}). "
        "The input must be an nf-core/rnaseq run that quantified gene expression."
    )


def find_igenomes_genome(input_params: dict) -> Optional[str]:
    """
    Find the iGenomes genome an input dataset was run against, if there was one.
    """

    # nf-core/rnaseq 3 keeps a separate genome selection per aligner, so the
    # branch to read depends on the aligner that was chosen
    reference_genome = (
        input_params
        .get("aligner_and_reference_genome", {})
        .get("reference_genome", {})
    )
    aligner = reference_genome.get("aligner")
    if aligner in RNASEQ_GENOME_SELECTIONS:
        selection = reference_genome.get(RNASEQ_GENOME_SELECTIONS[aligner], {})
        # No igenomes entry when the run used a Cirro genome index instead
        return selection.get("igenomes", {}).get("genome")

    return input_params.get("igenomes", {}).get("genome")


def set_genome(ds: PreprocessDataset):
    """
    Use the genome parameter which was selected for the input dataset.
    """

    # Get the metadata set up for this dataset, which
    # includes the params of the input dataset
    input_dataset = ds.metadata["inputs"][0]

    genome = find_igenomes_genome(input_dataset.get("params") or {})

    if genome is None:
        # Without a genome the pipeline annotates features from the count
        # matrix instead of a GTF, which is a usable if less detailed report
        ds.logger.warning(
            "No iGenomes reference found in the params of input dataset "
            f"{input_dataset.get('name')} ({input_dataset.get('id')}) - "
            "leaving the genome unset, so features will be annotated from the "
            "count matrix rather than a GTF"
        )
        return

    ds.add_param("genome", genome)


if __name__ == "__main__":
    ds = PreprocessDataset.from_running()
    make_samplesheet(ds)
    make_contrasts(ds)
    set_matrices(ds)
    set_genome(ds)
