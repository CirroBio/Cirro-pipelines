from typing import TYPE_CHECKING, Callable, Optional, Tuple

import pandas as pd

if TYPE_CHECKING:
    from cirro.helpers.preprocess_dataset import PreprocessDataset

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


def read_contrasts(path) -> pd.DataFrame:
    """Read the contrasts filled out in the form, keeping every value as text."""

    # The form field is marked required, but an empty table can still be submitted,
    # in which case Cirro either omits the param or writes an empty file
    if not path:
        raise ValueError("No contrasts were provided: add a row for each comparison")
    try:
        contrasts = pd.read_csv(path, dtype=str, keep_default_na=False)
    except pd.errors.EmptyDataError:
        contrasts = pd.DataFrame()
    if contrasts.empty:
        raise ValueError("The contrasts table is empty: add a row for each comparison")
    return contrasts


def read_sample_metadata(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t" if path.endswith(".tsv") else ",", dtype=str, keep_default_na=False)


def check_contrasts(contrasts: pd.DataFrame, sample_metadata: pd.DataFrame):
    """Every contrast must name columns and values of the sample metadata.

    The sample metadata is a free-form table picked in the form, so the contrasts
    table's schema cannot check its rows against it.
    """

    repeated = contrasts.loc[contrasts["id"].duplicated(), "id"].unique()
    if len(repeated) > 0:
        raise ValueError(f"Each contrast needs its own ID, but these are repeated: {', '.join(repeated)}")

    problems = []
    for contrast in contrasts.itertuples():
        if contrast.variable not in sample_metadata.columns:
            problems.append(f"{contrast.id}: no column {contrast.variable} in the sample metadata")
            continue
        values = set(sample_metadata[contrast.variable])
        for value in (contrast.reference, contrast.target):
            if value not in values:
                problems.append(f"{contrast.id}: no sample has {contrast.variable} = {value}")
        for column in filter(None, contrast.blocking.split(";")):
            if column not in sample_metadata.columns:
                problems.append(f"{contrast.id}: no blocking column {column} in the sample metadata")
    if problems:
        raise ValueError("The contrasts do not match the sample metadata: " + "; ".join(problems))


def find_matrices(rnaseq_results: str, exists: Callable[[str], bool]) -> Tuple[str, Optional[str]]:
    """
    The gene count matrix and, when published, the gene length matrix in the picked
    nf-core/rnaseq output. Which ones it has depends on the aligner and pseudo-aligner
    nf-core/rnaseq ran with.
    """

    # A whole dataset picked in the folder picker arrives as its data/ folder with a
    # trailing slash
    data_dir = rnaseq_results.rstrip("/")

    for prefix in RNASEQ_MATRIX_PREFIXES:
        matrix = f"{data_dir}/{prefix}.merged.gene_counts.tsv"
        if exists(matrix):
            lengths = f"{data_dir}/{prefix}.merged.gene_lengths.tsv"
            return matrix, lengths if exists(lengths) else None

    expected = ", ".join(f"{prefix}.merged.gene_counts.tsv" for prefix in RNASEQ_MATRIX_PREFIXES)
    raise ValueError(
        f"No gene count matrix found in {data_dir}/ (looked for {expected}). "
        "Pick the output of an nf-core/rnaseq run that quantified gene expression."
    )


def s3_object_exists(s3, uri: str) -> bool:
    from cirro.models.s3_path import S3Path

    path = S3Path(uri)
    listing = s3.list_objects_v2(Bucket=path.bucket, Prefix=path.key, MaxKeys=1)
    return any(obj["Key"] == path.key for obj in listing.get("Contents", []))


def find_rnaseq_dataset(inputs: list, rnaseq_results: str) -> dict:
    """The input dataset the nf-core/rnaseq output was picked from.

    The run's inputs are every dataset the form's picks point into, so the sample
    metadata may come from another one.
    """
    for dataset in inputs:
        if rnaseq_results.startswith(dataset["s3"].rstrip("/") + "/"):
            return dataset
    raise ValueError(f"None of the run's input datasets holds {rnaseq_results}")


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
        # The form keeps its default igenomes genome (GRCh38) even when the
        # run used a Cirro genome index, so only trust it for igenomes runs
        if selection.get("genome_source", "igenomes") != "igenomes":
            return None
        return selection.get("igenomes", {}).get("genome")

    return input_params.get("igenomes", {}).get("genome")


def set_matrices(ds: "PreprocessDataset"):
    import boto3

    s3 = boto3.client("s3")
    matrix, lengths = find_matrices(ds.params["rnaseq_results"], lambda uri: s3_object_exists(s3, uri))
    ds.logger.info(f"Using count matrix {matrix}")
    ds.add_param("matrix", matrix)
    if lengths is None:
        ds.logger.warning(f"No gene lengths found beside {matrix} - running without length normalization")
    else:
        ds.add_param("feature_length_matrix", lengths)


def set_genome(ds: "PreprocessDataset"):
    """
    Use the genome parameter which was selected for the input dataset.
    """

    # Get the metadata set up for this dataset, which
    # includes the params of the input datasets
    input_dataset = find_rnaseq_dataset(ds.metadata["inputs"], ds.params["rnaseq_results"])

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

    # Imported here so that test_preprocess.py runs without the SDK
    from cirro.helpers.preprocess_dataset import PreprocessDataset

    ds = PreprocessDataset.from_running()

    check_contrasts(read_contrasts(ds.params.get("contrasts")), read_sample_metadata(ds.params["input"]))
    set_matrices(ds)
    set_genome(ds)

    # Consumed above; not a parameter the workflow declares
    ds.remove_param("rnaseq_results")
