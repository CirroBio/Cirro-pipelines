import csv
from cirro.helpers.preprocess_dataset import PreprocessDataset
import json


def make_samplesheet(ds: PreprocessDataset):
    samplesheet = ds.samplesheet
    # Drop rows where sample column containing 'samtools'
    samplesheet = samplesheet[~samplesheet["sample"].str.contains("samtools")]

    variable = ds.params["variable"]
    if variable not in samplesheet:
        raise ValueError(f"Column {variable} not found in samplesheet")

    # Save to a file
    samplesheet.to_csv("samplesheet.csv", index=None)

    # Set up a workflow param pointing to that file (e.g., for nf-core/rnaseq)
    ds.logger.info(samplesheet.to_csv(index=None))


def make_contrasts(ds: PreprocessDataset):
    # Generate the contrasts.csv
    variable = ds.params["variable"]
    reference = ds.params["reference"]
    target = ds.params["target"]
    ds.remove_param("reference")
    ds.remove_param("target")
    ds.remove_param("variable")

    with open("contrasts.csv", "w") as f:
        csv.writer(f).writerow(["id", "variable", "reference", "target"])
        csv.writer(f).writerow(
            [f"{reference}_vs_{target}", variable, reference, target]
        )
    with open("contrasts.csv", "r") as f:
        ds.logger.info(f.readlines())


def set_genome(ds: PreprocessDataset):
    """
    Use the genome parameter which was selected for the input dataset.
    """
    input_params = ds.metadata["inputs"][0]["params"]
    genome = input_params.get("genome")
    if genome:
        ds.add_param("genome", genome)
        ds.logger.info(f"genome={genome}")
    else:
        ds.logger.info("No iGenomes genome param found in input dataset; skipping genome annotation")


def set_matrix_paths(ds: PreprocessDataset):
    """Set the count matrix paths based on the aligner used in the input rnaseq dataset."""
    input_params = ds.metadata["inputs"][0]["params"]
    aligner = input_params.get("aligner", "star_salmon")
    ds.logger.info(f"Input dataset aligner: {aligner}")

    # Map each aligner to the subdirectory where Salmon count matrices are written
    aligner_dir_map = {
        "star_salmon": "star_salmon",
        "bowtie2_salmon": "bowtie2_salmon",
    }
    aligner_dir = aligner_dir_map.get(aligner, "salmon")

    data_path = ds.metadata["inputs"][0]["dataPath"]
    ds.logger.info(f"Input dataset dataPath: {data_path}")

    matrix = f"{data_path}/{aligner_dir}/salmon.merged.gene_counts.tsv"
    transcript_length_matrix = f"{data_path}/data/{aligner_dir}/salmon.merged.gene_lengths.tsv"

    ds.logger.info(f"matrix={matrix}")
    ds.logger.info(f"transcript_length_matrix={transcript_length_matrix}")

    ds.add_param("matrix", matrix)
    ds.add_param("transcript_length_matrix", transcript_length_matrix)


if __name__ == "__main__":
    ds = PreprocessDataset.from_running()
    make_samplesheet(ds)
    make_contrasts(ds)
    set_genome(ds)
    set_matrix_paths(ds)
