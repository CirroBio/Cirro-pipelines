from typing import Optional

import pandas as pd
from cirro.helpers.preprocess_dataset import PreprocessDataset
import json

# The genome selection nf-core/rnaseq 3 stores for each of its aligners
RNASEQ_GENOME_SELECTIONS = {
    "star_salmon": "star_salmon_genome_selection",
    "star_rsem": "star_rsem_genome_selection",
    "hisat2": "hisat2_genome_selection",
    "bowtie2_salmon": "bowtie2_salmon_genome_selection",
}


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
    set_genome(ds)
