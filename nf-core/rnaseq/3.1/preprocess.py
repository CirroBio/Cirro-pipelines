#!/usr/bin/env python3

from cirro.helpers.preprocess_dataset import PreprocessDataset
import pandas as pd
import json
import urllib.request


def load_samplesheet(ds: PreprocessDataset) -> pd.DataFrame:
    """Read the samplesheet filled out in the form, inferring any blank strandedness."""

    # The form field is marked required, but an empty table can still be submitted,
    # in which case Cirro either omits the param or writes an empty file
    if not ds.params.get("input"):
        raise ValueError("No samplesheet was provided: add at least one row of FASTQ files")
    try:
        # Read every value as text, so that sample names like "01" are kept as-is
        samplesheet = pd.read_csv(ds.params["input"], dtype=str, keep_default_na=False)
    except pd.errors.EmptyDataError:
        raise ValueError("The samplesheet is empty: add at least one row of FASTQ files")

    if "strandedness" not in samplesheet.columns:
        samplesheet["strandedness"] = ""
    samplesheet["strandedness"] = samplesheet["strandedness"].replace("", "auto")

    return samplesheet


def log_table(ds: PreprocessDataset, title: str, df: pd.DataFrame):
    ds.logger.info(f"{title}:")
    for line in df.to_csv(index=None).split("\n"):
        ds.logger.info(line)


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
    for section in {**schema.get("$defs", {}), **schema.get("definitions", {})}.values():
        allowed.update(section.get("properties", {}).keys())

    ds.logger.info(f"Schema defines {len(allowed):,} parameters")

    for key in list(ds.params.keys()):
        if key not in allowed:
            ds.logger.info(f"Removing param not in schema: {key}")
            ds.remove_param(key, force=True)


if __name__ == "__main__":

    ds = PreprocessDataset.from_running()

    samplesheet = load_samplesheet(ds)
    log_table(ds, "Samplesheet", samplesheet)

    # Overwrite the file Cirro wrote to the dataset's config/ folder (mapped in process-input.json)
    samplesheet.to_csv(ds.params["input"], index=None)

    # Per-aligner form params to clean up — these are never passed to the pipeline directly
    aligner_specific_params = [
        "star_salmon_genome_source", "star_salmon_genome", "star_salmon_index", "star_salmon_salmon_index",
        "star_rsem_genome_source", "star_rsem_genome", "star_rsem_star_index",
        "hisat2_genome_source", "hisat2_genome",
        "bowtie2_salmon_genome_source", "bowtie2_salmon_genome", "bowtie2_salmon_salmon_index",
    ]
    # Pipeline index params — only the active aligner's is kept
    index_params_to_remove = ["rsem_index", "hisat2_index", "bowtie2_index", "salmon_index"]

    aligner = ds.params.get("aligner", "")
    genome_source = ds.params.get(f"{aligner}_genome_source")

    # If the user selected a Cirro dataset as the genome source, remove
    # igenomes_base and resolve fasta/gtf from the index dataset
    if genome_source == "dataset":
        ds.logger.info("genome_source=dataset: removing igenomes_base")
        ds.remove_param("igenomes_base", force=True)

        if aligner == "star_salmon":
            index_path = ds.params.get("star_salmon_index")
            ds.add_param("star_index", index_path)
            salmon_index = ds.params.get("star_salmon_salmon_index")
            if salmon_index:
                # The Salmon genome-index builder nests the index under a salmon_index/
                # subfolder and publishes transcriptome.fasta.gz beside it at the root.
                ds.add_param("salmon_index", f"{salmon_index}/salmon_index")
                ds.add_param("transcript_fasta", f"{salmon_index}/transcriptome.fasta.gz")
                index_params_to_remove.remove("salmon_index")
            else:
                ds.logger.warning("star_salmon with genome_source=dataset but salmon_index not set")
        elif aligner == "star_rsem":
            index_path = ds.params.get("star_rsem_star_index")
            ds.add_param("star_index", index_path)
            rsem_index = ds.params.get("rsem_index")
            if rsem_index:
                ds.logger.info(f"star_rsem: rsem_index={rsem_index}")
                index_params_to_remove.remove("rsem_index")
            else:
                ds.logger.warning("star_rsem with genome_source=dataset but rsem_index not set")
        elif aligner == "hisat2":
            index_path = ds.params.get("hisat2_index")
            index_params_to_remove.remove("hisat2_index")
        elif aligner == "bowtie2_salmon":
            index_path = ds.params.get("bowtie2_index")
            index_params_to_remove.remove("bowtie2_index")
            salmon_index = ds.params.get("bowtie2_salmon_salmon_index")
            if salmon_index:
                # See the star_salmon branch above.
                ds.add_param("salmon_index", f"{salmon_index}/salmon_index")
                ds.add_param("transcript_fasta", f"{salmon_index}/transcriptome.fasta.gz")
                index_params_to_remove.remove("salmon_index")
            else:
                ds.logger.warning("bowtie2_salmon with genome_source=dataset but salmon_index not set")
        else:
            index_path = None

        for param in index_params_to_remove:
            ds.remove_param(param, force=True)

        if index_path:
            ds.add_param("fasta", f"{index_path}/genome.fasta")
            ds.add_param("fasta_fai", f"{index_path}/genome.fasta.fai", overwrite=True)
            ds.add_param("gtf", f"{index_path}/genome.gtf")
        else:
            ds.logger.warning(f"genome_source=dataset but no index path found for aligner={aligner!r}")

    else:
        # iGenomes: promote the aligner-specific genome selection to the pipeline genome param
        genome = ds.params.get(f"{aligner}_genome")
        if genome:
            ds.add_param("genome", genome)
            ds.logger.info(f"genome_source=igenomes: genome={genome}")
        else:
            ds.logger.warning(f"genome_source=igenomes but no genome param found for aligner={aligner!r}")

        # No pre-built indices needed for iGenomes
        for param in index_params_to_remove:
            ds.remove_param(param, force=True)

    # Clean up all per-aligner form params before passing to the pipeline
    for param in aligner_specific_params:
        ds.remove_param(param, force=True)

    # Remove any parameters not defined in the nf-core/rnaseq nextflow_schema.json
    filter_params_by_schema(ds)
