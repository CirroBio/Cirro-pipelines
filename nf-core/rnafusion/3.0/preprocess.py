#!/usr/bin/env python3

from __future__ import annotations

from typing import TYPE_CHECKING

import pandas as pd

# The SDK is imported in __main__, so test_preprocess.py runs without it
if TYPE_CHECKING:
    from cirro.helpers.preprocess_dataset import PreprocessDataset


REF_SUBDIR = "rnafusion/3.0.2"


def read_samplesheet(path) -> pd.DataFrame:
    """Read the samplesheet filled out in the form, keeping every value as text."""

    # The form field is marked required, but an empty table can still be submitted,
    # in which case Cirro either omits the param or writes an empty file
    if not path:
        raise ValueError("No samplesheet was provided: add at least one row of FASTQ files")
    try:
        # Read every value as text, so that sample names like "01" are kept as-is
        samplesheet = pd.read_csv(path, dtype=str, keep_default_na=False)
    except pd.errors.EmptyDataError:
        samplesheet = pd.DataFrame()
    if samplesheet.empty:
        raise ValueError("The samplesheet is empty: add at least one row of FASTQ files")
    return samplesheet


def fill_strandedness(samplesheet: pd.DataFrame) -> pd.DataFrame:
    """Mark every row that leaves strandedness blank as unstranded.

    The pipeline's samplesheet schema requires a value in every row.
    """

    if "strandedness" not in samplesheet.columns:
        samplesheet = samplesheet.assign(strandedness="")
    return samplesheet.assign(strandedness=samplesheet["strandedness"].replace("", "unstranded"))


def update_params(ds, list, avail):
    for item in list:
        if item in avail:
            ds.add_param(item, True, overwrite=True)


def form_handling(ds: PreprocessDataset):
    """Process form handling"""

    ref_bucket = f"{ds.references_base}/{REF_SUBDIR}"
    params = ds.params

    # remember that genome is a placeholder pointing to pre-made s3 bucket.
    # rnafusion will never connect to iGenomes for its ref files. 
    genome = params.get('genome')
    if genome == 'GRCh38':
        genome_base = f"{ref_bucket}/"
        fasta = f"{ref_bucket}/ensembl/Homo_sapiens.GRCh38.102.all.fa"
        ds.add_param('genomes_base', genome_base, overwrite=True)
        ds.add_param('fasta', fasta, overwrite=True)
        ds.remove_param('genome')

        # workflows/rnafusion.nf opens these four at the top of the file, so every
        # run stages them whichever tools were picked. Each names a directory, and
        # HealthOmics reads an S3 value without a trailing separator as an object.
        for name, subdir in [
            ('ensembl_ref', 'ensembl'),
            ('starindex_ref', 'star'),
            ('starfusion_ref', 'starfusion/ctat_genome_lib_build_dir'),
            ('fusionreport_ref', 'fusion_report_db'),
        ]:
            ds.add_param(name, f"{ref_bucket}/{subdir}/", overwrite=True)


    ## parse tools list and replace with boolean params
    tools = params.get('tools')
    ds.remove_param('tools')
    avail = ['arriba', 'fusioncatcher', 'starfusion']

    update_params(ds, tools, avail)

    # workflow cannot access files within container. 
    # unpack latest github release, upl s3 and point to paths below. 
    if 'fusioncatcher' in tools:
        ds.add_param(
            'fusioncatcher_ref', f"{ref_bucket}/fusioncatcher/human_v102/", overwrite=True
        )

    if 'arriba' in tools:
        # arriba_ref names a directory and no process reads it -- nextflow.config
        # derives it from genomes_base and only the two files below are consumed.
        # Passing it would hand HealthOmics an S3 folder with no trailing separator.
        arriba_ref_blacklist = f"{ref_bucket}/arriba/blacklist_hg38_GRCh38_v2.4.0.tsv.gz"
        arriba_ref_protein_domains = f"{ref_bucket}/arriba/protein_domains_hg38_GRCh38_v2.4.0.gff3"
        ds.add_param('arriba_ref_blacklist', arriba_ref_blacklist, overwrite=True)
        ds.add_param('arriba_ref_protein_domains', arriba_ref_protein_domains, overwrite=True)

    # separate handling for report generation, boolean param. 
    report_generation = params.get('report')
    ds.remove_param('report')
    if report_generation:
        ds.add_param('fusionreport', True, overwrite=True)
    else:
        # FUSION_REPORT/INSPECTOR are toggled by the skip_viz param in modules.config. 
        ds.add_param('skip_vis', True, overwrite=True)

    # Trimming
    trim = params.get('trim')
    ds.remove_param('trim')
    if trim:
        ds.add_param('fastp_trim', True, overwrite=True)


if __name__ == "__main__":

    from cirro.helpers.preprocess_dataset import PreprocessDataset

    ds = PreprocessDataset.from_running()

    samplesheet = fill_strandedness(read_samplesheet(ds.params.get("input")))
    ds.logger.info("Samplesheet:")
    ds.logger.info(samplesheet.to_csv(index=None))

    # Overwrite the file Cirro wrote to the dataset's config/ folder (mapped in process-input.json)
    samplesheet.to_csv(ds.params["input"], index=None)

    ######### Process form handling ##########
    form_handling(ds)

    # log
    ds.logger.info(ds.params)
