#!/usr/bin/env python3

import pandas as pd
from cirro.helpers.preprocess_dataset import PreprocessDataset


ref_bucket = "s3://pubweb-references/rnafusion/3.0.2"

def make_manifest(ds: PreprocessDataset) -> pd.DataFrame:

    # Filter out any index files that may have been uploaded
    ds.files = ds.files.loc[
        ds.files.apply(
            lambda r: r.get('readType', 'R') == 'R',
            axis=1
        )
    ]

    # Make a wide manifest
    manifest = ds.wide_samplesheet(
        index=["sampleIndex", "sample", "lane"],
        columns="read",
        values="file",
        column_prefix="fastq_"
    )
    assert manifest.shape[0] > 0, "No files detected -- there may be an error with data ingest"

    # append metadata to file paths
    samples = ds.samplesheet.set_index("sample")
    manifest = manifest.assign(**{k: manifest["sample"].apply(v.get) for k, v in samples.items()})
    ordering = ['sample', 'fastq_1', 'fastq_2', 'strandedness']
    manifest = manifest.reindex(columns=ordering)

    ## check values for strandedness
    manifest = manifest.fillna({'strandedness':'unstranded'})
    strandedness = manifest.strandedness
    for i, row in strandedness.items():
        assert row in ['forward', 'reverse', 'unstranded'], print(f"ERROR: The column 'strandedness' must contain 'forward', 'reverse' or 'unstranded'.\nOffending entry:\n {manifest.iloc[i].to_frame().T}")

    return manifest


def update_params(ds, list, avail):
    for item in list:
        if item in avail:
            ds.add_param(item, True, overwrite=True)


def form_handling(ds: PreprocessDataset):
    """Process form handling"""

    params = ds.params

    # remember that genome is a placeholder pointing to pre-made s3 bucket.
    # rnafusion will never connect to iGenomes for its ref files. 
    genome = params.get('genome')
    if genome == 'GRCh38':
        genome_base = f"{ref_bucket}"
        fasta = f"{ref_bucket}/ensembl/Homo_sapiens.GRCh38.102.all.fa"
        ds.add_param('genomes_base', genome_base, overwrite=True)
        ds.add_param('fasta', fasta, overwrite=True)
        ds.remove_param('genome')


    ## parse tools list and replace with boolean params
    tools = params.get('tools')
    ds.remove_param('tools')
    avail = ['arriba', 'pizzly', 'squid', 'fusioncatcher', 'starfusion']

    update_params(ds, tools, avail)

    # workflow cannot access files within container. 
    # unpack latest github release, upl s3 and point to paths below. 
    if 'arriba' in tools:
        arriba_ref = f"{ref_bucket}/arriba"
        arriba_ref_blacklist = f"{ref_bucket}/arriba/blacklist_hg38_GRCh38_v2.4.0.tsv.gz"
        arriba_ref_protein_domain = f"{ref_bucket}/arriba/protein_domains_hg38_GRCh38_v2.4.0.gff3"
        ds.add_param('arriba_ref', arriba_ref, overwrite=True)
        ds.add_param('arriba_ref_blacklist', arriba_ref_blacklist, overwrite=True)
        ds.add_param('arriba_ref_protein_domain', arriba_ref_protein_domain, overwrite=True)

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
        ds.add_param('trimming', True, overwrite=True)


if __name__ == "__main__":

    ds = PreprocessDataset.from_running()

    manifest = make_manifest(ds)

    # Write to the dataset's config/ folder (mapped in process-input.json)
    manifest.to_csv(ds.params["input"], index=None)

    ######### Process form handling ##########
    form_handling(ds)

    # log
    ds.logger.info(ds.params)
