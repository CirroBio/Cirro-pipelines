#!/usr/bin/env python3

from cirro.helpers.preprocess_dataset import PreprocessDataset
import pandas as pd


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
            index=["sampleIndex", "sample"],
            columns="read",
            values="file"
        ).rename(
            columns=lambda i: f"fastq_{int(i)}"
        ).sort_index(
        ).sort_index(
            axis=1
        ).reset_index(
        ).drop(
            columns=["sampleIndex"]
        )

    # If the files weren't added with a samplesheet.csv, then we will
    # use different logic to construct the sample sheet.
    # This is intended to capture the scenario when there are multiple
    # pairs of FASTQs for any samples.
    else:

        manifest = []

        # Iterate over each file
        for sample_name, sample_files in ds.files.groupby("sample"):

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

    # The user must have specified a 'condition' for each sample
    msg = "ERROR: Must specify the 'condition' for each sample (e.g. treatment, control)"
    assert 'condition' in ds.samplesheet.columns.values, msg

    # Values must be provided for all samples
    msg = "ERROR: Samples are missing 'condition' information"
    missing = ds.samplesheet.loc[
        ds.samplesheet['condition'].isnull()
    ].index.values
    if len(missing) > 0:
        ds.logger.info("Missing condition information for :")
        ds.logger.info(", ".join(missing))
    assert len(missing) == 0, msg

    # Add the condition information to the samplesheet
    manifest = manifest.assign(
        condition=manifest["sample"].apply(
            ds.samplesheet.set_index("sample")["condition"].get
        )
    )

    return manifest


def make_contrasts(manifest: pd.DataFrame):
    """Make the contrasts table expected by rnasplice"""

    # Get all of the unique values for 'condition'
    conditions = manifest["condition"].drop_duplicates().to_list()

    # Iterate over each unique combination of 'condition'
    # and make a contrasts CSV

    return pd.DataFrame([
        dict(
            contrast=f"{treatment}_{control}",
            treatment=treatment,
            control=control
        )
        for treatment in conditions
        for control in conditions
        if control < treatment
    ])


def format_reference_paths(ds: PreprocessDataset):

    refmap = {
        'GRCh37': {
            "fasta": "Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome.fa",
            "bwa": "Homo_sapiens/Ensembl/GRCh37/Sequence/BWAIndex/version0.6.0/",
            "bowtie2": "Homo_sapiens/Ensembl/GRCh37/Sequence/Bowtie2Index/",
            "star": "Homo_sapiens/Ensembl/GRCh37/Sequence/STARIndex/",
            "bismark": "Homo_sapiens/Ensembl/GRCh37/Sequence/BismarkIndex/",
            "gtf": "Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf",
            "bed12": "Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.bed",
            "readme": "Homo_sapiens/Ensembl/GRCh37/Annotation/README.txt",
            "mito_name": "MT",
        },
        'GRCh38': {
            "fasta": "Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa",
            "bwa": "Homo_sapiens/NCBI/GRCh38/Sequence/BWAIndex/version0.6.0/",
            "bowtie2": "Homo_sapiens/NCBI/GRCh38/Sequence/Bowtie2Index/",
            "star": "Homo_sapiens/NCBI/GRCh38/Sequence/STARIndex/",
            "bismark": "Homo_sapiens/NCBI/GRCh38/Sequence/BismarkIndex/",
            "gtf": "Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.gtf",
            "bed12": "Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.bed",
            "mito_name": "chrM",
        },
        'GRCm38': {
            "fasta": "Mus_musculus/Ensembl/GRCm38/Sequence/WholeGenomeFasta/genome.fa",
            "bwa": "Mus_musculus/Ensembl/GRCm38/Sequence/BWAIndex/version0.6.0/",
            "bowtie2": "Mus_musculus/Ensembl/GRCm38/Sequence/Bowtie2Index/",
            "star": "Mus_musculus/Ensembl/GRCm38/Sequence/STARIndex/",
            "bismark": "Mus_musculus/Ensembl/GRCm38/Sequence/BismarkIndex/",
            "gtf": "Mus_musculus/Ensembl/GRCm38/Annotation/Genes/genes.gtf",
            "bed12": "Mus_musculus/Ensembl/GRCm38/Annotation/Genes/genes.bed",
            "readme": "Mus_musculus/Ensembl/GRCm38/Annotation/README.txt",
            "mito_name": "MT",
        },
        'TAIR10': {
            "fasta": "Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/WholeGenomeFasta/genome.fa",
            "bwa": "Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/BWAIndex/version0.6.0/",
            "bowtie2": "Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/Bowtie2Index/",
            "star": "Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/STARIndex/",
            "bismark": "Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/BismarkIndex/",
            "gtf": "Arabidopsis_thaliana/Ensembl/TAIR10/Annotation/Genes/genes.gtf",
            "bed12": "Arabidopsis_thaliana/Ensembl/TAIR10/Annotation/Genes/genes.bed",
            "readme": "Arabidopsis_thaliana/Ensembl/TAIR10/Annotation/README.txt",
            "mito_name": "Mt"
        },
        'EB2': {
            "fasta": "Bacillus_subtilis_168/Ensembl/EB2/Sequence/WholeGenomeFasta/genome.fa",
            "bwa": "Bacillus_subtilis_168/Ensembl/EB2/Sequence/BWAIndex/version0.6.0/",
            "bowtie2": "Bacillus_subtilis_168/Ensembl/EB2/Sequence/Bowtie2Index/",
            "star": "Bacillus_subtilis_168/Ensembl/EB2/Sequence/STARIndex/",
            "bismark": "Bacillus_subtilis_168/Ensembl/EB2/Sequence/BismarkIndex/",
            "gtf": "Bacillus_subtilis_168/Ensembl/EB2/Annotation/Genes/genes.gtf",
            "bed12": "Bacillus_subtilis_168/Ensembl/EB2/Annotation/Genes/genes.bed",
            "readme": "Bacillus_subtilis_168/Ensembl/EB2/Annotation/README.txt",
        },
        'UMD3.1': {
            "fasta": "Bos_taurus/Ensembl/UMD3.1/Sequence/WholeGenomeFasta/genome.fa",
            "bwa": "Bos_taurus/Ensembl/UMD3.1/Sequence/BWAIndex/version0.6.0/",
            "bowtie2": "Bos_taurus/Ensembl/UMD3.1/Sequence/Bowtie2Index/",
            "star": "Bos_taurus/Ensembl/UMD3.1/Sequence/STARIndex/",
            "bismark": "Bos_taurus/Ensembl/UMD3.1/Sequence/BismarkIndex/",
            "gtf": "Bos_taurus/Ensembl/UMD3.1/Annotation/Genes/genes.gtf",
            "bed12": "Bos_taurus/Ensembl/UMD3.1/Annotation/Genes/genes.bed",
            "readme": "Bos_taurus/Ensembl/UMD3.1/Annotation/README.txt",
            "mito_name": "MT"
        },
        'WBcel235': {
            "fasta": "Caenorhabditis_elegans/Ensembl/WBcel235/Sequence/WholeGenomeFasta/genome.fa",
            "bwa": "Caenorhabditis_elegans/Ensembl/WBcel235/Sequence/BWAIndex/version0.6.0/",
            "bowtie2": "Caenorhabditis_elegans/Ensembl/WBcel235/Sequence/Bowtie2Index/",
            "star": "Caenorhabditis_elegans/Ensembl/WBcel235/Sequence/STARIndex/",
            "bismark": "Caenorhabditis_elegans/Ensembl/WBcel235/Sequence/BismarkIndex/",
            "gtf": "Caenorhabditis_elegans/Ensembl/WBcel235/Annotation/Genes/genes.gtf",
            "bed12": "Caenorhabditis_elegans/Ensembl/WBcel235/Annotation/Genes/genes.bed",
            "mito_name": "MtDNA"
        },
        'CanFam3.1': {
            "fasta": "Canis_familiaris/Ensembl/CanFam3.1/Sequence/WholeGenomeFasta/genome.fa",
            "bwa": "Canis_familiaris/Ensembl/CanFam3.1/Sequence/BWAIndex/version0.6.0/",
            "bowtie2": "Canis_familiaris/Ensembl/CanFam3.1/Sequence/Bowtie2Index/",
            "star": "Canis_familiaris/Ensembl/CanFam3.1/Sequence/STARIndex/",
            "bismark": "Canis_familiaris/Ensembl/CanFam3.1/Sequence/BismarkIndex/",
            "gtf": "Canis_familiaris/Ensembl/CanFam3.1/Annotation/Genes/genes.gtf",
            "bed12": "Canis_familiaris/Ensembl/CanFam3.1/Annotation/Genes/genes.bed",
            "readme": "Canis_familiaris/Ensembl/CanFam3.1/Annotation/README.txt",
            "mito_name": "MT"
        },
        'GRCz10': {
            "fasta": "Danio_rerio/Ensembl/GRCz10/Sequence/WholeGenomeFasta/genome.fa",
            "bwa": "Danio_rerio/Ensembl/GRCz10/Sequence/BWAIndex/version0.6.0/",
            "bowtie2": "Danio_rerio/Ensembl/GRCz10/Sequence/Bowtie2Index/",
            "star": "Danio_rerio/Ensembl/GRCz10/Sequence/STARIndex/",
            "bismark": "Danio_rerio/Ensembl/GRCz10/Sequence/BismarkIndex/",
            "gtf": "Danio_rerio/Ensembl/GRCz10/Annotation/Genes/genes.gtf",
            "bed12": "Danio_rerio/Ensembl/GRCz10/Annotation/Genes/genes.bed",
            "mito_name": "MT"
        },
        'BDGP6': {
            "fasta": "Drosophila_melanogaster/Ensembl/BDGP6/Sequence/WholeGenomeFasta/genome.fa",
            "bwa": "Drosophila_melanogaster/Ensembl/BDGP6/Sequence/BWAIndex/version0.6.0/",
            "bowtie2": "Drosophila_melanogaster/Ensembl/BDGP6/Sequence/Bowtie2Index/",
            "star": "Drosophila_melanogaster/Ensembl/BDGP6/Sequence/STARIndex/",
            "bismark": "Drosophila_melanogaster/Ensembl/BDGP6/Sequence/BismarkIndex/",
            "gtf": "Drosophila_melanogaster/Ensembl/BDGP6/Annotation/Genes/genes.gtf",
            "bed12": "Drosophila_melanogaster/Ensembl/BDGP6/Annotation/Genes/genes.bed",
            "mito_name": "M"
        },
        'EquCab2': {
            "fasta": "Equus_caballus/Ensembl/EquCab2/Sequence/WholeGenomeFasta/genome.fa",
            "bwa": "Equus_caballus/Ensembl/EquCab2/Sequence/BWAIndex/version0.6.0/",
            "bowtie2": "Equus_caballus/Ensembl/EquCab2/Sequence/Bowtie2Index/",
            "star": "Equus_caballus/Ensembl/EquCab2/Sequence/STARIndex/",
            "bismark": "Equus_caballus/Ensembl/EquCab2/Sequence/BismarkIndex/",
            "gtf": "Equus_caballus/Ensembl/EquCab2/Annotation/Genes/genes.gtf",
            "bed12": "Equus_caballus/Ensembl/EquCab2/Annotation/Genes/genes.bed",
            "readme": "Equus_caballus/Ensembl/EquCab2/Annotation/README.txt",
            "mito_name": "MT"
        },
        'EB1': {
            "fasta": "Escherichia_coli_K_12_DH10B/Ensembl/EB1/Sequence/WholeGenomeFasta/genome.fa",
            "bwa": "Escherichia_coli_K_12_DH10B/Ensembl/EB1/Sequence/BWAIndex/version0.6.0/",
            "bowtie2": "Escherichia_coli_K_12_DH10B/Ensembl/EB1/Sequence/Bowtie2Index/",
            "star": "Escherichia_coli_K_12_DH10B/Ensembl/EB1/Sequence/STARIndex/",
            "bismark": "Escherichia_coli_K_12_DH10B/Ensembl/EB1/Sequence/BismarkIndex/",
            "gtf": "Escherichia_coli_K_12_DH10B/Ensembl/EB1/Annotation/Genes/genes.gtf",
            "bed12": "Escherichia_coli_K_12_DH10B/Ensembl/EB1/Annotation/Genes/genes.bed",
            "readme": "Escherichia_coli_K_12_DH10B/Ensembl/EB1/Annotation/README.txt",
        },
        'Galgal4': {
            "fasta": "Gallus_gallus/Ensembl/Galgal4/Sequence/WholeGenomeFasta/genome.fa",
            "bwa": "Gallus_gallus/Ensembl/Galgal4/Sequence/BWAIndex/version0.6.0/",
            "bowtie2": "Gallus_gallus/Ensembl/Galgal4/Sequence/Bowtie2Index/",
            "star": "Gallus_gallus/Ensembl/Galgal4/Sequence/STARIndex/",
            "bismark": "Gallus_gallus/Ensembl/Galgal4/Sequence/BismarkIndex/",
            "gtf": "Gallus_gallus/Ensembl/Galgal4/Annotation/Genes/genes.gtf",
            "bed12": "Gallus_gallus/Ensembl/Galgal4/Annotation/Genes/genes.bed",
            "mito_name": "MT"
        },
        'Gm01': {
            "fasta": "Glycine_max/Ensembl/Gm01/Sequence/WholeGenomeFasta/genome.fa",
            "bwa": "Glycine_max/Ensembl/Gm01/Sequence/BWAIndex/version0.6.0/",
            "bowtie2": "Glycine_max/Ensembl/Gm01/Sequence/Bowtie2Index/",
            "star": "Glycine_max/Ensembl/Gm01/Sequence/STARIndex/",
            "bismark": "Glycine_max/Ensembl/Gm01/Sequence/BismarkIndex/",
            "gtf": "Glycine_max/Ensembl/Gm01/Annotation/Genes/genes.gtf",
            "bed12": "Glycine_max/Ensembl/Gm01/Annotation/Genes/genes.bed",
            "readme": "Glycine_max/Ensembl/Gm01/Annotation/README.txt",
        },
        'Mmul_1': {
            "fasta": "Macaca_mulatta/Ensembl/Mmul_1/Sequence/WholeGenomeFasta/genome.fa",
            "bwa": "Macaca_mulatta/Ensembl/Mmul_1/Sequence/BWAIndex/version0.6.0/",
            "bowtie2": "Macaca_mulatta/Ensembl/Mmul_1/Sequence/Bowtie2Index/",
            "star": "Macaca_mulatta/Ensembl/Mmul_1/Sequence/STARIndex/",
            "bismark": "Macaca_mulatta/Ensembl/Mmul_1/Sequence/BismarkIndex/",
            "gtf": "Macaca_mulatta/Ensembl/Mmul_1/Annotation/Genes/genes.gtf",
            "bed12": "Macaca_mulatta/Ensembl/Mmul_1/Annotation/Genes/genes.bed",
            "readme": "Macaca_mulatta/Ensembl/Mmul_1/Annotation/README.txt",
            "mito_name": "MT"
        },
        'IRGSP-1.0': {
            "fasta": "Oryza_sativa_japonica/Ensembl/IRGSP-1.0/Sequence/WholeGenomeFasta/genome.fa",
            "bwa": "Oryza_sativa_japonica/Ensembl/IRGSP-1.0/Sequence/BWAIndex/version0.6.0/",
            "bowtie2": "Oryza_sativa_japonica/Ensembl/IRGSP-1.0/Sequence/Bowtie2Index/",
            "star": "Oryza_sativa_japonica/Ensembl/IRGSP-1.0/Sequence/STARIndex/",
            "bismark": "Oryza_sativa_japonica/Ensembl/IRGSP-1.0/Sequence/BismarkIndex/",
            "gtf": "Oryza_sativa_japonica/Ensembl/IRGSP-1.0/Annotation/Genes/genes.gtf",
            "bed12": "Oryza_sativa_japonica/Ensembl/IRGSP-1.0/Annotation/Genes/genes.bed",
            "mito_name": "Mt"
        },
        'CHIMP2.1.4': {
            "fasta": "Pan_troglodytes/Ensembl/CHIMP2.1.4/Sequence/WholeGenomeFasta/genome.fa",
            "bwa": "Pan_troglodytes/Ensembl/CHIMP2.1.4/Sequence/BWAIndex/version0.6.0/",
            "bowtie2": "Pan_troglodytes/Ensembl/CHIMP2.1.4/Sequence/Bowtie2Index/",
            "star": "Pan_troglodytes/Ensembl/CHIMP2.1.4/Sequence/STARIndex/",
            "bismark": "Pan_troglodytes/Ensembl/CHIMP2.1.4/Sequence/BismarkIndex/",
            "gtf": "Pan_troglodytes/Ensembl/CHIMP2.1.4/Annotation/Genes/genes.gtf",
            "bed12": "Pan_troglodytes/Ensembl/CHIMP2.1.4/Annotation/Genes/genes.bed",
            "readme": "Pan_troglodytes/Ensembl/CHIMP2.1.4/Annotation/README.txt",
            "mito_name": "MT"
        },
        'Rnor_5.0': {
            "fasta": "Rattus_norvegicus/Ensembl/Rnor_5.0/Sequence/WholeGenomeFasta/genome.fa",
            "bwa": "Rattus_norvegicus/Ensembl/Rnor_5.0/Sequence/BWAIndex/version0.6.0/",
            "bowtie2": "Rattus_norvegicus/Ensembl/Rnor_5.0/Sequence/Bowtie2Index/",
            "star": "Rattus_norvegicus/Ensembl/Rnor_5.0/Sequence/STARIndex/",
            "bismark": "Rattus_norvegicus/Ensembl/Rnor_5.0/Sequence/BismarkIndex/",
            "gtf": "Rattus_norvegicus/Ensembl/Rnor_5.0/Annotation/Genes/genes.gtf",
            "bed12": "Rattus_norvegicus/Ensembl/Rnor_5.0/Annotation/Genes/genes.bed",
            "mito_name": "MT"
        },
        'Rnor_6.0': {
            "fasta": "Rattus_norvegicus/Ensembl/Rnor_6.0/Sequence/WholeGenomeFasta/genome.fa",
            "bwa": "Rattus_norvegicus/Ensembl/Rnor_6.0/Sequence/BWAIndex/version0.6.0/",
            "bowtie2": "Rattus_norvegicus/Ensembl/Rnor_6.0/Sequence/Bowtie2Index/",
            "star": "Rattus_norvegicus/Ensembl/Rnor_6.0/Sequence/STARIndex/",
            "bismark": "Rattus_norvegicus/Ensembl/Rnor_6.0/Sequence/BismarkIndex/",
            "gtf": "Rattus_norvegicus/Ensembl/Rnor_6.0/Annotation/Genes/genes.gtf",
            "bed12": "Rattus_norvegicus/Ensembl/Rnor_6.0/Annotation/Genes/genes.bed",
            "mito_name": "MT"
        },
        'R64-1-1': {
            "fasta": "Saccharomyces_cerevisiae/Ensembl/R64-1-1/Sequence/WholeGenomeFasta/genome.fa",
            "bwa": "Saccharomyces_cerevisiae/Ensembl/R64-1-1/Sequence/BWAIndex/version0.6.0/",
            "bowtie2": "Saccharomyces_cerevisiae/Ensembl/R64-1-1/Sequence/Bowtie2Index/",
            "star": "Saccharomyces_cerevisiae/Ensembl/R64-1-1/Sequence/STARIndex/",
            "bismark": "Saccharomyces_cerevisiae/Ensembl/R64-1-1/Sequence/BismarkIndex/",
            "gtf": "Saccharomyces_cerevisiae/Ensembl/R64-1-1/Annotation/Genes/genes.gtf",
            "bed12": "Saccharomyces_cerevisiae/Ensembl/R64-1-1/Annotation/Genes/genes.bed",
            "mito_name": "MT"
        },
        'EF2': {
            "fasta": "Schizosaccharomyces_pombe/Ensembl/EF2/Sequence/WholeGenomeFasta/genome.fa",
            "bwa": "Schizosaccharomyces_pombe/Ensembl/EF2/Sequence/BWAIndex/version0.6.0/",
            "bowtie2": "Schizosaccharomyces_pombe/Ensembl/EF2/Sequence/Bowtie2Index/",
            "star": "Schizosaccharomyces_pombe/Ensembl/EF2/Sequence/STARIndex/",
            "bismark": "Schizosaccharomyces_pombe/Ensembl/EF2/Sequence/BismarkIndex/",
            "gtf": "Schizosaccharomyces_pombe/Ensembl/EF2/Annotation/Genes/genes.gtf",
            "bed12": "Schizosaccharomyces_pombe/Ensembl/EF2/Annotation/Genes/genes.bed",
            "readme": "Schizosaccharomyces_pombe/Ensembl/EF2/Annotation/README.txt",
            "mito_name": "MT"
        },
        'Sbi1': {
            "fasta": "Sorghum_bicolor/Ensembl/Sbi1/Sequence/WholeGenomeFasta/genome.fa",
            "bwa": "Sorghum_bicolor/Ensembl/Sbi1/Sequence/BWAIndex/version0.6.0/",
            "bowtie2": "Sorghum_bicolor/Ensembl/Sbi1/Sequence/Bowtie2Index/",
            "star": "Sorghum_bicolor/Ensembl/Sbi1/Sequence/STARIndex/",
            "bismark": "Sorghum_bicolor/Ensembl/Sbi1/Sequence/BismarkIndex/",
            "gtf": "Sorghum_bicolor/Ensembl/Sbi1/Annotation/Genes/genes.gtf",
            "bed12": "Sorghum_bicolor/Ensembl/Sbi1/Annotation/Genes/genes.bed",
            "readme": "Sorghum_bicolor/Ensembl/Sbi1/Annotation/README.txt",
        },
        'Sscrofa10.2': {
            "fasta": "Sus_scrofa/Ensembl/Sscrofa10.2/Sequence/WholeGenomeFasta/genome.fa",
            "bwa": "Sus_scrofa/Ensembl/Sscrofa10.2/Sequence/BWAIndex/version0.6.0/",
            "bowtie2": "Sus_scrofa/Ensembl/Sscrofa10.2/Sequence/Bowtie2Index/",
            "star": "Sus_scrofa/Ensembl/Sscrofa10.2/Sequence/STARIndex/",
            "bismark": "Sus_scrofa/Ensembl/Sscrofa10.2/Sequence/BismarkIndex/",
            "gtf": "Sus_scrofa/Ensembl/Sscrofa10.2/Annotation/Genes/genes.gtf",
            "bed12": "Sus_scrofa/Ensembl/Sscrofa10.2/Annotation/Genes/genes.bed",
            "readme": "Sus_scrofa/Ensembl/Sscrofa10.2/Annotation/README.txt",
            "mito_name": "MT"
        },
        'AGPv3': {
            "fasta": "Zea_mays/Ensembl/AGPv3/Sequence/WholeGenomeFasta/genome.fa",
            "bwa": "Zea_mays/Ensembl/AGPv3/Sequence/BWAIndex/version0.6.0/",
            "bowtie2": "Zea_mays/Ensembl/AGPv3/Sequence/Bowtie2Index/",
            "star": "Zea_mays/Ensembl/AGPv3/Sequence/STARIndex/",
            "bismark": "Zea_mays/Ensembl/AGPv3/Sequence/BismarkIndex/",
            "gtf": "Zea_mays/Ensembl/AGPv3/Annotation/Genes/genes.gtf",
            "bed12": "Zea_mays/Ensembl/AGPv3/Annotation/Genes/genes.bed",
            "mito_name": "Mt"
        },
        'hg38': {
            "fasta": "Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa",
            "bwa": "Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/version0.6.0/",
            "bowtie2": "Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/",
            "star": "Homo_sapiens/UCSC/hg38/Sequence/STARIndex/",
            "bismark": "Homo_sapiens/UCSC/hg38/Sequence/BismarkIndex/",
            "gtf": "Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf",
            "bed12": "Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.bed",
            "mito_name": "chrM",
        },
        'hg19': {
            "fasta": "Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa",
            "bwa": "Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/version0.6.0/",
            "bowtie2": "Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/",
            "star": "Homo_sapiens/UCSC/hg19/Sequence/STARIndex/",
            "bismark": "Homo_sapiens/UCSC/hg19/Sequence/BismarkIndex/",
            "gtf": "Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf",
            "bed12": "Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.bed",
            "readme": "Homo_sapiens/UCSC/hg19/Annotation/README.txt",
            "mito_name": "chrM",
        },
        'mm10': {
            "fasta": "Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa",
            "bwa": "Mus_musculus/UCSC/mm10/Sequence/BWAIndex/version0.6.0/",
            "bowtie2": "Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/",
            "star": "Mus_musculus/UCSC/mm10/Sequence/STARIndex/",
            "bismark": "Mus_musculus/UCSC/mm10/Sequence/BismarkIndex/",
            "gtf": "Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf",
            "bed12": "Mus_musculus/UCSC/mm10/Annotation/Genes/genes.bed",
            "readme": "Mus_musculus/UCSC/mm10/Annotation/README.txt",
            "mito_name": "chrM",
        },
        'bosTau8': {
            "fasta": "Bos_taurus/UCSC/bosTau8/Sequence/WholeGenomeFasta/genome.fa",
            "bwa": "Bos_taurus/UCSC/bosTau8/Sequence/BWAIndex/version0.6.0/",
            "bowtie2": "Bos_taurus/UCSC/bosTau8/Sequence/Bowtie2Index/",
            "star": "Bos_taurus/UCSC/bosTau8/Sequence/STARIndex/",
            "bismark": "Bos_taurus/UCSC/bosTau8/Sequence/BismarkIndex/",
            "gtf": "Bos_taurus/UCSC/bosTau8/Annotation/Genes/genes.gtf",
            "bed12": "Bos_taurus/UCSC/bosTau8/Annotation/Genes/genes.bed",
            "mito_name": "chrM"
        },
        'ce10': {
            "fasta": "Caenorhabditis_elegans/UCSC/ce10/Sequence/WholeGenomeFasta/genome.fa",
            "bwa": "Caenorhabditis_elegans/UCSC/ce10/Sequence/BWAIndex/version0.6.0/",
            "bowtie2": "Caenorhabditis_elegans/UCSC/ce10/Sequence/Bowtie2Index/",
            "star": "Caenorhabditis_elegans/UCSC/ce10/Sequence/STARIndex/",
            "bismark": "Caenorhabditis_elegans/UCSC/ce10/Sequence/BismarkIndex/",
            "gtf": "Caenorhabditis_elegans/UCSC/ce10/Annotation/Genes/genes.gtf",
            "bed12": "Caenorhabditis_elegans/UCSC/ce10/Annotation/Genes/genes.bed",
            "readme": "Caenorhabditis_elegans/UCSC/ce10/Annotation/README.txt",
            "mito_name": "chrM"
        },
        'canFam3': {
            "fasta": "Canis_familiaris/UCSC/canFam3/Sequence/WholeGenomeFasta/genome.fa",
            "bwa": "Canis_familiaris/UCSC/canFam3/Sequence/BWAIndex/version0.6.0/",
            "bowtie2": "Canis_familiaris/UCSC/canFam3/Sequence/Bowtie2Index/",
            "star": "Canis_familiaris/UCSC/canFam3/Sequence/STARIndex/",
            "bismark": "Canis_familiaris/UCSC/canFam3/Sequence/BismarkIndex/",
            "gtf": "Canis_familiaris/UCSC/canFam3/Annotation/Genes/genes.gtf",
            "bed12": "Canis_familiaris/UCSC/canFam3/Annotation/Genes/genes.bed",
            "readme": "Canis_familiaris/UCSC/canFam3/Annotation/README.txt",
            "mito_name": "chrM"
        },
        'danRer10': {
            "fasta": "Danio_rerio/UCSC/danRer10/Sequence/WholeGenomeFasta/genome.fa",
            "bwa": "Danio_rerio/UCSC/danRer10/Sequence/BWAIndex/version0.6.0/",
            "bowtie2": "Danio_rerio/UCSC/danRer10/Sequence/Bowtie2Index/",
            "star": "Danio_rerio/UCSC/danRer10/Sequence/STARIndex/",
            "bismark": "Danio_rerio/UCSC/danRer10/Sequence/BismarkIndex/",
            "gtf": "Danio_rerio/UCSC/danRer10/Annotation/Genes/genes.gtf",
            "bed12": "Danio_rerio/UCSC/danRer10/Annotation/Genes/genes.bed",
            "mito_name": "chrM"
        },
        'dm6': {
            "fasta": "Drosophila_melanogaster/UCSC/dm6/Sequence/WholeGenomeFasta/genome.fa",
            "bwa": "Drosophila_melanogaster/UCSC/dm6/Sequence/BWAIndex/version0.6.0/",
            "bowtie2": "Drosophila_melanogaster/UCSC/dm6/Sequence/Bowtie2Index/",
            "star": "Drosophila_melanogaster/UCSC/dm6/Sequence/STARIndex/",
            "bismark": "Drosophila_melanogaster/UCSC/dm6/Sequence/BismarkIndex/",
            "gtf": "Drosophila_melanogaster/UCSC/dm6/Annotation/Genes/genes.gtf",
            "bed12": "Drosophila_melanogaster/UCSC/dm6/Annotation/Genes/genes.bed",
            "mito_name": "chrM"
        },
        'equCab2': {
            "fasta": "Equus_caballus/UCSC/equCab2/Sequence/WholeGenomeFasta/genome.fa",
            "bwa": "Equus_caballus/UCSC/equCab2/Sequence/BWAIndex/version0.6.0/",
            "bowtie2": "Equus_caballus/UCSC/equCab2/Sequence/Bowtie2Index/",
            "star": "Equus_caballus/UCSC/equCab2/Sequence/STARIndex/",
            "bismark": "Equus_caballus/UCSC/equCab2/Sequence/BismarkIndex/",
            "gtf": "Equus_caballus/UCSC/equCab2/Annotation/Genes/genes.gtf",
            "bed12": "Equus_caballus/UCSC/equCab2/Annotation/Genes/genes.bed",
            "readme": "Equus_caballus/UCSC/equCab2/Annotation/README.txt",
            "mito_name": "chrM"
        },
        'galGal4': {
            "fasta": "Gallus_gallus/UCSC/galGal4/Sequence/WholeGenomeFasta/genome.fa",
            "bwa": "Gallus_gallus/UCSC/galGal4/Sequence/BWAIndex/version0.6.0/",
            "bowtie2": "Gallus_gallus/UCSC/galGal4/Sequence/Bowtie2Index/",
            "star": "Gallus_gallus/UCSC/galGal4/Sequence/STARIndex/",
            "bismark": "Gallus_gallus/UCSC/galGal4/Sequence/BismarkIndex/",
            "gtf": "Gallus_gallus/UCSC/galGal4/Annotation/Genes/genes.gtf",
            "bed12": "Gallus_gallus/UCSC/galGal4/Annotation/Genes/genes.bed",
            "readme": "Gallus_gallus/UCSC/galGal4/Annotation/README.txt",
            "mito_name": "chrM"
        },
        'panTro4': {
            "fasta": "Pan_troglodytes/UCSC/panTro4/Sequence/WholeGenomeFasta/genome.fa",
            "bwa": "Pan_troglodytes/UCSC/panTro4/Sequence/BWAIndex/version0.6.0/",
            "bowtie2": "Pan_troglodytes/UCSC/panTro4/Sequence/Bowtie2Index/",
            "star": "Pan_troglodytes/UCSC/panTro4/Sequence/STARIndex/",
            "bismark": "Pan_troglodytes/UCSC/panTro4/Sequence/BismarkIndex/",
            "gtf": "Pan_troglodytes/UCSC/panTro4/Annotation/Genes/genes.gtf",
            "bed12": "Pan_troglodytes/UCSC/panTro4/Annotation/Genes/genes.bed",
            "readme": "Pan_troglodytes/UCSC/panTro4/Annotation/README.txt",
            "mito_name": "chrM"
        },
        'rn6': {
            "fasta": "Rattus_norvegicus/UCSC/rn6/Sequence/WholeGenomeFasta/genome.fa",
            "bwa": "Rattus_norvegicus/UCSC/rn6/Sequence/BWAIndex/version0.6.0/",
            "bowtie2": "Rattus_norvegicus/UCSC/rn6/Sequence/Bowtie2Index/",
            "star": "Rattus_norvegicus/UCSC/rn6/Sequence/STARIndex/",
            "bismark": "Rattus_norvegicus/UCSC/rn6/Sequence/BismarkIndex/",
            "gtf": "Rattus_norvegicus/UCSC/rn6/Annotation/Genes/genes.gtf",
            "bed12": "Rattus_norvegicus/UCSC/rn6/Annotation/Genes/genes.bed",
            "mito_name": "chrM"
        },
        'sacCer3': {
            "fasta": "Saccharomyces_cerevisiae/UCSC/sacCer3/Sequence/WholeGenomeFasta/genome.fa",
            "bwa": "Saccharomyces_cerevisiae/UCSC/sacCer3/Sequence/BWAIndex/version0.6.0/",
            "bowtie2": "Saccharomyces_cerevisiae/UCSC/sacCer3/Sequence/Bowtie2Index/",
            "star": "Saccharomyces_cerevisiae/UCSC/sacCer3/Sequence/STARIndex/",
            "bismark": "Saccharomyces_cerevisiae/UCSC/sacCer3/Sequence/BismarkIndex/",
            "readme": "Saccharomyces_cerevisiae/UCSC/sacCer3/Annotation/README.txt",
            "mito_name": "chrM"
        },
        'susScr3': {
            "fasta": "Sus_scrofa/UCSC/susScr3/Sequence/WholeGenomeFasta/genome.fa",
            "bwa": "Sus_scrofa/UCSC/susScr3/Sequence/BWAIndex/version0.6.0/",
            "bowtie2": "Sus_scrofa/UCSC/susScr3/Sequence/Bowtie2Index/",
            "star": "Sus_scrofa/UCSC/susScr3/Sequence/STARIndex/",
            "bismark": "Sus_scrofa/UCSC/susScr3/Sequence/BismarkIndex/",
            "gtf": "Sus_scrofa/UCSC/susScr3/Annotation/Genes/genes.gtf",
            "bed12": "Sus_scrofa/UCSC/susScr3/Annotation/Genes/genes.bed",
            "readme": "Sus_scrofa/UCSC/susScr3/Annotation/README.txt",
            "mito_name": "chrM"
        }
    }

    # Get the values for the genome selected by the user
    selected = refmap.get(ds.params['genome'])

    msg = f"Invalid genome selection: {ds.params['genome']}"
    assert selected is not None, msg

    base = ds.params['igenomes_base']

    for kw in ["fasta", "gtf", "star_index"]:
        val = selected[kw.replace("_index", "")]
        ds.add_param(
            kw, f"{base}/{val}"
        )


def configure_miso(ds: PreprocessDataset):

    # If the user did not provide any genes for miso
    miso_genes = ds.params.get("miso_genes")
    if miso_genes is None or len(miso_genes) == 0:

        ds.logger.info("No genes provided for MISO -- turning off")

        # Turn off the option to run miso
        ds.add_param("sashimi_plot", False, overwrite=True)


if __name__ == "__main__":

    ds = PreprocessDataset.from_running()

    manifest = make_manifest(ds)

    contrasts = make_contrasts(manifest)

    # Write the manifest and contrasts to the dataset's config/ folder
    # (mapped in process-input.json)
    manifest.to_csv(ds.params["input"], index=None)
    print(f"Wrote out {manifest.shape[0]:,} lines to {ds.params['input']}")

    contrasts.to_csv(ds.params["contrasts"], index=None)
    print(f"Wrote out {contrasts.shape[0]:,} lines to {ds.params['contrasts']}")

    # Format the reference paths
    format_reference_paths(ds)

    # Parse the form inputs from miso
    configure_miso(ds)
