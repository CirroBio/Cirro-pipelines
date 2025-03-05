# KneadData (BioBakery)

The KneadData utility is used to preprocess metagenomic data for
microbiome analysis by removing contaminating host sequences and
running quality trimming.

![](https://huttenhower.sph.harvard.edu/wp-content/uploads/elementor/thumbs/kneaddata_workflow.drawio-pffevcxp8t0xbevc5vo6hv3zc4hlapsvv6p793lm80.png)

A selection of reference database is provided so that the appropriate
host genome can be used for decontamination.

### Outputs

The files output for each samples are organized into a top-level directory
with the sample name.

#### Sequences (FASTQ)

Within that folder, there are FASTQ files for each of the scenarios:

1. Both reads in the pair pass.
  - ${SAMPLE}\_paired_[1/2].fastq.gz
2. The read in the first mate passes, and the one in the second does not pass.
  - ${SAMPLE}\_unmatched_1.fastq.gz
  - ${SAMPLE}\_${DB}\_unmatched_2.fastq.gz
3. The read in the second mate passes, and the one in the first does not pass.
  - ${SAMPLE}\_${DB}\_unmatched_1.fastq.gz
  - ${SAMPLE}\_unmatched_2.fastq.gz
4. Both reads fail
  - ${SAMPLE}\_paired_[1/2].fastq.gz

#### Quality Control (MultiQC)

The FASTQ summary of all inputs and outputs are aggregated in a single
[MultiQC report](https://seqera.io/multiqc/) - `multiqc_report.html`
