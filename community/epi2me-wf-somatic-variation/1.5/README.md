# Somatic Variant Calling from Nanopore Sequencing

Somatic SNV/indel, SV and modified base calling from a tumor/normal pair or a tumor sample alone
([epi2me-labs/wf-somatic-variation](https://github.com/epi2me-labs/wf-somatic-variation))

Oxford Nanopore Technologies provides an analysis workflow for calling somatic variants from
long-read nanopore sequencing of cancer samples. This analysis workflow includes:

- somatic small variant calling (ClairS for a tumor/normal pair, ClairS-TO for tumor-only)
- somatic structural variant calling (Severus)
- modified base aggregation (modkit), with differential methylation (DSS) when a normal is provided

> Please note, the workflow is designed for whole-genome data and a complete analysis of a
60X/30X tumor/normal pair takes roughly 6h30m on 64 CPUs and 256 GB of memory.

Workflow repository: [CirroBio/wf-somatic-variation](https://github.com/CirroBio/wf-somatic-variation),
a fork of [epi2me-labs/wf-somatic-variation](https://github.com/epi2me-labs/wf-somatic-variation)

## Inputs

Nanopore sequencing data should be uploaded as either Unaligned BAM, or as Nanopore Sequencing POD5.
If POD5 files are uploaded, they should be processed by the Dorado basecalling pipeline prior to
analysis for variant calling. Both the tumor and the normal BAM must have been basecalled with the
same model.

The normal BAM is optional. Without it the workflow runs a reduced tumor-only analysis: small
variants are called with ClairS-TO, structural variants are called against a panel of normals and
annotated with segmental duplication regions, and differential methylation is skipped.

## Outputs

Output files may be aggregated including information for all samples or provided per sample. Per-sample files will be prefixed with respective aliases and represented below as {{ alias }}.

| Title | File path | Description | Per sample or aggregated |
|-------|-----------|-------------|--------------------------|
| Workflow alignment statistics report | {{ alias }}.wf-somatic-variation-readQC-report.html | Report of the alignment statistics for each tumor/normal paired sample. | per-sample |
| Workflow SNV report | {{ alias }}.wf-somatic-snv-report.html | Report of the SNV for each tumor/normal paired sample. | per-sample |
| Workflow SV report | {{ alias }}.wf-somatic-sv-report.html | Report of the SV for each tumor/normal paired sample. | per-sample |
| Workflow MOD report | {{ alias }}.wf-somatic-mod-report.html | Report of the modified bases for each tumor/normal paired sample. | per-sample |
| Somatic short variant VCF | {{ alias }}.wf-somatic-snv.vcf.gz | VCF file with the somatic SNVs for the sample. | per-sample |
| Somatic short variant VCF index | {{ alias }}.wf-somatic-snv.vcf.gz.tbi | The index of the resulting VCF file with the somatic SNVs for the sample. | per-sample |
| Somatic structural variant VCF | {{ alias }}.wf-somatic-sv.vcf.gz | VCF file with the somatic SVs for the sample. In tumor-only mode the VCF file is annotated with HG38 segmental duplication regions under the SegDup heading. | per-sample |
| Somatic structural variant VCF index | {{ alias }}.wf-somatic-sv.vcf.gz.tbi | The index of the resulting VCF file with the somatic SVs for the sample. | per-sample |
| Modified bases BEDMethyl (normal) | {{ alias }}.wf-somatic-mod.normal.bedmethyl.gz | BED file with the aggregated modification counts for the normal sample. | per-sample |
| Modified bases summary (normal) | {{ alias }}.normal.mod_summary.tsv | Summary modification stats for the normal sample. | per-sample |
| Single-change BEDMethyl (normal) | {{ alias }}/mod/{{ mod_type }}/bedMethyl/{{ mod_type }}.{{ alias }}.wf-somatic-mod.normal.bedmethyl.gz | BED file with the aggregated modification counts for a single modification type (e.g. 5mc) for the normal sample. | per-sample |
| Single-change DSS input file (normal) | {{ alias }}/mod/{{ mod_type }}/DSS/{{ mod_type }}.{{ alias }}_normal.dss.tsv | Input text file for DSS for a single modification type (e.g. 5mc) for the normal sample. | per-sample |
| Modified bases BEDMethyl (tumor) | {{ alias }}.wf-somatic-mod.tumor.bedmethyl.gz | BED file with the aggregated modification counts for the tumor sample. | per-sample |
| Modified bases summary (tumor) | {{ alias }}.normal.mod_summary.tsv | Summary modification stats for the tumor sample. | per-sample |
| Single-change BEDMethyl (tumor) | {{ alias }}/mod/{{ mod_type }}/bedMethyl/{{ mod_type }}.{{ alias }}.wf-somatic-mod.tumor.bedmethyl.gz | BED file with the aggregated modification counts for a single modification type (e.g. 5mc) for the tumor sample. | per-sample |
| Single-change DSS input file (tumor) | {{ alias }}/mod/{{ mod_type }}/DSS/{{ mod_type }}.{{ alias }}_tumor.dss.tsv | Input text file for DSS for a single modification type (e.g. 5mc) for the tumor sample. | per-sample |
| Differentially modified loci (DML) per change type | {{ alias }}/mod/{{ mod_type }}/DML/{{ alias }}.{{ mod_type }}.dml.tsv | Differentially modified loci from DSS for a single modification type (e.g. 5mc). | per-sample |
| Differentially modified regions (DMR) per change type | {{ alias }}/mod/{{ mod_type }}/DMR/{{ alias }}.{{ mod_type }}.dmr.tsv | Differentially modified regions from DSS for a single modification type (e.g. 5mc). | per-sample |
| Alignment file (normal) | {{ alias }}/bam/normal/reads.bam | BAM or CRAM file with the aligned reads for the normal sample. | per-sample |
| Alignment file index (normal) | {{ alias }}/bam/normal/reads.bam.bai | The index of the resulting BAM or CRAM file with the aligned reads for the normal sample. | per-sample |
| Alignment file (tumor) | {{ alias }}/bam/tumor/reads.bam | BAM or CRAM file with the aligned reads for the tumor sample. | per-sample |
| Alignment file index (tumor) | {{ alias }}/bam/tumor/reads.bam.bai | The index of the resulting BAM or CRAM file with the aligned reads for the tumor sample. | per-sample |
| Haplotagged alignment file (normal) | {{ alias }}_normal.ht.bam | BAM or CRAM file with the haplotagged reads for the normal sample. Only produced when the `--phase_normal` option is selected. | per-sample |
| Haplotagged alignment file index (normal) | {{ alias }}_normal.ht.bam.bai | The index of the resulting BAM or CRAM file with the haplotagged reads for the normal sample. Only produced when the `--phase_normal` option is selected. | per-sample |
| Haplotagged alignment file (tumor) | {{ alias }}_tumor.ht.bam | BAM or CRAM file with the haplotagged reads for the tumor sample. | per-sample |
| Haplotagged alignment file index (tumor) | {{ alias }}_tumor.ht.bam.bai | The index of the resulting BAM or CRAM file with the haplotagged reads for the tumor sample. | per-sample |

## Workflow source

`code.uri` points at the CirroBio fork, `CirroBio/wf-somatic-variation` at `v1.5.1-cirro.2`:
upstream v1.5.1 with processes changed to emit only the files they create (staged inputs are
no longer re-declared as outputs and re-uploaded per task), and `getParams` cacheable so
`-resume` works (epi2me-labs/wf-human-variation PR #298).

Upstream keeps `lib/` as a symlink into a `wf-human-variation` git submodule, which does not
survive a GitHub zipball. The fork replaces both with a real `lib/` directory, vendored as a git
subtree from the `cirro-somatic-lib` branch of CirroBio/wf-human-variation, so the tag is
self-contained.

## Container images

HealthOmics maps only `public.ecr.aws` and `quay.io` upstreams, so `process-compute.config` rewrites
the workflow's bare Docker Hub images into `public.ecr.aws/cirrobio`. Every tag pinned in the
workflow's `nextflow.config` must exist under `public.ecr.aws/cirrobio/ontresearch/`, staged as a
single amd64 manifest — a multi-arch index resolves to the wrong architecture through the
HealthOmics image cache.
