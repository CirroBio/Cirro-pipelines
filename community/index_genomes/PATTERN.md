# index_genomes — Shared Pattern

This folder holds one process per indexing tool, each of which builds a
reference index (bulk aligner or single-cell reference) from a FASTA. Every
method follows the same shape, and every method has a matching **ingest**
process under `intake/` that lets users upload a pre-built index of the same
type instead of building it here.

## Two ways to get the same data type

For every tool there are two processes that produce the same logical data type
(an "<Tool> Index" / "<Tool> Reference"):

| Tool | Build process (this folder) | Ingest process (`intake/`) |
|------|-----------------------------|----------------------------|
| bismark | `process-cirro-genome-index-bismark-1-0` | `genome_bismark_index` |
| bowtie2 | `process-cirro-genome-index-bowtie2-1-0` | `genome_bowtie2_index` |
| bwa | `process-cirro-genome-index-bwa-1-0` | `genome_bwa_index` |
| bwa-mem2 | `process-cirro-genome-index-bwa-mem2-1-0` | `genome_bwa-mem2_index` |
| hisat2 | `process-cirro-genome-index-hisat2-1-0` | `genome_hisat2_index` |
| kallisto | `process-cirro-genome-index-kallisto-1-0` | `genome_kallisto_index` |
| rsem | `process-cirro-genome-index-rsem-1-0` | `genome_rsem_index` |
| salmon | `process-cirro-genome-index-salmon-1-0` | `genome_salmon_index` |
| star | `process-cirro-genome-index-star-1-0` | `genome_star_index` |
| simpleaf (alevin-fry) | `process-cirro-genome-index-simpleaf-1-0` | `genome_simpleaf_index` |
| kb (kallisto\|bustools) | `process-cirro-genome-index-kb-1-0` | `genome_kb_index` |
| cellranger | `process-cirro-genome-index-cellranger-1-0` | `genome_cellranger_index` |
| cellranger-arc | `process-cirro-genome-index-cellranger-arc-1-0` | `genome_cellranger-arc_index` |
| cellranger-vdj | `process-cirro-genome-index-cellranger-vdj-1-0` | `genome_cellranger-vdj_index` |

The last five are single-cell reference builders ported from nf-core/scrnaseq.
The three Cell Ranger tools pin the container by selecting a version from an
enum on a Cumulus base image; simpleaf and kb take a full image string instead
(see the container conventions below).

- **Build** (`executor: NEXTFLOW`): consumes a `genome_fasta` dataset, runs the
  index-builder, and writes an index dataset. Users start from a FASTA already
  in Cirro.
- **Ingest** (`executor: INGEST`): recognizes a directory of pre-built index
  files by filename patterns (`fileMappingRules`) and registers it as a dataset.
  Users bring an index built elsewhere.

Both paths yield a dataset that downstream alignment pipelines can consume. The
build process is essentially the compute route to the same data type the ingest
process accepts by upload.

## Anatomy of a build method (`<tool>/1.0/`)

Each method is a versioned folder with the four standard process files:

- **process-definition.json** — registers the NEXTFLOW process. The parts that
  vary per tool are `id`, `dataType` (`"Genome Index (<Tool>)"`), `name`,
  `desc`, and `code.script`. `code.script` must be the exact entrypoint
  filename in `CirroBio/nf-index-genome`, which does not always match the folder
  name — e.g. the `simpleaf` folder runs `main_alevinfry.nf`, `cellranger-arc`
  runs `main_cellrangerarc.nf`, and `cellranger-vdj` runs `main_cellranger_vdj.nf`.
  The parts that are identical across all methods:
  - `parentProcessIds: ["genome_fasta"]` — input is always a genome FASTA
  - `childProcessIds: []` — downstream consumers are not wired here
  - `category: "Reference Data"`, `executor: "NEXTFLOW"`
  - `code.uri: "CirroBio/nf-index-genome"` — one shared repo; only the entry
    script differs per tool
  - `computeDefaults`, `preProcessScript` → the shared `process-compute.config`
    and `preprocess.py` at the top of this folder (see below)

- **process-form.json** — the input form. Every method has:
  - `fasta` (required) — the reference FASTA, `pathType: "dataset"`
  - a `<tool>_extra_args` string for pass-through CLI flags, usually with a
    `ui:placeholder` example (Cell Ranger tools name these per subcommand, e.g.
    `cellranger_mkref_args`, `cellranger_mkgtf_args`, `cellranger_mkvdjref_args`)
  - a way to pin the container (see the three conventions below)
  - Most also accept a `gtf` annotation (required for the single-cell tools).
    Tool-specific extras exist (e.g. bismark's `aligner`, salmon's
    `extra_fasta` / `transcriptome_source`, kb's `kb_workflow`, Cell Ranger's
    `*_reference_name`, cellranger-arc's `cellrangerarc_motifs`).

- **process-input.json** — maps form params to the workflow. Every method maps
  `fasta`, `outdir → $.dataset.dataPath`, and its extra-args field, and pins the
  container in one of three ways. The first two splice a `<tool>_version` field
  onto a fixed base image via the pipe syntax, so the tag is interpolated at
  runtime; the third takes a whole image string:
  - **Version-tag, free-text (bulk aligners)** — a required `<tool>_version`
    free-text field (default = a biocontainers tag; note bwa-mem2 uses
    `bwamem2_version`) on the biocontainers base image:
    ```
    "container": "quay.io/biocontainers/<tool>:|$.dataset.params.<tool>_version"
    ```
  - **Version-tag, enum (Cell Ranger tools)** — a required `<tool>_version`
    enum (a curated dropdown of released versions) on the Cumulus base image,
    e.g.:
    ```
    "container": "quay.io/cumulus/cellranger:|$.dataset.params.cellranger_version"
    ```
    Cell Ranger is proprietary 10x Genomics software distributed via the Cumulus
    images, which carry plain semver tags. This mirrors the `cirro/cellranger-*`
    processes (base image + selectable version). cellranger and cellranger-vdj
    share the `quay.io/cumulus/cellranger` image; cellranger-arc uses
    `quay.io/cumulus/cellranger-arc`.
  - **Full-image, free-text (simpleaf, kb)** — a required `container` free-text
    field holding the whole image reference, mapped directly:
    ```
    "container": "$.dataset.params.container"
    ```
    Used where there is no single stable base image / tag scheme to select from.

- **process-output.json** — identical across methods; a single
  `hot.Manifest` command to catalog the output files.

## Shared files (top of folder)

Both are referenced by every method's `process-definition.json` via
`s3://.../community/index_genomes/...`, so they are defined once:

- **process-compute.config** — common Nextflow resource block: retry on
  out-of-memory exit codes (137–140), scale cpus/memory by attempt, publish to
  `params.outdir`.
- **preprocess.py** — security guard run before launch. It scans every param
  whose name ends in `_args` against an allowlist regex and raises on shell
  metacharacters, so the free-text CLI-argument fields (`*_extra_args` and the
  Cell Ranger `*_mkref_args` / `*_mkgtf_args` / `*_mkvdjref_args`) cannot inject
  commands.

## Adding a new tool

1. Create `intake/genome_<tool>_index/process-definition.json` with
   `fileMappingRules` that match the tool's index files by name.
2. Create `community/index_genomes/<tool>/1.0/` with the four process files,
   copying an existing method and changing only: `id`/`name`/`dataType`/`desc`,
   `code.script` (the entrypoint filename in the workflow repo), the container
   field (version-tag or full-image, see above), and any tool-specific form
   fields. Keep the form param names in sync with the workflow's
   `nextflow.config` (`container`, `<tool>_extra_args`, etc.).
3. Add the matching entry script in the `CirroBio/nf-index-genome` repo and set
   `code.script` to its exact filename.
4. Reuse the shared `preprocess.py` and `process-compute.config` as-is. Keep any
   free-text CLI-argument field's name ending in `_args` so preprocess.py
   validates it.
