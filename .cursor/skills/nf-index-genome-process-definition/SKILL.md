---
name: nf-index-genome-process-definition
description: Create Cirro process-definition.json entries that wrap nf-index-genome main_* Nextflow entrypoints for building genome indices. Use when adding or updating Cirro processes that call these index-building workflows.
---

# nf-index-genome Process Definitions

This skill explains how to create a new Cirro `process-definition.json` that points to one of the `main_*.nf` entrypoints in the `nf-index-genome` repository (`CirroBio/nf-index-genome`).

The goal is to produce a consistent process definition that:
- Exposes genome-indexing tools (Bowtie2, STAR, STAR2, Bismark, BWA, HISAT2) through the Cirro catalog.
- Correctly references the `CirroBio/nf-index-genome` GitHub repo and chosen `main_*.nf` entrypoint.
- Uses the standard S3 configuration artifacts for parameters, forms, compute settings, and UI output.

## When to Use This Skill

Use this skill when:
- Adding a **new genome index** process to the Cirro catalog backed by `nf-index-genome`.
- Updating an existing process to point to a different `main_*.nf` entrypoint in `nf-index-genome`.
- Creating environment- or version-specific variants (e.g. Bowtie2 vs STAR, or a new tagged version).

---

## Quick Recipe: New nf-index-genome Process

Follow this checklist to create a new `process-definition.json` pointing at `nf-index-genome`.

All nf-index-genome processes created with this skill must be saved under the `cirro/nf-core/index_genome` directory in this repository (for example, in a versioned subdirectory such as `cirro/nf-core/index_genome/1.0/`).

1. **Choose the index tool and entrypoint**
   - Select the appropriate `main_*.nf` file from `nf-index-genome`:
     - `main_bowtie2.nf` – Bowtie2 genome index.
     - `main_star.nf` – STAR genome index (splice-aware; requires GTF).
     - `main_star2.nf` – Alias for STAR 2.x (same as `main_star.nf`).
     - `main_bismark.nf` – Bismark bisulfite genome preparation.
     - `main_bwa.nf` – BWA genome index.
     - `main_hisat2.nf` – HISAT2 genome index (optional GTF/splice sites).
   - All are defined in the `CirroBio/nf-index-genome` repo (see README in that repo for details).

2. **Define basic process metadata**
   - **`id`**:
     - Use a stable, structured ID such as:
       - `process-cirro-genome-index-bowtie2-1-0`
       - `process-cirro-genome-index-star-1-0`
     - Pattern: `process-<org>-<short-name>-<version>` (avoid spaces and punctuation).
   - **`name`**:
     - Human-readable catalog label, e.g.:
       - `Build Genome Index (Bowtie2)`
       - `Build Genome Index (STAR)`
   - **`desc`**:
     - One sentence describing what is indexed and by which tool, e.g.:
       - `Build a Bowtie2 index from a reference genome FASTA (and optional annotation).`
   - **`dataType`** (optional but recommended):
     - Describe the resulting index:
       - `Genome Index (Bowtie2)`
       - `Genome Index (STAR)`, etc.
   - **`category`**:
     - Use or introduce a sensible category such as:
       - `Reference Data` or `Genome Indexing` (consistent with existing catalog taxonomy).
   - **`documentationUrl`**:
     - Point at the internal docs page for this process (e.g. a “Reference data / genome index” section).

3. **Set topology (optional but recommended)**
   - **`parentProcessIds`**:
     - If there is an upstream reference-genome curation or download step, add its process ID here.
     - Otherwise, use an empty list `[]`.
   - **`childProcessIds`**:
     - Add downstream processes that consume these genome indices, if appropriate.
     - For example, aligner or quantification pipelines that require a pre-built index.

4. **Wire in code location (nf-index-genome)**
   - Use the following `code` block pattern:

   ```json
   "code": {
     "repository": "GITHUBPUBLIC",
     "uri": "CirroBio/nf-index-genome",
     "version": "main",
     "script": "main_bowtie2.nf"
   }
   ```

   - Adjust `script` to the chosen entrypoint:
     - `main_bowtie2.nf`, `main_star.nf`, `main_star2.nf`, `main_bismark.nf`, `main_bwa.nf`, or `main_hisat2.nf`.
   - If pinning to a tag instead of a branch, set `version` to the tag (e.g. `v0.1`).

5. **Execution and compute defaults**
   - Always set:
     - `"executor": "NEXTFLOW"`
   - Use the standard compute configuration pattern:

   ```json
   "computeDefaults": [
     {
       "executor": "NEXTFLOW",
       "json": "s3://<RESOURCES_BUCKET>/<PROCESS_DIRECTORY>/process-compute.config",
       "name": "Default"
     }
   ]
   ```

   - Place any process-specific compute tuning in `process-compute.config`.

6. **Configuration artifacts (S3)**
   - Assume the standard S3 layout under `<PROCESS_DIRECTORY>`:
     - **`paramMapJson`**:
       - `s3://<RESOURCES_BUCKET>/<PROCESS_DIRECTORY>/process-input.json`
       - This file should define mappings for:
         - `--fasta` (required for all tools).
         - `--gtf` (required for STAR/STAR2; optional for HISAT2).
         - `--splicesites` (optional for HISAT2).
         - `--outdir` (optional; default `./results`).
     - **`formJson`**:
       - `s3://<RESOURCES_BUCKET>/<PROCESS_DIRECTORY>/process-form.json`
       - UI fields corresponding to the parameters above:
         - Reference FASTA upload / selection.
         - Optional GTF / splice site inputs when needed.
         - Optional outdir label or name.
     - **`webOptimizationJson`**:
       - `s3://<RESOURCES_BUCKET>/<PROCESS_DIRECTORY>/process-output.json`
       - Configure how index outputs are visualized or linked in Cirro.
     - **`preProcessScript`** (optional):
       - `s3://<RESOURCES_BUCKET>/<PROCESS_DIRECTORY>/preprocess.py`
       - Use this when you need to transform form inputs into Nextflow parameters, validate paths, or normalize genome names.

7. **Input behavior flags**
   - For genome indexing, typical defaults are:
     - `"allowMultipleSources": false` (usually one reference per run).
     - `"usesSampleSheet": false` (indices are usually not driven by per-sample metadata).
   - If a more advanced workflow accepts multiple genome builds or batch indices, you can set:
     - `"allowMultipleSources": true`

8. **(Optional) File mapping rules**
   - If you want Cirro to auto-map index outputs to data types or reference entities, define `fileMappingRules` similar to other processes:

   ```json
   "fileMappingRules": [
     {
       "description": "Bowtie2 Index Files",
       "min": 1,
       "fileNamePatterns": [
         {
           "exampleName": "genome.1.bt2",
           "description": "Bowtie2 Index",
           "sampleMatchingPattern": ".*/(?P<sampleName>[\\S ]*)\\.1\\.bt2"
         }
       ]
     }
   ]
   ```

   - Adjust description, `exampleName`, and regex pattern (`sampleMatchingPattern`) according to the actual index naming scheme.

---

## Example: Bowtie2 Genome Index Process Definition

This is a minimal example for a Bowtie2 index process backed by `main_bowtie2.nf`:

```json
{
  "id": "process-cirro-genome-index-bowtie2-1-0",
  "parentProcessIds": [],
  "childProcessIds": [],
  "dataType": "Genome Index (Bowtie2)",
  "name": "Build Genome Index (Bowtie2)",
  "desc": "Build a Bowtie2 index from a reference genome FASTA.",
  "executor": "NEXTFLOW",
  "category": "Reference Data",
  "documentationUrl": "https://<DOCS_SITE>/pipelines/catalog-reference/#genome-index-bowtie2",

  "code": {
    "repository": "GITHUBPUBLIC",
    "uri": "CirroBio/nf-index-genome",
    "version": "main",
    "script": "main_bowtie2.nf"
  },

  "allowMultipleSources": false,
  "usesSampleSheet": false,

  "computeDefaults": [
    {
      "executor": "NEXTFLOW",
      "json": "s3://<RESOURCES_BUCKET>/<PROCESS_DIRECTORY>/process-compute.config",
      "name": "Default"
    }
  ],

  "paramMapJson": "s3://<RESOURCES_BUCKET>/<PROCESS_DIRECTORY>/process-input.json",
  "formJson": "s3://<RESOURCES_BUCKET>/<PROCESS_DIRECTORY>/process-form.json",
  "webOptimizationJson": "s3://<RESOURCES_BUCKET>/<PROCESS_DIRECTORY>/process-output.json"
}
```

To create a STAR index process, copy this structure and change:

- `id`, `name`, `desc`, `dataType`, `documentationUrl` to reflect STAR.
- `code.script` to `main_star.nf` (or `main_star2.nf`).
- Ensure `process-input.json` and `process-form.json` expose both `--fasta` and `--gtf` as required fields.

---

## Summary

- Process definitions in this repo share a consistent JSON structure with stable IDs, clear metadata, Nextflow code references, standard S3 config paths, and optional topology/file-mapping metadata.
- `nf-index-genome` entrypoints are thin wrappers around single workflows that validate parameters, run the index workflow, and publish index files under `params.outdir`.
- This skill provides a repeatable template and checklist for wiring new Cirro `process-definition.json` files to any of the `main_*.nf` entrypoints in `nf-index-genome`, ensuring consistency across all genome-indexing processes.

