# Cirro Pipelines — Claude Reference

## process-form.json Extended Syntax

Cirro's `process-form.json` files use standard JSON Schema draft-07 as a base, extended with custom properties for UI rendering, file selection, and structured inputs. The top-level structure is:

```json
{
  "form": { ... },
  "ui": { ... },
  "structuredInputs": { ... }
}
```

---

### Selecting Files and Folders

Set `format` on a `string` property to let the user pick a file or folder from the
project's datasets (the same values as nf-schema):

| `format` | The user can select |
|----------|---------------------|
| `file-path`, `file-path-pattern` | A file |
| `directory-path` | A folder (including a whole dataset) |
| `path` | A file or a folder |

| Property | Purpose |
|----------|---------|
| `pattern` | Regex that filters the files offered, matched against `data/<path>`. Anchor it at the end (`"\\.bam$"`), never the start: the value the form stores is a `cirro::` URI, and the pattern is validated against it |
| `multiple` | `true` allows several picks. The value is a JSON array, so map it with `[*]` (`"$.dataset.params.fastq[*]"`) to pass a comma-separated string |

The workflow receives the full S3 path; folders have no trailing slash. `default` does not
pre-select anything in this picker.

**Example:**
```json
"gtf": {
  "title": "Annotation File",
  "type": "string",
  "format": "file-path",
  "pattern": "\\.gtf(\\.gz)?$"
}
```

---

### One File Picker for the Whole Form

The portal picks a mode for the whole form, not per field:

- **Path mode** — the form has a `structuredInputs` entry, or a `format` path field anywhere
  (including inside `dependencies` and `if/then` branches). There is no "Datasets to use"
  selector: the run's input datasets (`$.inputs`) are whichever datasets the form's picks
  point into. Legacy `file` globs use the new picker too, and their `default` stops resolving.
- **Dataset mode** — everything else. The user selects input datasets first, and `file`
  globs offer files from those datasets only.

Convert a pipeline as a whole: replace every `file` glob with `format: file-path` and a
`pattern`, and set `usesSampleSheet` to `false`. Path mode loads no samples, so a pipeline
that still requires them can never be run.

---

### Legacy File Keywords

| Property | Values | Purpose |
|----------|--------|---------|
| `pathType` | `"dataset"`, `"reference"` (also spelled `"references"`) | `dataset` without `file` picks a whole dataset (its S3 path). `reference` picks a file from the project's references, filtered by the `file` glob |
| `process` | array of process ID strings | Limits a `pathType: "dataset"` picker to datasets produced by these processes (pipeline or intake process IDs). The new picker has no equivalent, so dataset selectors stay on `pathType` |
| `file` | glob pattern (e.g. `"**/*.f*"`) | Legacy file picker; use `format: file-path` instead |
| `mm` | `{ "matchBase": true }` | Matching options for a `file` glob |
| `multiple` | `true` / `false` | Allow several selections; the value is a JSON array |
| `useS3Path` | `true` / `false` | Ignored by the portal; drop it when converting a form |

**Example:**
```json
"genome_index": {
  "title": "STAR Index",
  "type": "string",
  "pathType": "dataset",
  "process": ["process-cirro-genome-index-star-1-0", "genome_star_index"]
}
```

---

### Structured Inputs (Sample Sheets)

A `structuredInputs` entry describes a file the user fills out as a table in the form, and
Cirro writes it to the dataset's `config/` folder at launch. This replaces building a
samplesheet in `preprocess.py` from `ds.files` and `ds.samplesheet`.

```json
"form": {
  "type": "object",
  "required": ["input"],
  "properties": {
    "input": {"type": "string", "schema": "pipeline_samplesheet"}
  }
},
"structuredInputs": {
  "pipeline_samplesheet": {
    "fileFormat": "CSV",
    "schema": {
      "title": "Samplesheet",
      "type": "array",
      "items": {
        "type": "object",
        "properties": {
          "sample": {"type": "string", "title": "Sample", "pattern": "^\\S+$"},
          "fastq_1": {"type": "string", "format": "file-path", "title": "FASTQ (R1)"}
        },
        "required": ["sample", "fastq_1"]
      }
    }
  }
}
```

Map the file's path in `process-input.json`; the form field itself never reaches
`$.dataset.params`:

```json
"input": "$.structuredInputs.pipeline_samplesheet.path"
```

- **Naming** — the key becomes `config/<key>.csv` (or `.tsv`/`.json`), plus a raw copy at
  `config/<key>.data.json` that "Copy parameters" reloads. Never use `samplesheet`, `files`,
  `params`, `metadata` or `transform-params`: nothing stops them from colliding with Cirro's
  own config files. Use `pipeline_samplesheet`.
- **Columns** — one per property, in the order listed. A path `format` gives a file or folder
  picker (one value per cell; use `fastq_1`/`fastq_2`, never an array), `enum` gives a
  dropdown, and `integer`/`number`/`boolean` are checked. `object` and `array` columns are
  not supported. `pattern`, `minimum`/`maximum`, `minLength`/`maxLength` and a string
  `errorMessage` are checked against each row, and a file column's `pattern` sees the
  dataset-relative path. Column `default`, `enumNames` and `oneOf` labels are ignored, and
  array-level `minItems`/`maxItems` are not enforced.
- **Auto-fill** — users can import rows from a project Sheet, a CSV in a dataset, or a
  dataset's file names. The generated rows are matched to columns by name: `sample` takes
  the sample name, and `fastq_1`/`fastq_2` take paired reads. When the dataset has I1/I2
  index reads, they become rows of their own, so exclude them in the FASTQ columns' `pattern` (see
  `nf-core/rnaseq/3.1/process-form.json`).
- **Validation is advisory** — invalid rows don't block Run, nothing re-checks them on the
  server, and `required` on the form field doesn't mean at least one row. With no rows, Cirro
  writes an empty file (no header), or omits the param if the field was never touched.
  `preprocess.py` should fail clearly in both cases (see `cirro/cellranger-gex/1.1/preprocess.py`).
- **File contents** — header row first, columns in schema order, blank cells empty, lines
  ending in `\r\n`. File and folder cells hold absolute S3 URIs (an access-point alias for
  datasets shared from another project).

To change the file before the run, read and overwrite it at its mapped param. Read with
`dtype=str, keep_default_na=False` so that sample names such as `01` stay as text:

```python
samplesheet = pd.read_csv(ds.params["input"], dtype=str, keep_default_na=False)
samplesheet.to_csv(ds.params["input"], index=None)
```

---

### UI Object Properties

The top-level `"ui"` object mirrors the `form.properties` structure and contains rendering hints keyed by property name.

| Property | Values | Purpose |
|----------|--------|---------|
| `ui:widget` | `"checkboxes"`, `"radio"` | Override default widget (e.g. multi-select checkboxes instead of a dropdown) |
| `ui:placeholder` | string | Placeholder/hint text shown in empty inputs |

**Example:**
```json
"ui": {
  "extra_args": {
    "ui:placeholder": "--genomeSAindexNbases 11"
  },
  "tools": {
    "ui:widget": "checkboxes"
  }
}
```

---

### Enum Display

| Property | Values | Purpose |
|----------|--------|---------|
| `enumNames` | array of strings (parallel to `enum`) | Human-readable display labels for enum values |

**Example:**
```json
"genome": {
  "type": "string",
  "enum": ["GRCh38", "GRCh37"],
  "enumNames": ["Homo sapiens (GRCh38)", "Homo sapiens (GRCh37)"]
}
```

---

### Conditional / Dependency Logic

Standard JSON Schema conditionals are used for dynamic form behavior:

- **`dependencies`** — show/require fields based on another field's value (using `oneOf`)
- **`allOf` + `if/then`** — show different sets of fields based on enum selections

**Example:**
```json
"allOf": [
  {
    "if": { "properties": { "analysis_type": { "const": "Germline" } } },
    "then": { "properties": { "tools": { ... } } }
  }
]
```

---

### Top-Level `form` Object

The `form` object is a JSON Schema `object` type. It may include:
- `name`: a string identifier for the form
- Standard `properties`, `required`, `title`, `description`

---

## process-input.json — Container Version Pattern

Container image tags are interpolated from form params using a pipe syntax:

```json
"container": "quay.io/biocontainers/salmon:|$.dataset.params.salmon_version"
```

The `|` separates the static image prefix from a JSONPath expression that resolves the tag at runtime from the user-supplied parameter.

---

## process-input.json — JSONPath Roots

The mapping resolves against the launch payload, whose roots are `$.dataset`,
`$.project`, `$.process`, `$.inputs`, `$.structuredInputs` and `$.referencesBucketUri`.
Use `$.dataset.*` — the older `$.params.dataset.paramJson.*` / `$.params.dataset.s3` /
`$.params.inputs[*]` spellings still resolve but are not the documented form.

| Use | Path |
|-----|------|
| A form parameter | `$.dataset.params.<field>` (nested groups: `$.dataset.params.<group>.<field>`) |
| Output directory | `$.dataset.dataPath` — always this, never `$.dataset.s3\|/data/` |
| Dataset root on S3 | `$.dataset.s3` |
| An input dataset | `$.inputs[0].dataPath`, `$.inputs[0].s3`, or `$.inputs[*].dataPath` for all |
| A structured input file | `$.structuredInputs.<key>.path` |
| Cirro's shared references bucket | `$.referencesBucketUri` (e.g. `$.referencesBucketUri\|/igenomes/`) |

`$.dataset.s3` is still right for the dataset's other folders — `|/config/<name>`,
`|/web/`, `|/logs/`. Only the output directory has a dedicated field, and only that
field is matched by name when a run is normalized for HealthOmics.

---

## preprocess.py — Staging Files to the Dataset

A `preprocess.py` must **not** write workflow inputs to its working directory. Nothing
stages that directory for the run, so a bare filename cannot resolve on HealthOmics,
and the file is lost rather than recorded alongside the results. Write to the dataset's
own `config/` folder instead.

1. Declare the path in `process-input.json`:

   ```json
   "input": "$.dataset.s3|/config/manifest.csv"
   ```

2. Write to that resolved path in `preprocess.py`. No `add_param` call is needed — the
   mapping has already put the value in `ds.params`:

   ```python
   # Write to the dataset's config/ folder (mapped in process-input.json)
   manifest.to_csv(ds.params["input"], index=None)
   ```

**Never use these four names** — `config/` is also where Cirro writes the dataset's own
config, and `PreprocessDataset.from_running()` reads all four at startup:

| Reserved | Read into |
|----------|-----------|
| `config/samplesheet.csv` | `ds.samplesheet` |
| `config/files.csv` | `ds.files` |
| `config/params.json` | `ds.params` |
| `config/metadata.json` | `ds.metadata` |

Where a pipeline's own input is conceptually a samplesheet, `config/pipeline_samplesheet.csv`
is the convention. Prefer a [structured input](#structured-inputs-sample-sheets) over writing
one here: building a samplesheet from `ds.files` and `ds.samplesheet` is deprecated.

**Conditional writes:** the mapping populates the param on every run, so a file that is
only written under some conditions must have its param removed otherwise — a param
pointing at a missing object fails the run:

```python
if probe_barcodes is None:
    ds.remove_param("probe_barcodes")
```

**Non-pandas writes:** `pandas.to_csv` accepts an S3 URI, but `open()` does not. Use
`boto3` for YAML, JSON or plain text (see `community/pgap/1.0/preprocess.py`).

**Files referenced from inside a staged file** must carry the full S3 URI, not a relative
name — derive the folder from the param that is mapped:

```python
submol_uri = ds.params["sample_sheet"].rsplit("/", 1)[0] + "/submol.yaml"
```

Cromwell (`"executor": "CROMWELL"`) processes are the exception: `inputs.json` and
`options.json` are picked up from the working directory by convention and have no
`process-input.json` param to redirect.

---

## process-definition.json Structure

Key fields:

| Field | Purpose |
|-------|---------|
| `id` | Unique process identifier (e.g. `"process-cirro-genome-index-salmon-1-0"`) |
| `parentProcessIds` | Data type(s) this process accepts as input (e.g. `["genome_fasta"]`) |
| `childProcessIds` | Process IDs that can consume this process's output |
| `dataType` | Human-readable output data type label |
| `executor` | `"NEXTFLOW"` for pipeline processes, `"INGEST"` for intake processes |
| `executionModes` | Optional. Where the process can run: `["STANDARD"]` (AWS Batch), `["HEALTHOMICS"]` (AWS HealthOmics), or both. Omit to declare Standard only |
| `usesSampleSheet` | Deprecated. `true` makes the user select Cirro samples; it must be `false` once the form is in [path mode](#one-file-picker-for-the-whole-form) |
| `category` | UI grouping (e.g. `"Reference Data"`) |
| `code.repository` | `"GITHUBPUBLIC"` or `"NA"` |
| `code.uri` | GitHub repo (e.g. `"CirroBio/nf-index-genome"`) |
| `code.script` | Entry point Nextflow script (e.g. `"main_salmon.nf"`) |
| `paramMapJson` | S3 path to `process-input.json` (use `s3://<RESOURCES_BUCKET>/<PROCESS_DIRECTORY>/...`) |
| `formJson` | S3 path to `process-form.json` |
| `webOptimizationJson` | S3 path to `process-output.json` |

---

## intake/process-definition.json Structure

Intake processes use `"executor": "INGEST"` and define `fileMappingRules` instead of code/compute fields.

### fileMappingRules

Each rule specifies a required file group:

```json
"fileMappingRules": [
  {
    "description": "Human-readable label",
    "min": 1,
    "isSample": false,
    "fileNamePatterns": [
      {
        "exampleName": "genome.sa",
        "description": "Short description",
        "sampleMatchingPattern": "[\\S ]*.sa"
      }
    ]
  }
]
```

| Field | Purpose |
|-------|---------|
| `min` | Minimum number of matching files required |
| `isSample` | `false` means files are not per-sample (use for reference/index files) |
| `fileNamePatterns` | Array of regex patterns; a file matching any one satisfies the rule |
| `sampleMatchingPattern` | Regex; use `(?P<sampleName>...)` capture group when `isSample: true`, omit it when `isSample: false` |

Use `"fileMappingRules": []` when no specific files are required.
