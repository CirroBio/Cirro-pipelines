# Cirro Pipelines — Claude Reference

## process-form.json Extended Syntax

Cirro's `process-form.json` files use standard JSON Schema draft-07 as a base, extended with custom properties for UI rendering, file resolution, and cloud path support. The top-level structure is:

```json
{
  "form": { ... },
  "ui": { ... }
}
```

---

### File / Path Resolution Properties

These appear on `string`-type properties to wire them up to file pickers.

| Property | Values | Purpose |
|----------|--------|---------|
| `pathType` | `"dataset"`, `"references"` | Where to resolve files from — a user dataset or a reference library |
| `file` | glob pattern (e.g. `"**/*.f*"`) | File pattern for the file picker |
| `folder` | glob pattern | Like `file` but selects a directory |
| `useS3Path` | `true` / `false` | Allow direct S3 path entry in addition to the picker UI |
| `multiple` | `true` / `false` | Allow selection of multiple files |
| `process` | array of process ID strings | Filter the dataset picker to only show datasets produced by the listed process IDs (pipeline process IDs or intake process IDs) |
| `mm` | `{ "matchBase": true }` | Controls file matching behavior (matchBase ignores leading path) |
| `schema` | file path string | Path to external JSON Schema for validating selected file contents |
| `mimetype` | MIME type string (e.g. `"text/plain"`) | Restricts accepted file types |

**Example:**
```json
"fasta": {
  "title": "Reference FASTA",
  "type": "string",
  "pathType": "dataset",
  "file": "**/*.f*",
  "description": "Reference genome FASTA."
}
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

## process-definition.json Structure

Key fields:

| Field | Purpose |
|-------|---------|
| `id` | Unique process identifier (e.g. `"process-cirro-genome-index-salmon-1-0"`) |
| `parentProcessIds` | Data type(s) this process accepts as input (e.g. `["genome_fasta"]`) |
| `childProcessIds` | Process IDs that can consume this process's output |
| `dataType` | Human-readable output data type label |
| `executor` | `"NEXTFLOW"` for pipeline processes, `"INGEST"` for intake processes |
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
