# Spatial Data Studio

Point this at one or more spatial datasets. It finds the spatial data inside each one,
loads it with the reader for its format, runs that format's preprocessing recipes, and
publishes a SpatialData checkpoint per dataset along with a browsable viewer over the
whole run.

Documentation: [cirrobio.github.io/spatial-data-studio](https://cirrobio.github.io/spatial-data-studio/)

## Data types

| Data type | Recognized by | Preprocessing |
| --------- | ------------- | ------------- |
| 10x Xenium | `experiment.xenium` | QC, filter, normalize, log1p, PCA, neighbors, UMAP (2D + 3D), Leiden, markers, cellular neighborhoods |
| 10x Visium HD | `binned_outputs/`, `segmented_outputs/`, `*feature_slice.h5` | scanpy Visium path (with highly-variable-gene selection), then cellular neighborhoods |
| 10x Visium | `spatial/scalefactors_json.json` plus a count matrix | as Visium HD |
| Vizgen MERSCOPE | `detected_transcripts.csv` plus `cell_by_gene.csv`/`cell_metadata.csv` | as Xenium |
| NanoString CosMx | `*exprMat_file.csv`, `*metadata_file.csv`, `*fov_positions_file.csv` | as Xenium |
| Curio Seeker | `anndata.h5ad` plus `Metrics.csv`/`cluster_assignment.txt` | as Visium |
| Steinbock | `cells.h5ad` plus `ome/` and a masks folder | z-score, PCA, neighbors, Leiden, UMAP, then neighborhoods |
| MCMICRO | `quantification/` plus `markers.csv` and `registration/`/`dearray/` | as Steinbock |

**Only the Xenium defaults have been run against real data.** Every other format's
recognition patterns and recipe are derived from `spatialdata_io`'s own format constants
and the standard analysis path for the modality, but no run has proven them. Treat the
parameter defaults for those types as a starting point.

## Input

Attach one or more datasets. Each becomes one root of the search, and its results are
published under its dataset name, so several datasets land in a single organized output.

Each root is walked and every folder matching one of the formats above becomes a dataset
(turn off **Recurse into subfolders** if each attached dataset is itself a single spatial
dataset). Matching is greedy: a folder that is recognized is not descended into, and when
two formats match the same folder the more specific one wins, so a Visium HD run does not
also register as the Visium-shaped matrix it contains.

## Output

```
index.html                         the viewer, at the output root
index.json                         lists every checkpoint below
results/<dataset>/<subfolders>/
    results-<hash>.sdata.zarr.zip    the full checkpoint
    lowres-<hash>.sdata.zarr.zip     same analysis, image pyramid capped
    plots/<NN>_<recipe>/figure.svg   static figures (also PDF)
    results.log                      always written
multiqc/multiqc_report.html        one report over the whole run
```

Both `.zarr.zip` checkpoints open in the Spatial Data Studio app (New Session, then Load)
when you want to compute on them further.

**A dataset that fails to load does not fail the run.** Its log is published where its
checkpoint would have gone, it shows in the MultiQC Datasets table as `failed`, and the
other datasets carry on. Likewise a recipe step that fails is recorded in the checkpoint's
history and the following steps still run, so the dataset is published with the part that
worked.

## Notes

- The analysis parameters apply per data type; each description lists which types use it.
  A parameter is ignored by a type that does not use it.
- Analysis tasks install their pinned Python dependencies at startup, and the viewer is
  downloaded from the latest Spatial Data Studio release, so the compute environment
  needs outbound network access.
