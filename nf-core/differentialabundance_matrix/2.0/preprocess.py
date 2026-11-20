#!/usr/bin/env python3
"""Set up nf-core/differentialabundance 2.0 to test one matrix of a combined dataset.

A combined dataset (written by the Cirro Multi-omics app) holds samples.csv
(sample_id, subject_id, then one column per variable), combined.json (which
records each variable's type) and, per modality, <modality>/raw.csv and
<modality>/prepared.csv: features in rows, a feature_id column, then one
column per sample that has the modality.

From these and the form, this writes the pipeline's samplesheet, a one-contrast
contrasts file and a copy of the matrix to the dataset's config/ folder, and
sets the method params that nf-core's analysis profiles would otherwise set.
Whitespace in sample IDs is replaced with '_' in both the samplesheet and the
matrix copy, since the pipeline rejects it.
"""

from __future__ import annotations

import json
import re
from typing import TYPE_CHECKING

import pandas as pd

# The SDK is imported in __main__, so test_preprocess.py runs without it
if TYPE_CHECKING:
    from cirro.helpers.preprocess_dataset import PreprocessDataset

# The samplesheet column the pipeline matches to the matrix columns. dream reads
# it as "sample" whatever observations_id_col says.
SAMPLE_COLUMN = "sample"

# The matrix's feature ID column. dream reads the matrix's IDs from "gene_id"
# whatever features_id_col says, so the matrix is copied with this header.
FEATURE_COLUMN = "gene_id"

# The params of nf-core's rnaseq, rnaseq_dream, generic_matrix and
# generic_matrix_dream profiles (conf/profile/ in 2.0.0) that differ between
# them, by (values, whether there is a random effect)
METHOD_PARAMS = {
    ("counts", False): {
        "paramset_name": "rnaseq",
        "study_type": "rnaseq",
        "study_abundance_type": "counts",
        "differential_method": "deseq2",
        "exploratory_assay_names": "raw,normalised,variance_stabilised",
        "exploratory_final_assay": "variance_stabilised",
        "exploratory_log2_assays": "raw,normalised",
        "filtering_min_abundance": 1,
    },
    ("counts", True): {
        "paramset_name": "rnaseq_dream",
        "study_type": "rnaseq",
        "study_abundance_type": "counts",
        "differential_method": "dream",
        "dream_apply_voom": True,
        "exploratory_assay_names": "raw,normalised",
        "exploratory_final_assay": "normalised",
        "exploratory_log2_assays": "raw",
        "filtering_min_abundance": 1,
    },
    ("log_abundance", False): {
        "paramset_name": "generic_matrix",
        "study_type": "generic_matrix",
        "study_abundance_type": "abundance",
        "differential_method": "limma",
        "limma_use_voom": False,
        "exploratory_assay_names": "raw",
        "exploratory_final_assay": "raw",
        # An empty list: the values are already on a log scale
        "exploratory_log2_assays": "[]",
        # Log-scale values can be negative, so no minimum applies
        "filtering_min_abundance": False,
    },
    ("log_abundance", True): {
        "paramset_name": "generic_matrix_dream",
        "study_type": "generic_matrix",
        "study_abundance_type": "abundance",
        "differential_method": "dream",
        "dream_apply_voom": False,
        "exploratory_assay_names": "raw",
        "exploratory_final_assay": "raw",
        "exploratory_log2_assays": "[]",
        "filtering_min_abundance": False,
    },
}

# DESeq2's vst() fits its trend on this many features (deseq2_vst_nsub's default)
# and stops when the matrix has fewer
DESEQ2_VST_NSUB = 1000

# Form params that only this script reads
FORM_PARAMS = [
    "combined_dataset",
    "matrix_in_dataset",
    "values",
    "variable",
    "reference",
    "target",
    "covariates",
    "random_effect",
]

# A column the R model formula can name without quoting
R_NAME = re.compile(r"^[A-Za-z][A-Za-z0-9._]*$")

# The pipeline's samplesheet schema (assets/schema_input.json) requires sample to match ^\S+$
WHITESPACE = re.compile(r"\s+")


def parse_columns(value) -> list[str]:
    """The column names in a comma-separated form field, which may be unset."""
    if not value:
        return []
    return [name.strip() for name in str(value).split(",") if name.strip()]


def check_model_columns(samples: pd.DataFrame, columns: list[str]):
    """Fail clearly on a model column that samples.csv lacks or R can't name in a formula."""
    for column in columns:
        if column not in samples.columns:
            raise ValueError(f"samples.csv has no column '{column}'. Its columns are: {', '.join(samples.columns)}")
        if not R_NAME.match(column):
            raise ValueError(
                f"Column '{column}' can't be used in a model: a name must start with a letter "
                "and hold only letters, digits, '.' and '_'"
            )


def make_samplesheet(
    samples: pd.DataFrame,
    matrix_samples: list[str],
    model_columns: list[str],
    variable: str,
    reference: str,
    target: str,
) -> pd.DataFrame:
    """
    The samples that have this matrix and a value for every model column, with
    sample_id renamed to the pipeline's sample column. The pipeline stops on a
    sample missing from the matrix or a missing value in the design.
    """
    if SAMPLE_COLUMN in samples.columns:
        raise ValueError(f"samples.csv has a variable named '{SAMPLE_COLUMN}', which the pipeline reserves for sample IDs")

    kept = samples[samples["sample_id"].isin(matrix_samples)]
    for column in model_columns:
        kept = kept[kept[column] != ""]

    for level in (reference, target):
        if not (kept[variable] == level).any():
            raise ValueError(
                f"No sample with this matrix and a value for every model column ({', '.join(model_columns)}) "
                f"has {variable} = '{level}'"
            )

    return kept.rename(columns={"sample_id": SAMPLE_COLUMN})


def pipeline_sample_names(sample_ids: list[str]) -> dict[str, str]:
    """Each sample ID with its whitespace replaced, as the pipeline's samplesheet requires."""
    names = {sample_id: WHITESPACE.sub("_", sample_id) for sample_id in sample_ids}
    if len(set(names.values())) < len(names):
        raise ValueError("Two sample IDs are the same once their spaces are replaced with '_'")
    return names


def r_name(text: str) -> str:
    """R's make.names() for ASCII text, which the pipeline applies to design coefficients."""
    name = re.sub(r"[^A-Za-z0-9._]", ".", text)
    if not re.match(r"[A-Za-z]|\.(?![0-9])", name):
        name = "X" + name
    return name


def model_term(column: str, variable_types: dict) -> str:
    """A covariate as a formula term: a factor unless combined.json records it as numeric."""
    if variable_types.get(column) == "numeric":
        return column
    return f"factor({column})"


def make_contrast(
    variable: str,
    reference: str,
    target: str,
    covariates: list[str],
    random_effect: str | None,
    variable_types: dict,
) -> dict:
    """
    One entry of the YAML contrasts file. A plain comparison factors every
    blocking variable, so a numeric covariate or a random effect needs a formula.
    """
    contrast_id = re.sub(r"[^A-Za-z0-9_-]+", "_", f"{variable}_{target}_vs_{reference}")

    numeric_covariates = [c for c in covariates if variable_types.get(c) == "numeric"]
    if not random_effect and not numeric_covariates:
        contrast = {"id": contrast_id, "comparison": [variable, reference, target]}
        if covariates:
            contrast["blocking_factors"] = covariates
        return contrast

    # With no intercept and the variable first, every level of the variable is
    # a coefficient, so the contrast is the difference of two of them
    terms = [f"factor({variable})"] + [model_term(c, variable_types) for c in covariates]
    if random_effect:
        terms.append(f"(1 | {random_effect})")
    return {
        "id": contrast_id,
        "formula": "~ 0 + " + " + ".join(terms),
        "make_contrasts_str": (
            f"{r_name(f'factor({variable}){target}')} - {r_name(f'factor({variable}){reference}')}"
        ),
    }


def contrasts_yaml(contrast: dict) -> str:
    # JSON is valid YAML, and spares the quoting of values such as '~ 0 + ...'
    return json.dumps({"contrasts": [contrast]}, indent=2) + "\n"


def prepare_matrix(matrix: pd.DataFrame, sample_names: dict[str, str]) -> pd.DataFrame:
    """The matrix with its feature IDs under the pipeline's column, limited to these samples and renamed."""
    if matrix.columns[0] != "feature_id":
        raise ValueError(f"The matrix's first column must be feature_id, not '{matrix.columns[0]}'")
    return (
        matrix[["feature_id"] + list(sample_names)]
        .rename(columns={"feature_id": FEATURE_COLUMN, **sample_names})
    )


def method_params(values: str, repeated: bool, n_features: int) -> dict:
    """The params that choose and configure the method, for this kind of values and design."""
    params = dict(METHOD_PARAMS[(values, repeated)])
    if params["differential_method"] == "deseq2" and n_features < DESEQ2_VST_NSUB:
        # Small count matrices, such as microbial taxa, are stabilised with rlog instead
        params["deseq2_vs_method"] = "rlog"
    return params


def read_text(uri: str) -> str:
    # Imported here so that test_preprocess.py runs without the SDK or boto3
    import boto3
    from cirro.models.s3_path import S3Path

    path = S3Path(uri)
    return boto3.client("s3").get_object(Bucket=path.bucket, Key=path.key)["Body"].read().decode()


def write_text(uri: str, content: str):
    import boto3
    from cirro.models.s3_path import S3Path

    path = S3Path(uri)
    boto3.client("s3").put_object(Bucket=path.bucket, Key=path.key, Body=content.encode())


def run(ds: PreprocessDataset):
    root = ds.params["combined_dataset"].rstrip("/")
    values = ds.params["values"]
    variable = ds.params["variable"]
    reference = str(ds.params["reference"])
    target = str(ds.params["target"])
    covariates = parse_columns(ds.params.get("covariates"))
    random_effect = ds.params.get("random_effect") or None

    if values not in ("counts", "log_abundance"):
        raise ValueError(f"Values must be counts or log_abundance, not '{values}'")

    # Read every value as text, so that sample IDs such as "01" are kept as-is
    samples = pd.read_csv(f"{root}/samples.csv", dtype=str, keep_default_na=False)
    model_columns = [variable] + covariates + ([random_effect] if random_effect else [])
    check_model_columns(samples, model_columns)

    matrix_uri = f"{root}/{ds.params['matrix_in_dataset']}"
    ds.logger.info(f"Testing {matrix_uri}")
    matrix = pd.read_csv(matrix_uri, dtype={"feature_id": str})
    samplesheet = make_samplesheet(samples, list(matrix.columns[1:]), model_columns, variable, reference, target)
    ds.logger.info(f"{len(samplesheet):,} of {len(samples):,} samples have this matrix and a value for every model column")
    sample_names = pipeline_sample_names(list(samplesheet[SAMPLE_COLUMN]))
    samplesheet = samplesheet.assign(**{SAMPLE_COLUMN: samplesheet[SAMPLE_COLUMN].map(sample_names)})

    combined = json.loads(read_text(f"{root}/combined.json"))
    variable_types = {v["name"]: v.get("type") for v in combined.get("variables", [])}
    contrast = make_contrast(variable, reference, target, covariates, random_effect, variable_types)
    ds.logger.info(f"Contrast: {contrast}")

    # Write to the dataset's config/ folder (mapped in process-input.json)
    samplesheet.to_csv(ds.params["input"], index=None)
    prepare_matrix(matrix, sample_names).to_csv(
        ds.params["matrix"], sep="\t", index=None, na_rep="NA"
    )
    write_text(ds.params["contrasts"], contrasts_yaml(contrast))

    for name, value in method_params(values, random_effect is not None, len(matrix)).items():
        ds.add_param(name, value, overwrite=True)

    # nf-schema rejects params the pipeline doesn't define
    for name in FORM_PARAMS:
        ds.remove_param(name, force=True)


if __name__ == "__main__":
    from cirro.helpers.preprocess_dataset import PreprocessDataset

    run(PreprocessDataset.from_running())
