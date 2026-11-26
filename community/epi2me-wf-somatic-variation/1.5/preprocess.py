#!/usr/bin/env python3
"""Check the parameter combinations passed to epi2me-labs/wf-somatic-variation.

The form allows options which the workflow either ignores or rejects depending
on which subworkflows are enabled, whether a matched normal was provided, and
which reference genome is used. Fatal combinations raise here so the run fails
fast; ignored options only warn.
"""
import re
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from cirro.helpers.preprocess_dataset import PreprocessDataset

SUBWORKFLOWS = ["snv", "sv", "mod"]

# Options which have no effect unless the named subworkflow is enabled. Only
# options the user has to opt into are listed; warning about a form default
# being carried through would fire on every run.
SUBWORKFLOW_OPTIONS = {
    "snv": [
        "override_basecaller_cfg", "fast_mode", "phase_normal", "liquid_tumor",
        "skip_haplotype_filter", "normal_vcf", "hybrid_mode_vcf",
        "genotyping_mode_vcf",
    ],
    "sv": ["severus_args", "tr_bed", "pon_file"],
    "mod": ["force_strand", "modkit_args"],
}

# Options which only mean something when a matched normal BAM is provided
PAIRED_ONLY_OPTIONS = ["phase_normal", "normal_vcf"]

HG38_PATTERN = re.compile(r"(hg38|grch38)", re.IGNORECASE)
HG19_PATTERN = re.compile(r"(hg19|grch37|hs37)", re.IGNORECASE)


def is_set(value) -> bool:
    return value not in (None, "", False)


def infer_reference_build(ref: str):
    """Return 'hg38', 'hg19', or None based on the reference file name."""
    if HG38_PATTERN.search(ref):
        return "hg38"
    if HG19_PATTERN.search(ref):
        return "hg19"
    return None


def drop_empty_params(ds: "PreprocessDataset"):
    """Unset optional form fields arrive as empty strings; remove them so the
    workflow does not receive e.g. --bed ''."""
    for key, value in list(ds.params.items()):
        if value == "":
            ds.logger.info(f"Removing empty parameter: {key}")
            ds.remove_param(key, force=True)


def check_subworkflows(ds: "PreprocessDataset") -> list:
    enabled = [wf for wf in SUBWORKFLOWS if ds.params.get(wf)]
    if not enabled:
        raise ValueError(
            "No analysis subworkflow selected. Enable at least one of: "
            + ", ".join(SUBWORKFLOWS)
        )
    ds.logger.info(f"Enabled subworkflows: {', '.join(enabled)}")

    for wf, options in SUBWORKFLOW_OPTIONS.items():
        if wf in enabled:
            continue
        for opt in options:
            if is_set(ds.params.get(opt)):
                ds.logger.warning(
                    f"'{opt}' is set but the {wf} subworkflow is disabled; "
                    "it will have no effect"
                )
    return enabled


def check_tumor_only(ds: "PreprocessDataset", enabled: list) -> bool:
    """The normal BAM is optional; without it the workflow runs a reduced
    tumor-only analysis with different callers and skipped stages."""
    tumor_only = not is_set(ds.params.get("bam_normal"))
    if not tumor_only:
        return False

    ds.logger.info(
        "No normal BAM provided; running in tumor-only mode (ClairS-TO for "
        "small variants, Severus against a panel of normals, no differential "
        "methylation). 'normal_min_coverage' and 'germline' are unused."
    )
    for opt in PAIRED_ONLY_OPTIONS:
        if is_set(ds.params.get(opt)):
            ds.logger.warning(
                f"'{opt}' only applies to a tumor/normal pair and is ignored "
                "in tumor-only mode"
            )
    if "mod" in enabled and ds.params.get("diff_mod"):
        ds.logger.warning(
            "Differential modification analysis needs a matched normal; only "
            "the tumor bedMethyl will be produced"
        )
    if "snv" in enabled and ds.params.get("liquid_tumor"):
        ds.logger.warning(
            "ClairS-TO has no liquid tumor presets; 'liquid_tumor' will have "
            "no effect in tumor-only mode"
        )
    return True


def check_reference(ds: "PreprocessDataset", enabled: list, tumor_only: bool):
    ref = ds.params["ref"]
    build = infer_reference_build(ref)
    ds.logger.info(f"Reference build inferred from file name: {build or 'unknown'}")

    if ds.params.get("annotation") and build is None:
        ds.logger.warning(
            "SnpEff annotation is only performed for hg19/hg38 references; "
            f"could not confirm the build from '{ref}', so annotation may be skipped"
        )

    if "sv" in enabled:
        if build != "hg38" and not is_set(ds.params.get("pon_file")):
            ds.logger.warning(
                "The bundled Severus panel of normals is hg38 only; the "
                f"reference appears to be {build or 'a non-hg38 build'}, so SV "
                "calls will not be filtered against expected germline events"
            )
        if tumor_only:
            ds.logger.info(
                "Tumor-only SV calls are annotated with hg38 segmental "
                "duplication regions and need review of the SegDup flag"
            )


def check_option_combinations(ds: "PreprocessDataset", enabled: list):
    p = ds.params

    if "snv" in enabled and is_set(p.get("hybrid_mode_vcf")) and is_set(p.get("genotyping_mode_vcf")):
        raise ValueError(
            "Provide either a hybrid mode VCF or a genotyping mode VCF, not both"
        )

    if is_set(p.get("normal_vcf")):
        ds.logger.info(
            "The pre-computed normal VCF must be accompanied by its .tbi index "
            "in the same dataset"
        )

    if p.get("fast_mode") and not p.get("germline", True):
        ds.logger.warning(
            "'fast_mode' only changes how germline calls are made, but "
            "germline calling is disabled"
        )

    if "mod" in enabled:
        ds.logger.info(
            "Modified-base calling requires MM/ML tags in the input BAMs; this "
            "cannot be checked before the run starts"
        )


if __name__ == "__main__":

    # Imported here so that test_preprocess.py runs without the SDK
    from cirro.helpers.preprocess_dataset import PreprocessDataset

    ds = PreprocessDataset.from_running()

    drop_empty_params(ds)
    enabled = check_subworkflows(ds)
    tumor_only = check_tumor_only(ds, enabled)
    check_reference(ds, enabled, tumor_only)
    check_option_combinations(ds, enabled)

    # The HealthOmics pre-process Lambda reads params.json from the working directory and
    # fails the run when that file is absent. The SDK writes the file only when a parameter
    # changes, so a run that removes nothing must write it explicitly.
    ds.keep_params(list(ds.params.keys()))

    ds.logger.info("Final parameters:")
    for k, v in ds.params.items():
        ds.logger.info(f"{k}: {v}")
