#!/usr/bin/env python3
"""Check the parameter combinations passed to epi2me-labs/wf-human-variation.

The form allows options which the workflow either ignores or rejects depending
on which subworkflows are enabled and which reference genome is used. Fatal
combinations raise here so the run fails fast; ignored options only warn.
"""
import re
from cirro.helpers.preprocess_dataset import PreprocessDataset

SUBWORKFLOWS = ["sv", "snp", "cnv", "str", "mod"]

# Options which have no effect unless the named subworkflow is enabled
SUBWORKFLOW_OPTIONS = {
    "snp": ["GVCF", "override_basecaller_cfg"],
    "sv": ["tr_bed", "sniffles_args"],
    "cnv": ["use_qdnaseq", "qdnaseq_bin_size", "spectre_args"],
    "mod": ["force_strand", "modkit_args"],
}

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


def drop_empty_params(ds: PreprocessDataset):
    """Unset optional form fields arrive as empty strings; remove them so the
    workflow does not receive e.g. --bed ''."""
    for key, value in list(ds.params.items()):
        if value == "":
            ds.logger.info(f"Removing empty parameter: {key}")
            ds.remove_param(key, force=True)


def drop_inactive_dependents(ds: PreprocessDataset):
    """Fields shown only when a toggle is on are meaningless when it is off."""
    for toggle, dependent in [
        ("downsample_coverage", "downsample_coverage_target"),
        ("use_qdnaseq", "qdnaseq_bin_size"),
    ]:
        if not ds.params.get(toggle) and dependent in ds.params:
            ds.logger.info(f"Removing {dependent} because {toggle} is not enabled")
            ds.remove_param(dependent, force=True)


def check_subworkflows(ds: PreprocessDataset) -> list:
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


def check_reference(ds: PreprocessDataset, enabled: list):
    ref = ds.params["ref"]
    build = infer_reference_build(ref)
    ds.logger.info(f"Reference build inferred from file name: {build or 'unknown'}")

    if "str" in enabled and build != "hg38":
        msg = "STR genotyping (Straglr) requires an hg38 reference"
        if build is None:
            ds.logger.warning(
                f"{msg}; could not confirm the build from '{ref}'"
            )
        else:
            raise ValueError(f"{msg}, but the reference appears to be {build}")

    if "cnv" in enabled and build is None:
        ds.logger.warning(
            "CNV calling requires an hg19 or hg38 reference; could not "
            f"confirm the build from '{ref}'"
        )

    if ds.params.get("annotation") and build is None:
        ds.logger.warning(
            "SnpEff annotation is only performed for hg19/hg38 references; "
            f"could not confirm the build from '{ref}', so annotation may be skipped"
        )


def check_option_combinations(ds: PreprocessDataset, enabled: list):
    p = ds.params

    if p.get("phased") and "snp" not in enabled:
        raise ValueError("Phasing requires the SNP subworkflow to be enabled")

    if p.get("override_basecaller_cfg") == "custom":
        raise ValueError(
            "override_basecaller_cfg='custom' requires a custom Clair3 model, "
            "which this pipeline configuration does not expose"
        )

    if p.get("annotation") and "snp" not in enabled and "sv" not in enabled:
        ds.logger.warning(
            "SnpEff annotation is enabled but neither the SNP nor SV "
            "subworkflow is selected; there are no variant calls to annotate"
        )

    if is_set(p.get("partner")):
        if not p.get("annotation"):
            ds.logger.warning(
                f"Partner export ({p['partner']}) expects annotated variant "
                "calls, but annotation is disabled"
            )
        if "snp" not in enabled:
            ds.logger.warning(
                f"Partner export ({p['partner']}) is intended for small-variant "
                "output, but the SNP subworkflow is disabled"
            )

    if p.get("use_qdnaseq") and is_set(p.get("spectre_args")):
        ds.logger.warning(
            "spectre_args is set but use_qdnaseq is enabled; Spectre will not "
            "run and the arguments will be ignored"
        )

    if is_set(p.get("sex")) and not ({"str", "cnv"} & set(enabled)):
        ds.logger.warning(
            "'sex' is only used by the STR and CNV subworkflows, neither of "
            "which is enabled"
        )

    if is_set(p.get("alignment_report_coverage_threshold")) and not is_set(p.get("bed")):
        ds.logger.info(
            "alignment_report_coverage_threshold only applies when a target "
            "regions BED is provided"
        )

    if p.get("downsample_coverage"):
        target = p.get("downsample_coverage_target")
        min_cov = p.get("bam_min_coverage")
        if target is None:
            ds.logger.warning(
                "downsample_coverage is enabled without a target; the workflow "
                "default will be used"
            )
        elif min_cov is not None and target < min_cov:
            ds.logger.warning(
                f"downsample_coverage_target ({target}) is below "
                f"bam_min_coverage ({min_cov}); the downsampled BAM will fail "
                "the minimum coverage check"
            )

    if "mod" in enabled:
        ds.logger.info(
            "Modified-base calling requires MM/ML tags in the input BAM; this "
            "cannot be checked before the run starts"
        )


if __name__ == "__main__":
    ds = PreprocessDataset.from_running()

    drop_empty_params(ds)
    drop_inactive_dependents(ds)
    enabled = check_subworkflows(ds)
    check_reference(ds, enabled)
    check_option_combinations(ds, enabled)

    # The HealthOmics pre-process Lambda reads params.json from the working directory and
    # fails the run when that file is absent. The SDK writes the file only when a parameter
    # changes, so a run that removes nothing must write it explicitly.
    ds.keep_params(list(ds.params.keys()))

    ds.logger.info("Final parameters:")
    for k, v in ds.params.items():
        ds.logger.info(f"{k}: {v}")
