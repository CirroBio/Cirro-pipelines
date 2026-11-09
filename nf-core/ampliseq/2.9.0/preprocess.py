#!/usr/bin/env python3

from __future__ import annotations

from typing import TYPE_CHECKING

import pandas as pd

# The SDK is imported in __main__, so test_preprocess.py runs without it
if TYPE_CHECKING:
    from cirro.helpers.preprocess_dataset import PreprocessDataset


def read_samplesheet(path) -> pd.DataFrame:
    """Read the tab-separated samplesheet filled out in the form, keeping every value as text."""

    # The form field is marked required, but an empty table can still be submitted,
    # in which case Cirro either omits the param or writes an empty file
    if not path:
        raise ValueError("No samplesheet was provided: add one row per sample")
    try:
        samplesheet = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    except pd.errors.EmptyDataError:
        samplesheet = pd.DataFrame()
    if samplesheet.empty:
        raise ValueError("The samplesheet is empty: add one row per sample")
    return samplesheet


def check_unique_sample_ids(samplesheet: pd.DataFrame):
    repeated = samplesheet.loc[samplesheet["sampleID"].duplicated(), "sampleID"].unique()
    if len(repeated) > 0:
        raise ValueError(f"Each sampleID can be listed only once, but these are repeated: {', '.join(repeated)}")


def is_single_end(samplesheet: pd.DataFrame, only_use_forward_reads: bool) -> bool:
    """Whether to analyze every sample as single-end.

    nf-core/ampliseq runs either all samples paired-end or all single-end, and stops on a
    paired-end row without reverse reads.
    """
    return bool(only_use_forward_reads) or bool((samplesheet["reverseReads"] == "").any())


def read_metadata(path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def check_metadata_ids(samplesheet: pd.DataFrame, metadata: pd.DataFrame):
    """nf-core/ampliseq checks only that the metadata file exists, and QIIME2 fails late in
    the run when it lacks a sample, so check that every sampleID has a row in its first
    (ID) column."""

    missing = sorted(set(samplesheet["sampleID"]) - set(metadata.iloc[:, 0].str.strip()))
    if missing:
        raise ValueError(
            "The sample metadata file needs a row for every sampleID in the samplesheet, "
            f"with the ID in its first column. Missing: {', '.join(missing)}"
        )


# DADA2 reference taxonomies staged under the Cirro references bucket. The
# pipeline otherwise fetches these from zenodo when the run starts, which
# HealthOmics blocks under networkingMode RESTRICTED. Keyed by the
# dada_ref_taxonomy value they replace; anything not listed still downloads.
STAGED_DADA_REF = {
    "silva=138": (
        "ampliseq/silva-138.1/silva_nr99_v138.1_wSpecies_train_set.fa.gz",
        "ampliseq/silva-138.1/silva_species_assignment_v138.1.fa.gz",
    )
}


def use_staged_reference(ds: PreprocessDataset):
    """Swap a staged copy in for the selected DADA2 reference taxonomy.

    dada_ref_tax_custom takes precedence over dada_ref_taxonomy in the workflow,
    so setting both custom params and dropping dada_ref_taxonomy selects the
    staged files. These are the same DADA2-preformatted files the pipeline would
    have downloaded, so the custom path -- which skips FORMAT_TAXONOMY -- reads
    them as-is. A reference the user picked themselves always wins.
    """
    if ds.params.get("dada_ref_tax_custom"):
        return

    staged = STAGED_DADA_REF.get(ds.params.get("dada_ref_taxonomy"))
    if staged is None:
        return

    assign_tax, add_species = staged
    ds.logger.info(f"Using staged reference for {ds.params['dada_ref_taxonomy']}")
    ds.add_param("dada_ref_tax_custom", f"{ds.references_base}/{assign_tax}", overwrite=True)
    ds.add_param("dada_ref_tax_custom_sp", f"{ds.references_base}/{add_species}", overwrite=True)

    # Dead once a custom reference is set, and leaving it in overstates what ran
    ds.remove_param("dada_ref_taxonomy", force=True)


if __name__ == "__main__":

    from cirro.helpers.preprocess_dataset import PreprocessDataset

    ds = PreprocessDataset.from_running()

    samplesheet = read_samplesheet(ds.params.get("input"))
    check_unique_sample_ids(samplesheet)

    if is_single_end(samplesheet, ds.params.get("only_use_forward_reads", False)):
        ds.logger.info("Analyzing every sample as single-end, using only its forward reads")
        ds.add_param("single_end", True)

    # A form option only, which nf-core/ampliseq does not declare
    ds.remove_param("only_use_forward_reads", force=True)

    if ds.params.get("metadata"):
        check_metadata_ids(samplesheet, read_metadata(ds.params["metadata"]))

    # If the user set ignore_empty_input_files
    if ds.params.get("ignore_empty_input_files", False):
        # Also set the other ignore flags
        ds.add_param("ignore_failed_filtering", True, overwrite=True)
        ds.add_param("ignore_failed_trimming", True, overwrite=True)

    # The custom pickers are filled with the project's first reference as soon as they
    # are shown, and a run copied from before the box existed carries that value, so
    # they count only when the box is ticked
    if not ds.params.get("use_custom_taxonomy", False):
        ds.remove_param("dada_ref_tax_custom", force=True)
        ds.remove_param("dada_ref_tax_custom_sp", force=True)

    # A form option only, which nf-core/ampliseq does not declare
    ds.remove_param("use_custom_taxonomy", force=True)

    use_staged_reference(ds)

    # If a custom reference was provided
    if ds.params.get("dada_ref_tax_custom", False):
        # If the species-level reference was not provided
        if not ds.params.get("dada_ref_tax_custom_sp", False):
            # Turn off the species-level annotation
            ds.logger.info("Turning off species-level annotation, since no species reference was provided")
            ds.add_param("skip_dada_addspecies", True)

    # log
    ds.logger.info(ds.params)
