#!/usr/bin/env python3


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


def staged_reference(params: dict, references_base: str):
    """The custom-reference params that select a staged copy of the chosen DADA2
    reference taxonomy, or None when there is nothing to swap in.

    dada_ref_tax_custom takes precedence over dada_ref_taxonomy in the workflow,
    so setting both custom params and dropping dada_ref_taxonomy selects the
    staged files. These are the same DADA2-preformatted files the pipeline would
    have downloaded, so the custom path -- which skips FORMAT_TAXONOMY -- reads
    them as-is. A reference the user picked themselves always wins.
    """
    if params.get("dada_ref_tax_custom"):
        return None

    staged = STAGED_DADA_REF.get(params.get("dada_ref_taxonomy"))
    if staged is None:
        return None

    assign_tax, add_species = staged
    return {
        "dada_ref_tax_custom": f"{references_base}/{assign_tax}",
        "dada_ref_tax_custom_sp": f"{references_base}/{add_species}",
    }


if __name__ == "__main__":

    # Imported here so that test_preprocess.py runs without the SDK
    from cirro.helpers.preprocess_dataset import PreprocessDataset

    ds = PreprocessDataset.from_running()

    custom_reference = staged_reference(ds.params, ds.references_base)
    if custom_reference is not None:
        ds.logger.info(f"Using staged reference for {ds.params['dada_ref_taxonomy']}")
        for name, value in custom_reference.items():
            ds.add_param(name, value, overwrite=True)

        # Dead once a custom reference is set, and leaving it in overstates what ran
        ds.remove_param("dada_ref_taxonomy", force=True)

    # log
    ds.logger.info(ds.params)

    # Force params.json to be written: the HealthOmics pre-process Lambda fails the run
    # when the file is absent, and the SDK writes it only when a parameter changes.
    ds.keep_params(list(ds.params.keys()))
