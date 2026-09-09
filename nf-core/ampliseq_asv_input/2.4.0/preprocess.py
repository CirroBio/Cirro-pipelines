#!/usr/bin/env python3

from cirro.helpers.preprocess_dataset import PreprocessDataset


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

    ds = PreprocessDataset.from_running()

    use_staged_reference(ds)

    # log
    ds.logger.info(ds.params)
