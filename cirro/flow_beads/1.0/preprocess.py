#!/usr/bin/env python3

from cirro.helpers.preprocess_dataset import PreprocessDataset

ds = PreprocessDataset.from_running()


def _deduplicate_fcs(fcs_files: list) -> str:
    """
    Deduplicate a list of FCS files,
    making sure that the file names are all unique
    (even if some files are in different folders).

    :param fcs_files: List of FCS files.
    :return: Deduplicated comma-separated list of FCS files.
    """
    # This pattern is being used to differentiate the URI from the filename
    seen = set([])
    output = []

    for file in fcs_files:
        fn = file.split("/")[-1]
        if fn in seen:
            ds.logger.warning(f"Removing duplicate file: {file}")
        else:
            seen.add(fn)
            output.append(file)

    return ",".join(output)


if len(ds.params.get("subset") or []) > 0:
    n = len(ds.params['subset'])
    ds.logger.info(f"User specified a subset of {n:,} files to analyze")
    ds.add_param(
        "input_fcs",
        _deduplicate_fcs(ds.params["subset"]),
        overwrite=True
    )
else:
    ds.logger.info("No subset specified -- analyzing all files")

    input_fcs = [
        file
        for file in ds.files["file"].tolist()
        # Instruments write the extension in either case (e.g. Cytomics FC 500: .FCS)
        if file.lower().endswith(".fcs")
    ]
    ds.logger.info(f"Found {len(input_fcs):,} FCS files to analyze")
    ds.add_param(
        "input_fcs",
        _deduplicate_fcs(input_fcs)
    )

    input_anndata = [
        file
        for file in ds.files["file"].tolist()
        if file.endswith(".h5ad")
    ]
    if len(input_anndata) > 0:
        ds.logger.info(f"Found {len(input_anndata):,} anndata files to analyze")
        ds.add_param(
            "input_anndata",
            ",".join(input_anndata)
        )
    else:
        ds.logger.info("No anndata files found")

# log
ds.logger.info(ds.params)

# Force params.json to be written: the HealthOmics pre-process Lambda fails the run
# when the file is absent, and the SDK writes it only when a parameter changes.
ds.keep_params(list(ds.params.keys()))
