#!/usr/bin/env python3
from collections import defaultdict
import pandas as pd
from cirro.helpers.preprocess_dataset import PreprocessDataset


def main() -> None:
    ds = PreprocessDataset.from_running()
    samplesheet = parse_samplesheet(ds)
    # Write to the dataset's config/ folder (mapped in process-input.json)
    samplesheet.to_csv(ds.params["input"], sep="\t", index=None)


def parse_samplesheet(ds: PreprocessDataset) -> pd.DataFrame:

    # Fix an odd edge case where the metadata attribute 'file' has been
    # added to the sample table
    if 'file' in ds.samplesheet.columns.values:
        ds.samplesheet = ds.samplesheet.drop(
            columns=['file']
        )

    # Merge the file table with the sample metadata
    samplesheet = (
        ds.files
        .merge(
            ds.samplesheet,
            left_on="sample",
            right_on="sample"
        )
    )
    log_table("Starting table", ds, samplesheet)

    # If a subset of files were selected
    if len(ds.params.get("subset_files", "")) > 0:
        samplesheet = samplesheet.loc[
            samplesheet['file'].isin(ds.params["subset_files"].split(","))
        ]
        log_table("Subset of files selected", ds, samplesheet)
    else:
        ds.logger.info("All files will be processed - no subset selected")

    # Drop sampleIndex, dataset, and process
    for cname in ["sampleIndex", "dataset", "process"]:
        if cname in samplesheet.columns.values:
            samplesheet = samplesheet.drop(columns=[cname])

    # Map the required values
    samplesheet = (
        samplesheet
        .assign(**{
            "source name": samplesheet["sample"],
            "comment[file uri]": samplesheet["file"],
            "comment[data file]": samplesheet["file"].apply(
                lambda s: s.split("/")[-1]
            )
        })
        .drop(
            columns=['sample', 'file']
        )
    )
    # If the assay name is not provided, fall back on the source name
    if "assay name" not in samplesheet.columns.values:
        samplesheet = samplesheet.assign(**{
            "assay name": samplesheet["source name"]
        })

    # Drop any columns which are entirely missing
    missing = [
        cname
        for cname, cvals in samplesheet.items()
        if cvals.apply(pd.isnull).all()
    ]
    if len(missing) > 0:
        if len(missing) == samplesheet.shape[1]:
            raise Exception("All columns appear to be lacking data")
        ds.logger.info(
            f"Dropping columns which lack data: {', '.join(missing)}"
        )
        samplesheet = samplesheet.drop(columns=missing)

    log_table("Formatted source columns", ds, samplesheet)

    # If a user forgets the "comment[]" part of a name, fix it
    renamed_any = False
    for group, vals in [
        [
            "characteristics", [
                "age",
                "ancestry category",
                "biological replicate",
                "cell line",
                "cell type",
                "developmental stage",
                "disease",
                "individual",
                "organism part",
                "organism",
                "sex",
            ]
        ],
        [
            "comment", [
                "cleavage agent details",
                "data file",
                "file uri",
                "flow rate chromatogram",
                "fraction identifier",
                "fractionation method",
                "fragment mass tolerance",
                "gradient time",
                "instrument",
                "label",
                "modification parameters",
                "precursor mass tolerance",
                "technical replicate",
            ]
        ],
        [
            "factor value", [
                "flow rate chromatogram",
                "gradient time"
            ]
        ]
    ]:
        # e.g. 'age'
        for val in vals:
            # e.g. 'characteristics[age]'
            comb = f"{group}[{val}]"

            # If the short name is present, but not the long one
            if (
                val in samplesheet.columns.values and
                comb not in samplesheet.columns.values
            ):
                # Fix the name
                ds.logger.info(f"Renaming column {val} -> {comb}")
                samplesheet = samplesheet.rename(
                    columns={
                        val: comb
                    }
                )
                renamed_any = True

    if renamed_any:
        log_table("Renamed columns", ds, samplesheet)

    # Populate default values for required columns
    populated_any = False
    for cname, default in [
        ("characteristics[organism]", ds.params["organism"]),
        ("characteristics[organism part]", "not applicable"),
        ("characteristics[disease]", "not applicable"),
        ("characteristics[cell type]", "not applicable"),
        ("characteristics[biological replicate]", "not applicable"),
        ("comment[cleavage agent details]", "NT=Trypsin"),
        ("comment[instrument]", "NT=LTQ Orbitrap XL"),
        ("comment[label]", "NT=label free sample"),
        ("comment[fraction identifier]", "1")
    ]:
        if cname not in samplesheet.columns.values:
            ds.logger.info(f"Added {cname} = {default}")
            samplesheet = samplesheet.assign(
                **{cname: default}
            )
            populated_any = True

    for kw in ["organism"]:
        if kw in ds.params:
            ds.remove_param(kw)

    if populated_any:
        log_table("Populated default values", ds, samplesheet)

    # Populate the technical replicate
    if "comment[technical replicate]" not in samplesheet.columns.values:
        counter = defaultdict(int)

        replicate = []
        for sample in samplesheet["source name"].values:
            counter[sample] += 1
            replicate.append(int(counter[sample]))
        samplesheet = samplesheet.assign(**{
            "comment[technical replicate]": replicate
        })
        log_table("Populated technical replicate", ds, samplesheet)

    # Add the mass tolerances
    for tol in ["fragment", "precursor"]:
        kw = f"comment[{tol} mass tolerance]"
        val = ds.params[f"{tol}_mass_tolerance"]
        unit = ds.params[f"{tol}_mass_tolerance_unit"]
        ds.logger.info(f"Setting {kw} = {val} {unit}")
        samplesheet = samplesheet.assign(**{
            kw: f"{val} {unit}"
        })

    # Sort the column names
    samplesheet = samplesheet.reindex(
        columns=sort_columns(samplesheet.columns.values)
    )
    log_table("Sorted columns", ds, samplesheet)

    ds.logger.info("Done preparing samplesheet")

    return samplesheet


def sort_columns(columns):
    sorted = (
        ["source name"] +
        [cname for cname in columns if cname.startswith("characteristics")] +
        ["assay name"] +
        [cname for cname in columns if cname.startswith("comment")]
    )
    return sorted + [cname for cname in columns if cname not in sorted]


def log_table(msg: str, ds: PreprocessDataset, df: pd.DataFrame):
    ds.logger.info(msg)
    for line in df.to_csv(index=None).split("\n"):
        ds.logger.info(line)


if __name__ == "__main__":
    main()
