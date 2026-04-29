#!/usr/bin/env python3

from cirro.helpers.preprocess_dataset import PreprocessDataset
import pandas as pd
import re


def format_sample_id(s: str) -> str:
    s = str(s)
    s = re.sub('[^a-zA-Z0-9_]+', '_', s)
    s = s if re.match("^[a-zA-Z]+.*", s) else f"Sample_{s}"
    return s


def make_manifest(ds: PreprocessDataset) -> pd.DataFrame:

    ds.logger.info("Files annotated in the dataset:")
    ds.logger.info(ds.files)

    # Filter out any index files that may have been uploaded
    ds.files = ds.files.loc[
        ds.files.apply(
            lambda r: r.get('readType', 'R') == 'R',
            axis=1
        )
    ]

    # Make a wide samplesheet with the columns
    # sampleID, forwardReads, reverseReads
    samplesheet = pd.concat([
        pd.DataFrame(dict(
            sample=sample,
            **{
                f"R{read}": read_files['file'].sort_values().tolist()
                for read, read_files in sample_files.groupby("read")
            }
        ))
        for sample, sample_files in ds.files.assign(
            read=ds.files["read"].apply(int)
        ).query(
            "read <= 2"
        ).groupby("sample")
    ]).rename(
        columns={
            "sample": "sampleID",
            "R1": "forwardReads",
            "R2": "reverseReads"
        }
    ).reindex(
        columns=[
            "sampleID",
            "forwardReads",
            "reverseReads"
        ]
    ).fillna(
        ""
    )

    # Replace all disallowed characters
    # Make sure that every sampleID starts with a letter
    samplesheet["sampleID"] = (
        samplesheet["sampleID"]
        .apply(format_sample_id)
    )

    # Make sure that every sample name is unique
    samplename_vc = samplesheet["sampleID"].value_counts()
    if (samplename_vc > 1).any():
        ds.logger.info("Duplicate sample names detected")
        for samplename, count in samplename_vc.items():
            if count > 1:
                ds.logger.info(f"{samplename}: {count}")
        raise ValueError("Duplicate sample names detected")

    ds.logger.info("Formatted samplesheet:")
    ds.logger.info(samplesheet.to_csv(index=None))

    # If the user specified that only the reverse reads should be used
    if ds.params.get("only_use_forward_reads", False):
        ds.logger.info("User has requested that only forward reads are used")
        samplesheet = samplesheet.drop(columns=["reverseReads"])
        ds.add_param("single_end", True)

    # If any of the samples lack reverse reads
    elif (samplesheet["reverseReads"] == "").any():
        ds.logger.info("Dataset appears to be single-ended")
        samplesheet = samplesheet.drop(columns=["reverseReads"])
        ds.add_param("single_end", True)

    if 'only_use_forward_reads' in ds.params:
        ds.remove_param("only_use_forward_reads", force=False)

    msg = "No files detected -- there may be an error with data ingest"
    assert samplesheet.shape[0] > 0, msg

    return samplesheet


if __name__ == "__main__":

    ds = PreprocessDataset.from_running()

    samplesheet = make_manifest(ds)

    # Write out to a file
    samplesheet.to_csv("samplesheet.tsv", index=None, sep="\t")

    # Add the param for the samplesheet
    ds.add_param("input", "samplesheet.tsv")

    # If the user set ignore_empty_input_files
    if ds.params.get("ignore_empty_input_files", False):
        # Also set the other ignore flags
        ds.add_param("ignore_failed_filtering", True, overwrite=True)
        ds.add_param("ignore_failed_trimming", True, overwrite=True)

    # Write out the table of sample metadata as well
    metadata = ds.samplesheet.rename(
        columns=dict(
            sample="ID"
        )
    )

    if metadata.shape[1] > 1:

        metadata = metadata.reindex(
            columns=["ID"] + [
                cname
                for cname in metadata.columns.values
                if cname != "ID"
            ]
        )

        metadata["ID"] = metadata["ID"].apply(format_sample_id)

        ds.logger.info("Formatted metadata:")
        ds.logger.info(metadata)

        metadata.to_csv(
            "metadata.tsv",
            sep="\t",
            index=None
        )

        # Add the param for the metadata
        ds.add_param("metadata", "metadata.tsv")

    # If a custom reference was provided
    if ds.params.get("dada_ref_tax_custom", False):
        # If the species-level reference was not provided
        if not ds.params.get("dada_ref_tax_custom_sp", False):
            # Turn off the species-level annotation
            ds.logger.info("Turning off species-level annotation, since no species reference was provided")
            ds.add_param("skip_dada_addspecies", True)

    # log
    ds.logger.info(ds.params)
