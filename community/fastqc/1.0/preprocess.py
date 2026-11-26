#!/usr/bin/env python3


def input_folders(value) -> str:
    """The picked folders as the workflow's comma-separated --input.

    The workflow appends its glob (`**.fastq.gz`) straight onto each folder, so every
    folder must end in exactly one slash or the glob would also match sibling folders
    that share its name as a prefix.
    """
    if not value:
        raise ValueError("No input folders were provided: pick at least one dataset or folder")
    return ",".join(folder.rstrip("/") + "/" for folder in value.split(","))


if __name__ == "__main__":

    # Imported here so that test_preprocess.py runs without the SDK
    from cirro.helpers.preprocess_dataset import PreprocessDataset

    ds = PreprocessDataset.from_running()

    ds.add_param("input", input_folders(ds.params.get("input")), overwrite=True)
