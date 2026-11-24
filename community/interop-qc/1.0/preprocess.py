#!/usr/bin/env python3


def run_folder(value) -> str:
    """The picked run folder, without the trailing slash a whole-dataset pick has.

    The workflow joins `/*xml` and `/InterOp/*bin` onto it, and would otherwise glob a
    path containing `//`.
    """
    if not value:
        raise ValueError("No sequencing run folder was selected")
    return value.rstrip("/")


if __name__ == "__main__":

    # Imported here so that test_preprocess.py runs without the SDK
    from cirro.helpers.preprocess_dataset import PreprocessDataset

    ds = PreprocessDataset.from_running()

    ds.add_param("runDir", run_folder(ds.params.get("runDir")), overwrite=True)

    # Force params.json to be written: the HealthOmics pre-process Lambda fails the run
    # when the file is absent, and the SDK writes it only when a parameter changes.
    ds.keep_params(list(ds.params.keys()))
