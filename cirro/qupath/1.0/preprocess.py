#!/usr/bin/env python3


def empty_params(params: dict) -> list:
    """Names of the params the user left empty.

    args is optional, and a user who supplies none leaves an empty string behind.
    That is not a path and not alphanumeric, so HealthOmics rejects it; the
    workflow's own default applies once it is gone. A form has to be able to
    express "not set", so this is the point where it can be guaranteed.
    """
    return [
        param for param, value in params.items()
        if isinstance(value, str) and not value.strip()
    ]


if __name__ == "__main__":

    # Imported here so that test_preprocess.py runs without the SDK
    from cirro.helpers.preprocess_dataset import PreprocessDataset

    ds = PreprocessDataset.from_running()

    for param in empty_params(ds.params):
        ds.logger.info(f"Dropping empty parameter: {param}")
        ds.remove_param(param, force=True)

    ds.logger.info("Parameters for this run:")
    for param, value in ds.params.items():
        ds.logger.info(f"{param}: {value}")

    # Force params.json to be written: the HealthOmics pre-process Lambda fails the
    # run when the file is absent, and the SDK writes it only when a parameter changes.
    ds.keep_params(list(ds.params.keys()))
