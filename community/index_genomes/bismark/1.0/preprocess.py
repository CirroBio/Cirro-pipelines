#!/usr/bin/env python3

import re
from cirro.helpers.preprocess_dataset import PreprocessDataset

# Allowlist: flags, values, and paths used in legitimate CLI arguments.
# Blocks shell metacharacters: ; | & > < ` $ ( ) { } ! ~ * ? # ^ " ' \
_SAFE_ARGS_RE = re.compile(r'^[a-zA-Z0-9 \t\-=.,:/+@%_]*$')


def validate_extra_args(ds: PreprocessDataset, param: str):
    value = ds.params.get(param)
    if not value:
        return
    if not _SAFE_ARGS_RE.match(value):
        raise ValueError(
            f"Invalid characters detected in '{param}'. "
            "Only alphanumeric characters, spaces, hyphens, underscores, "
            "dots, equals signs, commas, colons, slashes, and plus signs are allowed."
        )
    ds.logger.info(f"Validated {param}: {value!r}")


if __name__ == "__main__":
    ds = PreprocessDataset.from_running()
    validate_extra_args(ds, "bismark_extra_args")
