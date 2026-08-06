#!/usr/bin/env python3
"""Log the parameters the run was launched with. Makes no changes to them."""

import json

from cirro.helpers.preprocess_dataset import PreprocessDataset


def main() -> None:
    ds = PreprocessDataset.from_running()
    ds.logger.info("Params:\n" + json.dumps(ds.params, indent=4, default=str))


if __name__ == "__main__":
    main()
