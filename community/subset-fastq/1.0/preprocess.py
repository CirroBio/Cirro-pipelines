#!/usr/bin/env python3

from cirro.helpers.preprocess_dataset import PreprocessDataset


def process_inputs(ds: PreprocessDataset):

    # Get the total list of files available
    manifest = ds.files

    ds.logger.info(manifest.to_csv(index=None))

    # Populate the input params appropriately
    ds.add_param(
        'input',
        ','.join(manifest['file'].tolist())
    )


if __name__ == '__main__':
    ds = PreprocessDataset.from_running()

    process_inputs(ds)
