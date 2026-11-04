#!/usr/bin/env python3

import re

import pandas as pd


def make_manifest(params: dict, dataset_name: str) -> pd.DataFrame:
    """
    nf-core/demux input requires this format:

    id,samplesheet,lane,flowcell
    DDMMYY_SERIAL_NUMBER_FC,/path/to/SampleSheet.csv,1,/path/to/sequencer/output
    DDMMYY_SERIAL_NUMBER_FC,/path/to/SampleSheet.csv,2,/path/to/sequencer/output

    The demux pipeline can accept multiple run directories, but for the purposes of
    Cirro, we will limit it to one run directory / flow cell.

    Same thing with the flow cell sample sheet. Generally the FC samplesheet is used
    to separate and route dataset, we assume that it will just be
    a separate analysis job / dataset as well.
    """

    # The ID reaches the demultiplexer's command line (and is Cell Ranger's --id), so
    # keep only the characters Cell Ranger accepts
    flowcell_id = re.sub(r'[^A-Za-z0-9_-]', '', params.get('flowcell_id') or dataset_name)
    if not flowcell_id:
        raise ValueError("Please provide a Flow Cell ID using letters, numbers, '-' or '_'")
    lanes = (params.get('lane') or '').split(',')
    samplesheet_path = params.get('samplesheet')

    run_folder = params.get('run_folder')
    if not run_folder:
        raise ValueError("Please select the sequencing run folder")

    manifest = pd.DataFrame.from_records([
        {
            'id': flowcell_id,
            'samplesheet': samplesheet_path,
            'lane': lane,
            'flowcell': run_folder
        }
        for lane in lanes
    ])
    return manifest


if __name__ == '__main__':

    # Imported here so that test_preprocess.py runs without the SDK
    from cirro.helpers.preprocess_dataset import PreprocessDataset

    ds = PreprocessDataset.from_running()
    samplesheet = make_manifest(ds.params, ds.metadata['dataset']['name'])
    ds.logger.info(samplesheet.to_csv())
    # Write to the dataset's config/ folder (mapped in process-input.json)
    samplesheet.to_csv(ds.params["input"], index=False)
