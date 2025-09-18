import pandas as pd

from cirro.helpers.preprocess_dataset import PreprocessDataset

def make_manifest(ds: PreprocessDataset) -> pd.DataFrame:
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

    flowcell_id = ds.params.get('flowcell_id') or ds.metadata['dataset']['name']
    lane = ds.params.get('lane')
    samplesheet_path = ds.params.get('samplesheet')

    # Try to find the correct run directory dataset
    # Since Samplesheet can be provided in a secondary dataset
    run_dir = next((
        dataset['dataPath'] for dataset in ds.metadata['inputs']
        if dataset['processId'] != 'files'
    ), None)

    if not run_dir:
        raise ValueError("Please provide at least one dataset with the sequencing run type")

    # TODO: support multiple lanes

    manifest = pd.DataFrame.from_records([
        {
            'id': flowcell_id,
            'samplesheet': samplesheet_path,
            'lane': lane,
            'flowcell': run_dir
        }
    ])
    return manifest


if __name__ == '__main__':
    ds = PreprocessDataset.from_running()
    samplesheet = make_manifest(ds)
    ds.logger.info(samplesheet.to_csv())
    samplesheet.to_csv('manifest.csv', index=False)
