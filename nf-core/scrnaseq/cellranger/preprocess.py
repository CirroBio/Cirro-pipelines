import pandas as pd

from cirro.helpers.preprocess_dataset import PreprocessDataset

def make_samplesheet(ds: PreprocessDataset) -> pd.DataFrame:
    """
    Create a samplesheet for my workflow
    """
    samplesheet = ds.pivot_samplesheet(
        metadata_columns=[],
        file_filter_predicate='readType == "R"'
    )
    return samplesheet


if __name__ == '__main__':
    ds = PreprocessDataset.from_running()
    samplesheet = make_samplesheet(ds)
    samplesheet.to_csv('samplesheet.csv', index=False)
