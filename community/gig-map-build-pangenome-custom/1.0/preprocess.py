
from cirro.helpers.preprocess_dataset import PreprocessDataset

# Instantiate the Cirro dataset object
ds = PreprocessDataset.from_running()

# Log the parameters present
for k, v in ds.params.items():
    ds.logger.info(f"{k}: {v}")
