
from cirro.helpers.preprocess_dataset import PreprocessDataset

# Instantiate the Cirro dataset object
ds = PreprocessDataset.from_running()

# Log the parameters present
for k, v in ds.params.items():
    ds.logger.info(f"{k}: {v}")

# Force params.json to be written: the HealthOmics pre-process Lambda fails the run
# when the file is absent, and the SDK writes it only when a parameter changes.
ds.keep_params(list(ds.params.keys()))
