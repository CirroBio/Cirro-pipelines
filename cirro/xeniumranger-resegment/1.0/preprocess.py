from cirro.helpers.preprocess_dataset import PreprocessDataset

ds = PreprocessDataset.from_running()

# Make sure that the user provided an input path
input_uri = ds.params.get("input", "")
msg = "Error: Must provide input dataset"
assert len(input_uri) > 0, msg

# Remove the '/analysis.zarr.zip' suffix
suffix = "/analysis.zarr.zip"
msg = f"Error: Expected input to end with '{suffix}', not {input_uri}"
assert input_uri.endswith(suffix), msg
ds.add_param(
    "input",
    input_uri[:-len(suffix)],
    overwrite=True
)

# If the user selected a modified expansion distance, then
# the resegment-nuclei parameter must also be set to True
if ds.params.get("expansion_distance") is not None:
    if ds.params.get("resegment_nuclei", False) is False:
        msg = "Setting --resegment-nuclei True (expansion-distance was set)"
        ds.logger.info(msg)
        ds.add_param("resegment_nuclei", True, overwrite=True)

# Consumed above; not a parameter the workflow declares.
ds.remove_param("resegment_nuclei", force=True)
