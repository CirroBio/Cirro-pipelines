from cirro.helpers.preprocess_dataset import PreprocessDataset
import pandas as pd

ds = PreprocessDataset.from_running()

# Samplesheet requirements
# The samplesheet specified in the input parameter should be a CSV file with the following columns

# image: [string] Path or URI to image to be processed
# convert: [boolean] Should the image be converted to a OME-TIFF
# he: [boolean] Is the image a H&E image
# minerva: [boolean] Should a Minerva story be generated
# miniature: [boolean] Should a Miniature thumbnail be generated
# id: optional [string] A custom identifier to replace image simpleName in output directory structure

ds.logger.info("User-provided (non-H&E) images:")
ds.logger.info(ds.params.get("images", ""))

ds.logger.info("User-provided H&E images:")
ds.logger.info(ds.params.get("he_images", ""))

samplesheet = pd.DataFrame([
    {
        "image": img,
        "convert": True,
        "he": he,
        "minerva": True,
        "miniature": True
    }
    for img_list, he in [
        (ds.params.get("images", "").split(","), False),
        (ds.params.get("he_images", "").split(","), True)
    ]
    for img in img_list
    if img
])

assert samplesheet.shape[0] > 0, "No images selected -- stopping"

ds.logger.info("Samplesheet:")
ds.logger.info(samplesheet.to_csv(index=None))
ds.logger.info("Writing samplesheet to input.csv")
samplesheet.to_csv("input.csv", index=False)
