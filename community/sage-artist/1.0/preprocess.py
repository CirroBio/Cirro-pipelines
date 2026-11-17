#!/usr/bin/env python3

import pandas as pd


def make_samplesheet(images, he_images) -> pd.DataFrame:
    """The samplesheet nf-artist reads from its `input` param.

    Its columns are:
    image: [string] Path or URI to image to be processed
    convert: [boolean] Should the image be converted to a OME-TIFF
    he: [boolean] Is the image a H&E image
    minerva: [boolean] Should a Minerva story be generated
    miniature: [boolean] Should a Miniature thumbnail be generated
    id: optional [string] A custom identifier to replace image simpleName in output directory structure
    """

    samplesheet = pd.DataFrame([
        {
            "image": img,
            "convert": True,
            "he": he,
            "minerva": True,
            "miniature": True
        }
        for img_list, he in [
            (images or [], False),
            (he_images or [], True)
        ]
        for img in img_list
        if img
    ])

    if samplesheet.empty:
        raise ValueError("No images selected -- stopping")
    return samplesheet


if __name__ == "__main__":

    # Imported here so that test_preprocess.py runs without the SDK
    from cirro.helpers.preprocess_dataset import PreprocessDataset

    ds = PreprocessDataset.from_running()

    ds.logger.info("User-provided (non-H&E) images:")
    ds.logger.info(ds.params.get("images", ""))

    ds.logger.info("User-provided H&E images:")
    ds.logger.info(ds.params.get("he_images", ""))

    samplesheet = make_samplesheet(ds.params.get("images"), ds.params.get("he_images"))

    ds.logger.info("Samplesheet:")
    ds.logger.info(samplesheet.to_csv(index=None))
    # Write to the dataset's config/ folder (mapped in process-input.json)
    ds.logger.info(f"Writing samplesheet to {ds.params['input']}")
    samplesheet.to_csv(ds.params["input"], index=False)
