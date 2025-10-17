# Convert VSI to OME-TIFF

Converts images in VSI format to the OME-TIFF format using QuPath.

## Notes

The VSI file format is actually a pointer to other files (with the .ets extension)
which contain image data for one or more series within the image.
This utility will create one OME-TIFF file for each of the series contained
within each .vsi file present at the top level of the input dataset.