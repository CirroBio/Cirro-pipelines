# Convert Image to OME-TIFF

Converts whole-slide and microscopy images to pyramidal OME-TIFF using QuPath.

## Notes

Every file under the input dataset carrying the selected suffix is converted, at any
depth, so a dataset which keeps its slides in subfolders works as well as a flat one.

Each image is converted in its own task. They run in parallel, stage only themselves
rather than the whole input dataset, and publish as they finish — so an interrupted run
keeps everything converted up to that point and a resume takes the rest.

Formats which pack several images into one file — `vsi`, `scn`, `dcm`, `dicom`, `nd2`,
`lif`, `czi` — write one OME-TIFF per series. The VSI format in particular is a pointer
to other files (with the `.ets` extension) which hold the image data for each series.
Everything else writes a single OME-TIFF per image.

Alongside the images the run writes `<image>.log.txt` per image and one
`conversion_timings.tsv`, giving seconds and input/output bytes per image.

## Compression

JPEG is the default. Images from brightfield slide scanners already hold lossy JPEG
tiles written by the scanner, so a lossless codec re-encodes them and inflates the
output without recovering any fidelity that was lost before the file was written.
Measured on one Aperio SVS, output size relative to JPEG:

| Compression | vs JPEG |
|---|---|
| JPEG | 1.0x |
| DEFAULT | 1.0x (it selects JPEG for RGB) |
| J2K_LOSSY | 5.2x |
| J2K | 7.1x |
| ZLIB | 8.9x |
| LZW | 12.2x |
| UNCOMPRESSED | 36.5x |

Choose a lossless option for acquisitions that were never lossily compressed, such as
fluorescence. `JPEG` requires an RGB or 8-bit image and `J2K` requires 8- or 16-bit;
QuPath refuses the conversion rather than writing something wrong.
