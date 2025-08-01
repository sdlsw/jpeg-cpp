# Encoding Examples
This folder contains sample JPEG files produced by the encoder at various quality levels.

# Files

## `uncompressed-original.png`
The original image. Note that this image was taken by a phone camera that does not offer access to RAW copies. The image was instead taken as a JPEG at 100% quality, then converted to PNG. Not perfect, but enough to showcase the algorithm.

## `uncompressed.bmp`
The original image, converted to BMP using ffmpeg.

## `qualitytestNNN.jpg`
Outputs from invoking the following command:

```sh
jpeg qualitytest uncompressed.bmp
```

These JPEG files range in quality level from 0 - 100, in increments of 10.

## `quality050-artifacts.png`
This image is a highly zoomed in (~600%) piece of `qualitytest050.jpg`, showing the 8x8 block artifacts characteristic of lower quality JPEGs.

# Sources
## `uncompressed-original.png`
- Own work.
- License: [CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/)