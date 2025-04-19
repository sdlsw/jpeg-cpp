# Encoding Examples
This folder contains sample JPEG files produced by the encoder at various quality levels.

# Files

## `uncompressed-original.png`
The original image.

## `uncompressed.bmp`
The original image, converted to BMP using ffmpeg. The BMP is also scaled to 75% of the original dimensions, so it fits into GitHub's file size limits.

## `qualitytestNNN.jpg`
Outputs from invoking the following command:

```sh
jpeg qualitytest uncompressed.bmp
```

These JPEG files range in quality level from 0 - 100, in increments of 10.

## `quality50-artifacts.png`
This image is a highly zoomed in (~600%) piece of `qualitytest050.jpg`, showing the 8x8 block artifacts characteristic of lower quality JPEGs.

# Sources
## `uncompressed-original.png`
- Title: Wintery landscape
- Author: Friedrich Mook
- Year: 1926
- Copyright: Public Domain
- Link: https://staedelmuseum.de/go/ds/sg1093z