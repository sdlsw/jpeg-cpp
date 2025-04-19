# jpeg-cpp

An implementation of JPEG image compression written in C++23. This is a learning project, but doubles as a code sample for portfolio purposes. I had a few goals in mind while writing this:

- Learn modern C++ using the most recent standard with reasonable support. As of writing this is C++23.
- Demonstrate ability to use technology I'm less familiar with to implement something nontrivial.
- Demonstrate ability to implement something described in a standards document.

I picked JPEG compression for the subject, since implementing it takes some work, it's standard, and most people are at least familiar with its existence. Plus, I'm fond of digital image algorithms in general.

In addition, I gave myself some rules to follow. First, no third party libraries. The goal was to familiarize myself with new language features (including stdlib), and having to write everything from scratch helped with that. Second, I did not refer to other JPEG implementations until I'd gotten mine working. I wanted to solve as many problems myself as possible.

## Building

`jpeg-cpp` uses [CMake](https://cmake.org/) (version 3.30 or higher required). Currently it has only been tested with MSVC. To build:

```sh
git clone https://github.com/sdlsw/jpeg-cpp.git
cd jpeg-cpp
mkdir cd build
cmake ..
cmake --build . --config Release
```

## Usage

The executable, `jpeg`, may be used to compress and decompress JPEG images. The only accepted uncompressed format is 24-bit color BMP files. The following sections are a brief description of the most common uses of the tool. To see a complete usage guide, run `jpeg --help`.

### Basic Encode
Compresses some raw image data into a JPEG file.

```sh
jpeg encode example.bmp -o example.jpeg
```

### Basic Decode
Decompresses a JPEG file to raw image data.

```sh
jpeg decode example.jpeg -o example.bmp
```

### File Scan
Reads a JPEG file and dumps all of its metadata to stdout. Useful for debugging.

```sh
jpeg scan example.jpeg
```

### Encode with Quality
`jpeg-cpp` supports a quality setting, from 0 to 100. Lower values yield smaller JPEG size at the cost of visual image quality. The default setting is 80.

```sh
jpeg encode example.bmp -o example.jpeg -q 50
```

## Examples
See [encode-examples](encode-examples) for some sample JPEGs produced by the encoder.

## Known Bugs/Issues

- JPEG encoding appears to slightly alter the colors of the image, to the point that there is a noticeable visual difference. When compared to other encoders (such as the one used in GIMP), this difference does not occur. It's still not known why this is the case.
- Progressive JPEG is unsupported. This will be fixed in a pending update.
