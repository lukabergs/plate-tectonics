Plate-tectonics
===============

[![Build Status](https://github.com/Mindwerks/plate-tectonics/actions/workflows/push.yml/badge.svg)](https://github.com/Mindwerks/plate-tectonics/actions/workflows/push.yml)

This is a library to simulate plate tectonics.
It is written in C++ and it has Python bindings (as part of this project), as well as Haskell bindings ([hplatec](http://github.com/ftomassetti/hplatec)).

How can I use it?
=================

Being a library you want probably to use it inside some larger program. From example [WorldEngine](https://github.com/Mindwerks/worldengine) (a world generator) is based on plate-tectonics.

You can also use the examples to just run the code of this library and generate a few maps. However the examples do not unleash the full power of this library. For running the examples check section _Running the examples (C++)_.

How it looks like
=================

The library offers an API to generate heightmaps and some other data about the world resulting from the simulation. The example permits also to generate maps like this one:

![](https://raw.githubusercontent.com/Mindwerks/plate-tectonics/master/screenshots/map_grayscale.png)

![](https://raw.githubusercontent.com/Mindwerks/plate-tectonics/master/screenshots/map_colors.png)

You can see a video of simulation based on an old version of this library: http://www.youtube.com/watch?v=bi4b45tMEPE#t=0

How to build plate-tectonics (C++)
==================================

We use [CMake](http://www.cmake.org/). Install it and then run the folowing commands

### Linux

```
mkdir -p build
cd build
cmake .. -G "Unix Makefiles"
make
```

### Mac OS-X

```
mkdir -p build
cd build
cmake ..
make
```

This should produce a library (libPlateTectonics.a) in the build directory.

### Windows

```
mkdir build
cd build
cmake .. -DCMAKE_TOOLCHAIN_FILE=C:/dev/vcpkg/scripts/buildsystems/vcpkg.cmake -DMSVC_RUNTIME=dynamic
cmake --build .
```

This project uses `libpng` for tests and examples. With the vcpkg toolchain enabled, the root `vcpkg.json` manifest installs that dependency automatically.

If you want to build also the examples run:

```
mkdir -p build
cd build
cmake .. -DWITH_EXAMPLES=ON
make
```

Note: All builds are now done in the `build/` directory to keep the source tree clean. The build directory is excluded from version control via `.gitignore`.

To compile on other platforms please run:

```
cmake --help
```

Running the examples (C++)
==========================

To run also the examples you need to install the library libpng.

From the root directory run:

```bash
mkdir -p build
cd build
cmake .. -DWITH_EXAMPLES=ON -DCMAKE_TOOLCHAIN_FILE=C:/dev/vcpkg/scripts/buildsystems/vcpkg.cmake -DMSVC_RUNTIME=dynamic
make
cd examples
./simulation --dim 600 400 --filename world
```

Metric heightmaps
=================

The example CLI now separates authoritative terrain data from preview rendering.

- Metric terrain domain: `0..65535`
- `0`: deepest seafloor sample
- `65535`: highest terrain sample
- `--sea-level-m`: coastline/crust-classification threshold in that metric domain
- Default data exports: `<base>.r16`, `<base>.r16.json`, `<base>.png` (16-bit grayscale PNG)
- Preview export: `<base>_preview.png`

Example workflows:

```bash
# Procedural world with explicit metric sea level
./simulation --dim 600 400 --sea-level-m 32000 --filename world

# Import a metric PNG16 heightmap
./simulation --input-png16 world.png

# Import a raw metric .r16 heightmap
./simulation --input-r16 world.r16
```

`.r16` files are little-endian `uint16_t`, row-major, one sample per pixel. Sidecar metadata lives in `<base>.r16.json` and stores width, height, `sea_level_m`, format, endianness, layout, and version.

How to run tests (C++)
======================

GoogleTest is automatically fetched by CMake using FetchContent, so no manual installation is required.

After building the library with CMake in the build directory:

```bash
cd build/test
make
./PlateTectonicsTests
```

Currently the test coverage is still poor (but improving!), tests are present only for new code and tiny portion of the old code that were refactored.

## Python bindings

Supported versions:
* Python >= 3.9

### Quick Start (Python)

Install using pip:

```bash
pip install PyPlatec
```

### Usage (Python)

The library is quite simple. `platec.create()` requires all 10 parameters and supports both positional and keyword arguments:

```python    
    import platec

    # Using positional arguments
    p = platec.create(3, 512, 512, 0.65, 60, 0.02, 1000000, 0.33, 2, 10)
    
    while platec.is_finished(p) == 0:
        platec.step(p)
    
    hm = platec.get_heightmap(p)
    platec.destroy(p)
```

The Python bindings also expose metric helpers on an existing simulation:

```python
platec.load_heightmap_u16(p, metric_samples, 32000)
sea_level_m = platec.get_sea_level_m(p)
```


With keyword arguments for clarity:

```python
    import platec
    
    # Using keyword arguments (recommended for readability)
    p = platec.create(
        seed=3,
        width=1000,
        height=800,
        sea_level=0.65,
        erosion_period=60,
        folding_ratio=0.02,
        aggr_overlap_abs=1000000,
        aggr_overlap_rel=0.33,
        cycle_count=2,
        num_plates=10
    )
    
    while platec.is_finished(p) == 0:
        platec.step(p)
    
    hm = platec.get_heightmap(p)
    platec.destroy(p)
```

### Building from Source (Python)

If you need to build from source instead of using pre-built wheels:

```bash
cd pybindings
python setup.py build
python setup.py install
```

For development:
```bash
pip install -e .
```

Code Quality and Linting
=========================

This project uses C++ linters to maintain code quality. See [LINTING.md](LINTING.md) for details.

Quick start:
```bash
# Run all linters
./run_linter.sh

# Run specific linter
./run_linter.sh clang-tidy
./run_linter.sh cppcheck

# Format code
./run_astyle.sh
```

Plans for the future
====================

* Improve the quality of the code and add some tests
* Support Google protocol buffer

Projects using plate-tectonics
==============================

[WorldEngine](http://github.com/Mindwerks/worldengine), a world generator  
[Widelands](https://www.widelands.org/maps/lost-islands/), a free, open source real-time strategy game  

Original project
================

A fork of platec http://sourceforge.net/projects/platec/ .
That project is part of a Bachelor of Engineering thesis in Metropolia University of Applied Sciences, Helsinki, Finland. The thesis is freely downloadable from http://urn.fi/URN:NBN:fi:amk-201204023993 .

Kudos to the original author: Lauri Viitanen!

