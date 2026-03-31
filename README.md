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
cmake ..
cmake --build .
```

On Windows, the repository uses `vcpkg.json` manifest mode to fetch `libpng`
automatically. A fresh configure like this is the intended path:

```powershell
cmake -S . -B build -DWITH_EXAMPLES=ON -DMSVC_RUNTIME=dynamic
cmake --build build --config Release
```

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

From the root directory run:

```bash
mkdir -p build
cd build
cmake .. -DWITH_EXAMPLES=ON -G "Unix Makefiles"
make
cd examples
./simulation
```

The example writes `.png` files.

The default CLI is deterministic. If you do not pass `-s` / `--seed`, it uses seed `1`.

You can also seed the simulation from an image heightmap:

```bash
./simulation --input-image input.png --filename transformed
./simulation --input-image input.png --grayscale --show-fault-lines --filename transformed
```

Input image rules:

* PNG input is supported
* RGB images are converted to grayscale using luminance
* Black = low terrain, white = high terrain
* `sea_level` still determines how much of the image becomes ocean vs land
* In image mode, the defaults are tuned down to `4` plates, `1` cycle, and `0` erosion

Useful knobs:

* `--num-plates N`
* `--cycle-count N`
* `--erosion-period N`
* `--sea-level F`
* `--folding-ratio F`
* `--aggr-overlap-abs N`
* `--aggr-overlap-rel F`
* `--progress-interval N`
* `--step N`
* `--show-fault-lines` writes additional `*_faults.png` files

`--cycle-count` is the number of tectonic lifetimes the simulation is allowed to
run. When plate motion stalls or the system cools down too much, the library can
restart with a fresh split of plates and continue from the updated world.

Fault line colors in `*_faults.png`:

* convergent = red
* divergent = cyan
* transform = magenta

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

