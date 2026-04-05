Plate Tectonics
===============

This fork is focused on the C++ library, the simulation tool, and the test suite. Python bindings are removed.

Layout
======

- `src/`: core library
- `tools/`: `simulation.exe` and shared heightmap I/O code
- `tests/`: GoogleTest suite
- `img/in/`: timestamped initial-state exports and source input images
- `img/out/`: timestamped final-state exports
- `img/gif/`: generated GIFs
- `img/frames/`: kept intermediate frames
- `x64-Debug/`, `x64-Release/`: the only intended local build folders

Build
=====

Use x64 Native Tools PowerShell for Visual Studio 2026:

```powershell
cmake --preset x64-debug
cmake --build --preset build-x64-debug
ctest --preset test-x64-debug
```

`x64-debug` intentionally uses `RelWithDebInfo`, not plain `Debug`, so runtime performance stays close to Release while PDB symbols remain available.

The Windows presets expect vcpkg at `C:/dev/vcpkg` and use `vcpkg.json` for `libpng`.

Manual commands are collected in `build_commands.txt`.

Executables
===========

- `simulation.exe`: the interactive/manual runner. It generates terrain, imports heightmaps, writes outputs, optionally saves frames, and optionally builds a GIF.
- `PlateTectonicsTests.exe`: the automated test binary used by CTest. It runs assertions against library behavior and writes regression artifacts for verification. It is not the simulator UI/tool.

Output Rules
============

Each run gets one id: `<YYMMDDHHMM>_<SEED-or-INPUTIMGNAME>`.

- Initial state: `img/in/<id>.r16` and `img/in/<id>.png`
- Final state: `img/out/<id>.r16` and `img/out/<id>.png`
- GIF: `img/gif/<id>.gif`
- Frames: `img/frames/<id>_<NUMOFFRAME>.png`
- Raw metric exports also write sidecars next to `.r16`: `.r16.json`
- `--export-heightmap-f32` writes `img/out/<id>.f32` and `.f32.json`

The `.r16` file is the canonical metric heightmap export. The `.png` file is the rendered preview image for the chosen display settings.

Simulation Usage
================

Example:

```powershell
.\x64-Debug\bin\simulation.exe --dim 600 400 -s 12345 --gif --step 25
```

Arguments
=========

- `-h`, `--help`: show help
- `-s SEED`: set the random seed
- `-i FILE`, `--input FILE`: load a legacy normalized PNG from the given path or from `img/in/`
- `--input-png16 FILE`: load a 16-bit grayscale metric PNG
- `--input-r16 FILE`: load a little-endian metric `.r16` heightmap
- `--dim WIDTH HEIGHT`: set dimensions when no input image is used
- `--sea-level-m N`: set metric coastline threshold in `[0, 65535]`
- `--colors`: render preview PNGs/frames in color
- `--grayscale`: render preview PNGs/frames in grayscale
- `--no-output-normalization`: disable first-frame normalization for preview rendering
- `--min-initial-height X`: remap the initial minimum height to `X` during preview normalization
- `--max-initial-height X`: remap the initial maximum height to `X` during preview normalization
- `--display-min X`: lower bound of the preview display window
- `--display-max X`: upper bound of the preview display window
- `--tone-map MODE`: preview tone mapping mode: `linear`, `log`, or `asinh`
- `--export-heightmap-f32`: export the final raw float32 heightmap and metadata
- `--filename NAME`: deprecated and ignored for on-disk naming
- `--cycles N`: number of simulation cycles
- `--plates N`: number of tectonic plates
- `--aggregation-overlap-abs N`: continent aggregation overlap threshold in pixels
- `--aggregation-overlap-rel X`: continent aggregation overlap threshold as a ratio
- `--folding-ratio X`: fraction of overlapping continental crust converted into uplift
- `--erosion-period N`: updates between erosion passes
- `--erosion-strength X`: erosion strength multiplier; `0` disables erosion
- `--rotation-strength X`: angular plate motion multiplier
- `--landmass-rotation X`: visible crust rotation multiplier; `0` disables it
- `--subduction-strength X`: oceanic crust removal multiplier during subduction
- `--step X`: save one intermediate frame every `X` simulation steps to `img/frames/`
- `--gif`: create a GIF in `img/gif/`; with `--step`, the GIF uses sampled frames, otherwise it uses every update
- `--no-steps`: after GIF creation, delete the generated frame PNGs instead of keeping them
- `--show-boundaries`: overlay convergent, divergent, and transform boundaries on preview PNGs, frames, and GIF frames

Testing
=======

```powershell
ctest --preset test-x64-debug
```

Code Quality
============

`.clang-format` and `.clang-tidy` are kept. The old shell wrappers were removed because they were stale, Unix-only, and tied to the previous layout.
