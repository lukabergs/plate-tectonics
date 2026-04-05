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
- `out/build/x64-Debug/`, `out/build/x64-Release/`: the only intended local build folders

Build
=====

See `BUILD_INSTRUCTIONS.md`. That file is the only supported build/test reference.

Executables
===========

- `simulation.exe`: the interactive/manual runner. It generates terrain, imports heightmaps, writes outputs, optionally saves frames, and optionally builds a GIF.
- `PlateTectonicsTests.exe`: the automated test binary used by CTest. It runs assertions against library behavior and writes regression artifacts for verification. It is not the simulator UI/tool.

Output Rules
============

Each run gets one id: `<YYMMDDHHMMSS>_<SEED-or-INPUTIMGNAME>`.

- Initial state: `img/in/<id>.r32` and `img/in/<id>.tiff`
- Final state: `img/out/<id>.r32` and `img/out/<id>.tiff`
- GIF: `img/gif/<id>.gif`
- Frames: `img/frames/<id>_<NUMOFFRAME>.png`
- Raw metric exports also write sidecars next to both `.r32` and `.tiff`: `.r32.json` and `.tiff.json`
- `--export-heightmap-f32` writes `img/out/<id>.f32` and `.f32.json`

The `.r32` file and the `.tiff` file are the authoritative 32-bit heightmap exports. The `.tiff` file stores grayscale float32 samples directly, not a preview render. Human preview images are only written to `img/frames/` and into GIFs.

Simulation Usage
================

Example:

```powershell
.\out\build\x64-Debug\bin\simulation.exe --dim 600 400 -s 12345 --gif --step 25
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
- `--min-initial-height X`: minimum metric height used when procedurally seeding the initial world; default `10000`
- `--max-initial-height X`: maximum metric height used when procedurally seeding the initial world; default `30000`
- `--display-min X`: lower bound of the preview display window
- `--display-max X`: upper bound of the preview display window
- `--tone-map MODE`: preview tone mapping mode: `linear`, `log`, or `asinh`
- `--export-heightmap-f32`: export the final raw float32 heightmap and metadata
- `--filename NAME`: deprecated and ignored for on-disk naming
- `--cycles N`: number of simulation cycles; default `2`
- `--plates N`: number of tectonic plates; default `10`
- `--aggregation-overlap-abs N`: absolute collision coverage needed to aggregate continents; default `max(64, WIDTH*HEIGHT/1000)`, which is `240` at `600x400`
- `--aggregation-overlap-rel X`: relative collision coverage needed to aggregate continents in `[0, 1]`; default `0.20`
- `--folding-ratio X`: fraction of overlapping continental crust converted into uplift in `[0, 1]`; default `0.08`
- `--erosion-period N`: simulation updates between erosion passes; default `60`
- `--erosion-strength X`: erosion strength multiplier; `0` disables erosion; default `1.0`
- `--rotation-strength X`: angular plate motion multiplier; default `1.0`
- `--landmass-rotation X`: visible crust rotation multiplier; `0` disables it; default `0.20`
- `--subduction-strength X`: oceanic crust removal strength in `[0, 1]`; `1.0` removes one full `OCEANIC_BASE` chunk per subduction event; default `1.0`
- `--step X`: save one intermediate frame every `X` simulation steps to `img/frames/`
- `--gif`: create a GIF in `img/gif/`; with `--step`, the GIF uses sampled frames, otherwise it uses every update
- `--no-steps`: after GIF creation, delete the generated frame PNGs instead of keeping them
- `--show-boundaries`: overlay convergent, divergent, and transform boundaries on preview PNGs, frames, and GIF frames

Default Tuning
==============

- `2` cycles and `10` plates keep the default run active enough to form varied interactions without making the out-of-the-box run slow.
- `0.08` folding, `0.20` relative aggregation, and the area-scaled absolute aggregation threshold make continents fuse after sustained contact instead of effectively disabling aggregation on larger maps.
- `60`/`1.0` erosion, `1.0` angular rotation, `0.20` landmass rotation, and `1.0` subduction keep all major processes visible by default without one effect overwhelming the others.

Code Quality
============

`.clang-format` and `.clang-tidy` are kept. The old shell wrappers were removed because they were stale, Unix-only, and tied to the previous layout.
