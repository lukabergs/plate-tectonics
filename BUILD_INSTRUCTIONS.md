# Build Instructions

This repository must use the committed presets from `CMakePresets.json`.

Rules
=====

- Supported build folders are `out/build/x64-Debug/` and `out/build/x64-Release/`.
- Do not configure into repo-root `x64-Debug/`, repo-root `x64-Release/`, `build/...`, or any other ad-hoc folder.
- `x64-debug` intentionally maps to `RelWithDebInfo`, not plain `Debug`.
- The Windows toolchain expects vcpkg at `C:/dev/vcpkg`.

Command Line
============

Run these from x64 Native Tools PowerShell or x64 Native Tools Command Prompt for Visual Studio 2026:

```powershell
cmake --preset x64-debug
cmake --build --preset build-x64-debug
ctest --preset test-x64-debug
.\out\build\x64-Debug\bin\simulation.exe --dim 600 400 -s 12345
```

Release:

```powershell
cmake --preset x64-release
cmake --build --preset build-x64-release
ctest --preset test-x64-release
```

Visual Studio 2026
==================

- Open the repo as a CMake project or open the solution workspace.
- Select the `x64-Debug` or `x64-Release` preset.
- Build, test, and run from that preset only.
- If Visual Studio starts targeting any other folder, switch back to the preset-based configuration and reconfigure.

VS Code
=======

- The workspace setting in `.vscode/settings.json` forces CMake Tools to use presets.
- In the VS Code terminal, use the same `cmake --preset ...` and `cmake --build --preset ...` commands shown above.

Cleaning
========

If stale folders were already created, remove them explicitly:

```powershell
Remove-Item -Recurse -Force .vs, out, build, x64-Debug, x64-Release -ErrorAction SilentlyContinue
```

One Source Of Truth
===================

If any LLM or manual edit changes the build location, restore it so that:

- `x64-debug` configures into `C:\dev\plate-tectonics\out\build\x64-Debug`
- `x64-release` configures into `C:\dev\plate-tectonics\out\build\x64-Release`
