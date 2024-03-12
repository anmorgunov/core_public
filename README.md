# Calculations of Core Excitations

[![Checked with mypy](https://www.mypy-lang.org/static/mypy_badge.svg)](https://mypy-lang.org/)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

This is part of the code extracted from a private repository.

## Structure

> Please note that each module has its own README file. (TODO: update outdated readmes for calculations/analysis)

- [basis_sets](basis_sets): module to create custom basis sets
- [geom_optimization](geom_optimization): module to optimize geometries
- [calculations](calculations): module to actually run calculations
- [analysis/parser](analysis/parser): module that allows you to parse all experimental results and store in an Excel file
- [analysis/regression](analysis/regression): module that builds linear regressions and compares values to experiment

## Getting Started

Clone repository into a folder `core_excitations`:

```bash
git clone git@github.com:anmorgunov/cebe_prediction.git
```

Rename `PRIVATE(example).py` to `PRIVATE.py`. Replace contents in `BASE_PATH` to specify path where your `core_excitations` folder is located. I use absolute paths for clarity and to avoid saving files exported by one of the modules in another.

In this project, the path is always specified as a list of folder names and/or file name. For example, `/Users/user/path/file.xlsx` is represented by `["", "Users", "user", "path", "file.xlsx"]` This is my way of making scripts cross-platform: UNIX-based systems (macOS, Linux) use `/` as separator in path, while Windows uses `\`. `os.sep` can give you the proper separator and I use it in `tools.join_path` to create a path with the proper separator.

## General Procedure

1. Get your geometries (see `geom_optimization`). Collect them as `.xyz` files in a folder.
2. Place geometries in `calculations/mols/`.
3. Create new (or modify existing) `calculations/mom_atom.sh`
4. Copy files to server and `bash mom_atom.sh`
5. Bring back calculation folders from server and place them in `calculations/mom/`
