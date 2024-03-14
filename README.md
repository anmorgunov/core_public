# CEBE Prediction with MP2 Composite Schemes

[![Checked with mypy](https://www.mypy-lang.org/static/mypy_badge.svg)](https://mypy-lang.org/)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

## Abstract

X-ray photoelectron spectroscopy (XPS) measures core-electron binding energies (CEBEs) to reveal element-specific insights into chemical environment and bonding. Accurate theoretical CEBE prediction aids XPS interpretation but requires proper modeling of orbital relaxation and electron correlation upon core-ionization. This work systematically investigates basis set selection for extrapolation to the complete basis set (CBS) limit of CEBEs from ΔMP2 and ΔCC energies across 94 K-edges in diverse organic molecules. We demonstrate that an alternative composite scheme using ΔMP2 in a large basis corrected by ΔCC-ΔMP2 difference in a small basis can quantitatively recover optimally extrapolated ΔCC CEBEs within 0.02 eV. Unlike ΔCC, MP2 calculations do not suffer from convergence issues and are computationally cheaper, and, thus, the composite ΔMP2/ΔCC scheme balances accuracy and cost, overcoming limitations of solely using either method. We conclude by providing a comprehensive analysis of the choice of small and large basis sets for the composite schemes and provide practical recommendations for highly accurate (within 0.10-0.15 eV MAE) ab initio prediction of XPS spectra.

## Preprint

You can find the preprint for this work on [arXiv](https://arxiv.org/abs/2403.06364).

## Structure

> Please note that each module has its own README file

- [basis_sets](basis_sets): module to create custom basis sets
- [geom_optimization](geom_optimization): module to optimize geometries
- [calculations](calculations): module to actually run calculations
- [Analysis](Analysis): module that allows you to parse calculation results, extrapolate them to CBS, and reproduce all figures/tables from the paper
- [Data](Data): all relevant calculation results, XYZ files for geometries used and all figures/tables

## General Procedure for CEBE Calculation

TODO: Update the procedure... (the [Calculations README](calculations/README.md) is up to date though!)

1. Get your geometries (see `geom_optimization`). Collect them as `.xyz` files in a folder.
2. Place geometries in `calculations/mols/`.
3. Create new (or modify existing) `calculations/mom_atom.sh`
4. Copy files to server and `bash mom_atom.sh`
5. Bring back calculation folders from server and place them in `calculations/mom/`

## Tables/Figures Reproduction

In [perform_analysis.py](perform_analysis.py) you can find all the calls that correspond to:

1. Parsing of Calculation Results
2. Extrapolation of Calculation Results to the CBS
3. Creation of all Tables in the preprint
4. Creation of all Figures in the preprint
