# Basis Sets Module

> Scripts written to convert basis set description from CFOUR to NWChem (suitable for PySCF).

For their calculations, [Zheng and Cheng](https://pubs.acs.org/doi/abs/10.1021/acs.jctc.9b00568) used correlation-consistent basis sets recontracted for the SFX2C-1e scheme. I found those basis sets on the CFOUR.de website, but they were in "GENBAS" format while PySCF takes basis sets in NWChem format. These scripts were written to accomplish that.

The way I did this:

1. Place CFOUR files in `cfour` folder
2. Run `convert.py`
3. The converted files are in `nwchem` folder, which, because of `__init__.py` serves as a submodule. You can `scp` the whole folder to the directory from which you run core calculations and do `import nwchem` in your launch `.py` files.
