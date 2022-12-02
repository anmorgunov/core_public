# Parser Module

## Some general ideas

- Molecule names are determined during calculations step (i.e., they're determined by the folder name in which the results are stored)
- The parser takes everything contained in folders specified by `METHOD_TO_FOLDERS` dictionary in `perform_parsing()` function. If `METHOD_TO_FOLDERS = {'mom': ['cl', 'o']}`, then you should have calculation results in `calculations/mom/cl/` and `calculations/mom/o/`

## Procedure

1. Ensure you have a CEBE_Data.xlsx file (even if its empty)[^1]. You could simply rename `CEBE_Data(sample).xlsx` to `CEBE_Data.xlsx`.
2. Run `parse.py` to update CEBE_Data.xlsx file.
3. Use `extract.py` to store results of any desired molecules in a separate file[^2].

## Current dependencies

- `parse.py` reads an atom to be excited from the last line in `submit.sh`

## Legacy files

- `legacy/parse_single.py` is a small script I wrote when I needed to put a few results in a table. It's from back in the day when `parse.py` wasn't as universal as it is now.
- `legacy/parse_ccsd.py` is a script I wrote to take FC-CVS-EOM results from QChem. Didn't use it as much, but maybe it'll be useful one day.

[^1]: I've implemented it in the way that you can leave notes next to a result in column `S` (as of aug 2022) and those notes will be saved whenever you parse new results.
[^2]: My idea was to make a `CEBE_data.xlsx` a giant file with all results stored in one place. If you want to do something with those results, you could just parse it. If you want to focus on a subset of molecules, you just extract them into a separate `excel` file.
