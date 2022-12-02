# Visualizer module

## General comments

As in parser module, molecule names are determined by the folder names in which the results of calculations are stored

## Procedure

1. Add molecule list to `my_constants.py`
2. Run `analyze.py`
3. Run `evaluate.py`

## Analyze.py

Open `my_constants.py` and in dictionary `ATOM_TO_MOLS` add a new key `{atom}` which maps to a string of molecules that you have Core IP calculations for. The script (analyze.py) will search for those results in `calculations/{algorithm}/{atom}.lower()/`.

Analyze.py runs over all key:val pairs in ATOM_TO_MOLS, and for each creates an object with all CEBE values. It then can build linear regressions based on those values. Function `parser()` specifies key parameters:

```python
methods = 'UHF MP2 MP3'.split()
method2 = 'CCSD'
bases = 'D T'
algorithm = 'mom'
```

`method2` is what we compare against. `methods` is what we compare. I.e. with this setup, you'll get graphs for `UHF(x) vs CCSD (y)`, `MP2(x) vs CCSD (y)`, and `MP3(x) vs CCSD (y)`. `algorithm` is where the results in `calculations/{algorithm}` are stored. `bases` is equivalent to `zetas` from `calculations` launch scripts.

The graphs are stored in three formats: `html`, `svg`, and `jpg`. `html` plots are interactive, so you can point over any given data point and see what molecule is that.

## Evaluate.py

Similarly to `analyze.py`, `evaluate.py` runs over all key:val pairs in `my_constants.ATOM_TO_MOLS` and parses all Core IP results, gets experimental data from `analysis/parser/experimental.xlsx` and does the comparison. The summary is printed as an excel file (via `self.export_deviations_to_excel()` in `main()`) that is stored in `exports/evaluation/summary.xlsx`.

You also get a few statistical metrics in `exports/deviations_summary.md` (via `self.summary_of_deviations()`)

Currently, evaluate.py has also code for using one regression equation against another dataset, but I'll write documentation for it on some other day.

## To-do

1. Alter `evaluate.py` so that it can handle many methods (differentiate between them)
