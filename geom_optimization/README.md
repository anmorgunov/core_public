# Geometry Optimization Module

> Scripts to automate geometry optimization

For core excitations, I try to use experimental geometries whenever they're available. I usually check the [CCCBDB](https://cccbdb.nist.gov/expgeom1x.asp). Even if they don't have an experimental geometry, they have a huge DB of calculated geometries, and I usully pick MP2(Full)/cc-pVT/QZ.

In those cases when CCCBDB has neither experimental nor calculated geometries, I optimize those structures on my own with ORCA.

## Procedure

1. Save initial guesses as regular `.xyz` files and place them in a folder within `inits`. Your folder should be named `run-{suffix}` with any desired `{suffix}` (e.g. `reopt`)
2. Open `create_orca_inp.py` and modify `ATOM` variable at the very top to be equal to the chosen `{suffix}`
3. Run `create_orca_inp.py`. The submit files for geometries will be in `orca/run-{suffix}/` folder.
4. Copy them to cluster, and launch. You may find `run_geom.sh` helpful.
5. ORCA will save the latest (optimized) geometry in `geom.xyz` file.
6. In `calculations/mols/` create a folder in which you're going to place your geometries.

The general structure for your `.xyz` file in `calculations/mols`:

```text
numatoms
X A Y
geometry...
```

where `X` is the number of electrons[^1], `A` is the symbol of an element to be excited (this decides which element gets the `cc-pCVXZ` basis set, while others get `cc-pVXZ`), and `Y` is the index of orbital to be excited (counting from 0).

[^1]: okay, to be honest, you don't really need it, since PySCF can get that number. I don't know how I started to put them there, but I just did. It just happened that I had already written a python script that could calculate the number of electrons based on chemical formula.

## Naming convention

Proper filename for an `.xyz` file in `calculations/mols/` is essential.

1. Filenames should be unique.
2. In `analysis/parser/experimental.xlsx` you'll find experimental values for Core IP and filenames in an adjacent column. The scripts in `analysis` use that file to get experimental values, and they search for proper molecule by filename.
3. Because you specify which atom gets the `C` basis set (cc-p**C**VXZ) in an `.xyz` file, you obviously need separate files for excitations of different atoms. Here's how I handle that:
    - For the excitation of the heaviest atom (say Cl in CCl4) just write the molecular formula `ccl4.xyz`
    - For all other atoms, separate the atom of excitation by a dash. Say, if you want to excite C in CCl4: `c-cl4.xyz`. Excite O in N2O? `nno.xyz`. Excite middle nitrogen in N2O? `n-n-o.xyz`. Excite leftmost nitrogen? `n-no.xyz`.

## Logbook

I keep track of which geometries I used in `geometriesDB.xlsx`.

## ORCA vs Q-Chem

There are scripts to create input files for both Q-Chem and ORCA. I found Q-Chem to be surprisingly slower than ORCA. Here's the running time (in minutes) for optimization at RI-MP2 level (cc-pVQZ cc-PVQZ/C basis sets). To surprise you even more, Q-Chem calculations were allotted 16 cores and ORCA had only 8.

While RIJCOSX (with def2/J) speeds up calculations dramatically, per Prof. Van Voorhis suggestion, I used the results from calculations without it.

|           | QChem | ORCA | ORCA (RIJCOSX) |
|-----------|-------|------|----------------|
| c5h5cl    | 612   | 107  | 17             |
| c5h5f     | 524   | 69   | 15             |
| c6h11cl   | 520   | 316  | 45             |
| ch3cosch3 |       | 101  | 14             |
| ch3cosh   | 171   | 46   | 11             |
| ch3cssch3 | 710   | 159  | 20             |
| ch3cssh   | 260   | 41   | 11             |
| ch3nhcl   | 74    | 14   | 5              |
| ch3nhcocl | 274   | 70   | 14             |
| ch3nhcof  | 231   | 48   | 7              |
| ch3nhf    | 62    | 9    | 4              |
| ch3of     | 34    | 9    | 3              |
| hcosch3   | 256   | 48   | 8              |
| hcosh     | 70    | 11   | 5              |
| hcssch3   | 192   | 68   | 8              |
| hcssh     | 57    | 23   | 5              |
| hoch2cl   | 69    | 16   | 5              |
| nh2ch2cl  | 90    | 25   | 6              |
| nh2ch2f   | 46    | 12   | 3              |
