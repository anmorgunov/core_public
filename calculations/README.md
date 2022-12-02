# Calculations Module

## Launch Script

I usually keep separate launch scripts for each atom-series. I don't mind having 10 files and if I only want to run a few more molecules, I don't have to worry about ensuring I modified all the paths and parameters properly.

### General Structure

```bash
run_path=$PWD
geom_path="$run_path/mols/exp/f" #where you store geometries
out_path="$run_path/mom/f/" # where output files will be created

run_time="100:00:00" # time limit
nproc=8 # num of CPUs (I usually use 8)
nodes=1
mem0="8000" # memory per core (!)
mem=`echo $nproc $mem0 | awk '{print $1 * $2}'`

zetas="D T"
methods="mom"
pyscriptsuffix="frozen"

mols="c2h3f c2h5f c3h7f c5h5f c6h5f c6h11f c6h11of ch3cof ch3nhcof ch3nhf ch3of foh hccf hoch2f nh2ch2f nh2f ch3f"
```

At the end, the launch script will create a `submit.sh` file that'll call `run_{method}_{pyscriptsuffix}.py` file. That's the purpose of `methods` and `pyscriptsuffix` parameters.

A few words on `zetas`. By default, if you don't have `doChengBasis=True` or `doSpecialBasis=True` (see below), the launch script will instruct to run calculation with `cc-pCV{zeta}Z` on atoms of the type you excite (i.e. on all O atoms if you excite O in CO2) and `cc-pV{zeta}Z` on other atoms. So by iterating over `zetas` you get the `zeta`.

Outputs will be stored in `$out_path/{moleculename}/{basis}/`.

## Calculations Launch Files

Currently, you have two of them: `run_mom.py`. Both have the same set of parameters:

| Parameter | Possible Values | Notes |
| - | - | - |
| `doMP` | `True/False` | if set to True, both MP2 and MP3 will be performed |
| `doCCSD` | `True/False` | if set to True, CCSD and CCSD(T) will be performed. `doMP` must be True if CCSD is to be used |
| `doSFX` | `True/False` | if set to True, SFX2C-1e correction will be applied[^1] |
| `toFreeze` | `True/False` | if set to True, the orbital to be excited will be frozen |
| `doChengBasis` | `True/False` | if set to True, the SFX2C-1e recontracted Dunnig Triple zeta basis will be used |
| `doSpecialBasis` | `True/False` | if set to True, the special basis from `externalBasisLibrary` will be used. Must be accompanied by `basisLibKey` |
| `basisLibKey` | `str` | a key for the special basis |
| `toPrintDensity` | `True/False` | if set to True, a cube file with difference in density between ionized and ground states will be written |

Recommended parameters:

```python
doMP = True
doCCSD = True
doSFX = True 
toFreeze = True
doChengBasis = False
doDeconT = False
toPrintDensity = False # cube files are rather large, so I don't print them unless necessary
```

### Custom Basis Sets

You can run calculation with any other custom basis set if you store it in a `nwchem` folder and adjust launch files accordingly. (see also [basis_sets](../basis_sets/) module).

## Current dependencies

If you do MOM calculations, it's best (actually, you have to, otherwise you'll need to tweak scripts) to store all results in a `mom` folder which has `{atom}` subfolders.

[^1]: I haven't seen any (literally ANY) effect introduced by SFX2C-1e. This, of course, could mean I'm doing something wrong.

## To-do

Make sure localization is generalized
Check population analysis for C-localization
