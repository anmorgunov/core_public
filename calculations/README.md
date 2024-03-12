# Calculations Module

## General Overview

To start CEBE calculations you need to run `run_mom.py` specifying a series of command-line arguments. Here's a list of **required** arguments:

- `geom=/path/to/geom.xyz` - path to the geometry (an XYZ file) of the molecule of interest
- `atom=O` - you need to specify an element for which CEBE is calculated
- `orbital=0` - you need to specify the index (0-based) of the molecular orbital from which you want to extract an electron
- `corebasis=cc-pCVTZ` - basis set that will be placed on element for which CEBE is calculated
- `regularbasis=cc-pVTZ` - basis set that will be placed on all other elements

In addition, you can specify **optional** parameters:

- `localize=True` (defaults to *False*) - whether Boys localization will be run on the ground-state wavefunction. This is required for accurate prediction of CEBE for molecules that have symmetry beyond C1.
- `scf_maxcycles=100` (defaults to *1000*) - the upper limit on the number of SCF iterations. If the SCF equations do not converge after `scf_maxcycles`, the calculation will be terminated.
- `diis_start=5` (defaults to *0*) - iteration of UHF on which DIIS will be enabled
- `diis_space=20` (defaults to *50*) - size of the DIIS space during UHF
- `level_shift=0.005` (in a.u., defaults to *0*) - a level shift to a core orbital could be applied (sometimes helps with convergence)

An important **optional** parameter is the specification of the general config:

- `config=TEST` - which config, defined in [configs/setups.py](configs/setups.py), will be used. Available keys: `TEST`, `MP3`, `LOCAL`. If nothing is specified, the general config will be used.

## Preparing Calculation Scripts

It's pretty straightforward to prepare calculation scripts for any large number of molecules. Let's say you want to calculate CEBE for O-based core orbitals. Place geometries (as `.xyz` files) of all molecules in `mols/exp/o/`, and on the second line (the comments line) you should write `X O Y`, where `X` used to be a way of specifying total number of electrons[^2], `O` in this case specifies element that will be treated with `corebasis`, and `Y` specifies the index (0-based) of the orbital that should be ionized. After that, create a new block in [launch_parameters.txt](launch_parameters.txt):

```bash
START
atom: o
bases: STO-3G STO-6G 
mols: h2o c2h5oh
mols(loc): co2
time: 999
memory: 200
cores: 24
partition: high
END
```

- The `atom` determines the `mols/exp/${atom}` path and also the `mom/${atom}` where calculation results will be written.
- `bases` specifies basis sets (separated by space) in which CEBEs will be calculated. Available options: `STO-3G STO-6G 3-21G 4-31G 6-31G def2svp def2svpd` - if one of these is selected, both `regularbasis` and `corebasis` will have the same value. If, however, any of `D T Q 5` are specified, the `cc-pVXZ/cc-pCVXZ` (with `X=D,T,Q,5`) will be used.
- `mols` - a list of molecules for which CEBE will be calculated. The name should be the same as that of the geometry file `mols/exp/${atom}/${mol}.xyz`.
- `mols(loc)` - similar to `mols`, but for all of these molecules `localize=True` will be set.
- `time` - time limit for the SLURM job
- `memory` (in GB) - memory that will be allocated for the job
- `cores` - number of cores that will be allocated for the job
- `partition` - with some clusters, you may need to specify the partition in your SLURM submission, so this is how you do it.

To create submission scripts, simply run `generate_bash.py`.

When you're with a given set of molecules, you can add a hashtag before `#START` and the whole block will be ignored on later calls to `generate_bash.py`.

## More on Configs

In general, the CEBE calculations will proceed with the following config (the `GEN_CONFIG` in `configs/setups.py`):

```python
GEN_CONFIG = dict(
    doMP2 = True,
    doMP3 = True,
    doCCSD = True,
    doTriples = True,
    doSFX = True,
    toFreeze = True,
    doChengBasis = False,
    doSpecialBasis = False,
    toPrintDensity = False,
)
```

where doMP2, doMP3, doCCSD are self-explanatory, doTriples controls whether perturbative triples correction to CCSD is calculated, doSFX controls whether SFX2C-1e correction is applied, toFreeze controls whether ionized orbital should be frozen for CCSD calculations (helps with convergence)[^1]. If toPrintDensity is set to True, a cube file with difference in density between ionized and ground states will be written.

If calculation in quadruple-zeta basis (`cc-pVQZ`) is requested, the following config is applied:

```python
Q_CONFIG = dict(
    doMP3 = False,
)
```

For every parameter not specified in `Q_CONFIG`, or any other custom config, the value will be taken from `GEN_CONFIG`. In other words, it's enough to specify what's different. MP3 calculations are disabled (by default) for quadruple bases because MP3 code is not optimized for memory usage and may easily use >300 GB of disk space even for a simple molecule.

Similarly, only MP2 calculations are executed in pentuple basis `cc-pV5Z`:

```python
PENTUPLE_CONFIG = dict(
    doMP3 = False,
    doCCSD = False,
    doTriples = False
)
```

For easy debugging of the code on a local machine, the following config can be called:

```python
LOCAL_CONFIG = dict(
    doMP2=True,
    doMP3=False,
    doCCSD=False,
    doTriples=False,
    doSFX=True,
    toFreeze=True,
)
```

Of course you can modify any of the configs above or add new configs in [configs/setups.py](configs/setups.py).

### Custom Basis Sets

You can run calculation with any other custom basis set if you store it in a `nwchem` folder and adjust launch files accordingly. (see also [basis_sets](../basis_sets/) module).

[^1]: doChengBasis and basisLibKey were implemented to suppport custom basis sets (at some point, I tried to use the Dunning basis sets recontracted for SFX2C1e scheme, used by Lan Cheng in his 2019 paper). It may be worth rewriting the custom basis support in better way, so I'll keep these parameters as *legacy* for now.
[^2]:  this is really legacy and isn't used anymore, I will remove it at some point. All existing geometries, however, do specify a correct value of X.
