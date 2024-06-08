import os
import sys

sys.path.append(
    "../core_excitations"
)  # this allows us to import files from parent directories
import tools

import PRIVATE

ATOM = "test"

SUBMODULE = "geom_optimization"  # the name of the current folder
BASE = PRIVATE.BASE_PATH + [
    SUBMODULE,
]


def create_qchem_input(path, geometry):
    rem = """$rem
   JOBTYPE Opt
   METHOD riMP2
   BASIS cc-pVQZ
   AUX_BASIS rimp2-cc-pVQZ
$end\n\n"""
    geom = f"""$molecule
0 1
{geometry}$end
"""
    with open(tools.join_path(path + ["qchem.in"]), "w") as f:
        f.write(rem)
        f.write(geom)


def create_slurm_script(path):
    script = """#!/bin/sh

#SBATCH -t 100:00:00
#SBATCH -o sbatch.out
#SBATCH -e sbatch.err
#SBATCH -c 16
#SBATCH -N 1
#SBATCH --mem=64000

export OMP_NUM_THREADS=8
export MKL_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=8
export NUMEXPR_NUM_THREADS=8

echo `date` >> "time.txt"
qchem.latest qchem.in qchem.out
echo `date` >> "time.txt"
"""
    with open(tools.join_path(path + ["submit.sh"]), "w") as f:
        f.write(script)


if f"run-{ATOM}" not in os.listdir(
    tools.join_path(
        BASE
        + [
            "qchem",
        ],
        endWithSlash=True,
    )
):
    os.mkdir(tools.join_path(BASE + ["qchem", f"run-{ATOM}"], endWithSlash=True))

MOLS = os.listdir(tools.join_path(BASE + ["inits", f"run-{ATOM}"], endWithSlash=True))
query = " ".join([mol.split(".")[0] for mol in MOLS])
print(query)
for molExt in MOLS:
    path = tools.join_path(BASE + ["qchem", f"run-{ATOM}"], endWithSlash=True)
    mol = molExt.split(".")[0]
    if mol not in os.listdir(path):
        os.mkdir(
            tools.join_path(BASE + ["qchem", f"run-{ATOM}", mol], endWithSlash=True)
        )
    with open(tools.join_path(BASE + ["inits", f"run-{ATOM}", molExt]), "r") as f:
        lines = f.read().split("\n")
        geometry = "\n".join(lines[2:])

    create_qchem_input(BASE + ["qchem", f"run-{ATOM}", mol], geometry)
    create_slurm_script(BASE + ["qchem", f"run-{ATOM}", mol])
