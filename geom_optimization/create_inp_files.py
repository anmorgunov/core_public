from pathlib import Path

base_path = Path(__file__).resolve().parent
folder_name = "test"

qchem_path = base_path / "qchem" / folder_name
qchem_path.mkdir(exist_ok=True)


def create_qchem_input(path: Path, geometry: str) -> None:
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
    with open(path / "qchem.in", "w") as f:
        f.write(rem)
        f.write(geom)


def create_slurm_script(path: Path) -> None:
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
    with open(path / "submit.sh", "w") as f:
        f.write(script)


molecules = [f.name for f in (base_path / "inits" / folder_name).glob("*.xyz")]
query = " ".join([mol.split(".")[0] for mol in molecules])
print(query)
for mol_file in molecules:
    mol = mol_file.split(".")[0]
    molPath = qchem_path / mol
    molPath.mkdir(exist_ok=True)

    with open(base_path / "inits" / folder_name / mol_file, "r") as f:
        lines = f.read().split("\n")
        geometry = "\n".join(lines[2:])

    create_qchem_input(qchem_path / mol, geometry)
    create_slurm_script(qchem_path / mol)
