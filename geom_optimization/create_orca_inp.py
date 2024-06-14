from pathlib import Path

base_path = Path(__file__).resolve().parent
folder_name = "test"

orca_path = base_path / "orca" / folder_name
orca_path.mkdir(exist_ok=True)


def create_orca_input(path: Path, geometry: str, molecule: str) -> None:
    keywords = """! RI-MP2 cc-pVQZ cc-PVQZ/C Opt
%maxcore 4000
%pal nprocs 16 end

"""
    geom = f"""* xyz 0 1
{geometry}*
"""
    with open(path / f"{molecule}.inp", "w") as f:
        f.write(keywords)
        f.write(geom)


def create_slurm_script(path: Path, molecule: str) -> None:
    script = f"""#!/bin/sh

#SBATCH -t 100:00:00
#SBATCH -o sbatch.out
#SBATCH -e sbatch.err
#SBATCH -J geom-{molecule}
#SBATCH -n 16
#SBATCH -N 1
#SBATCH --mem=64000

echo `date` >> "time.txt"
/home/morgunov/orca/orca {molecule}.inp > {molecule}.out
echo `date` >> "time.txt"
"""
    with open(path / "submit.sh", "w") as f:
        f.write(script)


molecules = [f.name for f in (base_path / "inits" / folder_name).glob("*.xyz")]
query = " ".join([mol.split(".")[0] for mol in molecules])
print(query)
for mol_file in molecules:
    mol = mol_file.split(".")[0]
    molPath = orca_path / mol
    molPath.mkdir(exist_ok=True)

    with open(base_path / "inits" / folder_name / mol_file, "r") as f:
        lines = f.read().split("\n")
        geometry = "\n".join(lines[2:])

    create_orca_input(orca_path / mol, geometry, mol)
    create_slurm_script(orca_path / mol, mol)

# cf3cooh ch3cooch3
# pnh2pyridi-n-e me2-n-cho pohpyridine 2pyrido-n-e onh2pyrdi-n-e mnh2pyridi-n-e 2pyrido-n-earo pfpyridine ofpyridine ch3sc-n
# cf3ocf3 cf3cch cf3ocf3tw cf3chch2
