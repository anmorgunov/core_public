import os
import sys
sys.path.append('../core_excitations') # this allows us to import files from parent directories
import PRIVATE
import tools 

ATOM='newf'

SUBMODULE = "geom_optimization" #the name of the current folder
BASE = PRIVATE.BASE_PATH + [SUBMODULE,]

def create_orca_input(path, geometry, molecule):
    keywords = '''! RI-MP2 cc-pVQZ cc-PVQZ/C Opt
%maxcore 4000
%pal nprocs 16 end

'''
    geom = f'''* xyz 0 1
{geometry}*
'''
    with open(tools.join_path(path+[f'{molecule}.inp',]), 'w') as f:
        f.write(keywords)
        f.write(geom)

def create_slurm_script(path, molecule):
    script = f'''#!/bin/sh

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
'''
    with open(tools.join_path(path+['submit.sh',]), 'w') as f:
        f.write(script)

if f"run-{ATOM}" not in os.listdir(tools.join_path(BASE + ["orca",], endWithSlash=True)): os.mkdir(tools.join_path(BASE + ["orca", f"run-{ATOM}"], endWithSlash=True))

MOLS = os.listdir(tools.join_path(BASE+["inits", f"run-{ATOM}"], endWithSlash=True))
query = " ".join([mol.split('.')[0] for mol in MOLS])
print(query)
for molExt in MOLS:
    path = tools.join_path(BASE + ["orca", f"run-{ATOM}"], endWithSlash=True)
    mol = molExt.split('.')[0]
    if mol not in os.listdir(path): os.mkdir(tools.join_path(BASE + ["orca", f"run-{ATOM}", mol], endWithSlash=True))
    with open(tools.join_path(BASE + ["inits", f"run-{ATOM}", molExt]), 'r') as f:
        lines = f.read().split('\n')
        geometry = '\n'.join(lines[2:])

    create_orca_input(BASE + ["orca", f"run-{ATOM}", mol], geometry, mol)
    create_slurm_script(BASE + ["orca", f"run-{ATOM}", mol], mol)

# cf3cooh ch3cooch3
# pnh2pyridi-n-e me2-n-cho pohpyridine 2pyrido-n-e onh2pyrdi-n-e mnh2pyridi-n-e 2pyrido-n-earo pfpyridine ofpyridine ch3sc-n
# cf3ocf3 cf3cch cf3ocf3tw cf3chch2