import load_params
import os 
from pathlib import Path
import configs
import json
from typing import Dict

import re
def get_zeta(basis_set_name:str)->int:
    zeta_dict:Dict[str,int] = {
        'DZ': 2,
        'TZ': 3,
        'QZ': 4,
        '5Z': 5
    }
    match = re.search(r'pcX-(\d+)|cc-p?V?(DZ|TZ|QZ|5Z)', basis_set_name)
    if match:
        # Check which group got a match
        if match.group(1):
            # This is from 'pcX-n'
            return int(match.group(1))
        elif match.group(2):
            # This is from 'cc-pVXX'
            return zeta_dict[match.group(2)]
    raise ValueError(f"Could not determine zeta value for basis set {basis_set_name}")


class SubmitGenerator:
    """This is an object with which I create submission scripts
    """

    def __init__(self):
        pp = load_params.ParamParser()
        self.launchParams = pp.main()
        # self.BASE = PRIVATE.BASE_PATH + ["calculations"]
        self.atomToFname = {
            'c': 'carbon',
            'n': 'nitrogen',
            'o': 'oxygen',
            'f': 'fluorine',
            's': 'sulfur',
            'cl': 'chlorine',
        }
        self.fnames = set()
        self.clusterToFnames = {}
        self.cntr = 0

    def process_bases(self, bases):
        """In launch_parameters.txt you can use keyword "bases" and provide either a single letter {D, T, Q, 5} or a full name of a basis set.

        If you give X in {D, T, Q, 5}, this function interpretes it as: use cc-pCVXZ on core atom, cc-pVXZ on other atoms
        If you give any other basis set, it'll use that basis set on all atoms.

        Args:
            bases (str): names of basis sets separated by spaces

        Returns:
            tuple: a tuple of strings that denote bash arrays
        """
        zetaToBasis = {'D': ('cc-pVDZ', 'cc-pCVDZ'), # first element is a regular basis set, second element is something to place on a core atom
                       'T': ('cc-pVTZ', 'cc-pCVTZ'),
                       'Q': ('cc-pVQZ', 'cc-pCVQZ'),
                       '5': ('cc-pV5Z', 'cc-pCV5Z'),
        }
        basesList = bases.split()
        regList = [zetaToBasis[basis][0] if basis in zetaToBasis else basis for basis in basesList]
        coreList = [zetaToBasis[basis][1] if basis in zetaToBasis else basis for basis in basesList]
        markList = [basis if basis in zetaToBasis else 'None' for basis in basesList ]
        regStr = '(' + str(json.dumps(regList))[1:-1].replace(',', '') + ')'
        coreStr = '(' + str(json.dumps(coreList))[1:-1].replace(',', '') + ')'
        markStr = '(' + str(json.dumps(markList))[1:-1].replace(',', '') + ')'
        return regStr, coreStr, markStr

    def process_base_pairs(self, basespairs:str):
        regBases, coreBases, markList = [], [], []
        for basis_pair in basespairs.replace(', ', ',').split():
            reg, core = basis_pair.replace('(', '').replace(')', '').split(',')
            regBases.append(reg)
            coreBases.append(core)
            markList.append(core)
        return regBases, coreBases, markList
    
    def create_scripts(self, params):
        """Creates a bash script with specified parameters."""
        atom = params['atom']
        mols = params.get('mols')
        molsLoc = params.get('mols(loc)')
        memory = params.get('memory', configs.parameters.GENERAL_RAM_MEMORY)
        time = params.get('time', configs.parameters.GENERAL_TIME_LIMIT)
        extraParams = ''
        for extraParam in {'diis_start', 'diis_space', 'max_scfcycles', 'level_shift'}:
            if extraParam in params:
                extraParams += f"{extraParam}={params[extraParam]}"
        for extraParam in {'atoms_to_localize'}:
            if extraParam in params:
                extraParams += f"{extraParam}=\"{params[extraParam]}\""

        if "basespairs" in params:
            regBases, coreBases, markList = self.process_base_pairs(params['basespairs'])
            extraParams += "doSpecialBasis=true"
        else:
            bases = params.get('bases', configs.parameters.GENERAL_BASES)
            regBases, coreBases, markList = self.process_bases(bases)
        cores = params.get('cores', configs.parameters.GENERAL_CORES)
        partition = params.get('partition', '')

        memToCluster = configs.parameters.MEMORY_TO_CLUSTER
        cluster = memToCluster[max([mem for mem in memToCluster.keys() if int(memory) >= mem])] #select a cluster that corresponds to a largest key s.t. the key is smaller than requirement
        memory = f'{int(memory) * 1000}'  # Convert GB to MB for script compatibility
        partition = f'\n#SBATCH -p {partition}' if cluster == 'telemachus' else ''

        suffix = {True: '_loc', False: ''}
        scripts_directory = Path(f'submitscripts{os.sep}{cluster}')
        scripts_directory.mkdir(parents=True, exist_ok=True)

        for molList, doLoc in ((mols, False), (molsLoc, True)):
            if not molList:
                continue
            script = self.generate_bash_script(atom, time, cores, memory, regBases, coreBases, markList, molList, partition, doLoc, extraParams, suffix[doLoc])
            filename = scripts_directory / f'generated_{atom}_{self.cntr}{suffix[doLoc]}.sh'
            with open(filename, 'w') as file:
                file.write(script)
            self.fnames.add(str(filename))
            self.clusterToFnames.setdefault(cluster, set()).add(str(filename))
            self.cntr += 1

    def generate_bash_script(self, atom, time, cores, memory, regBases, coreBases, markList, molList, partition, doLoc, extraParams, suffix):
        regBases_str = ' '.join(f'"{item}"' for item in regBases)
        coreBases_str = ' '.join(f'"{item}"' for item in coreBases)
        markList_str = ' '.join(f'"{item}"' for item in markList)
        mols_str = ' '.join(f'"{item}"' for item in molList.split())
        localize_flag = 'true' if doLoc else 'false'
        return f"""#!/bin/bash

run_path=$PWD
geom_path="$run_path/mols/exp/{atom}"
out_path="$run_path/mom/{atom}/"

run_time="{time}:00:00"
nproc={cores}
nodes=1
mem={memory}

regBases=({regBases_str})
coreBases=({coreBases_str})
markList=({markList_str})
methods="mom"

mols=({mols_str})

for method in $methods
do
    mkdir -p $out_path
    cd $out_path
    
    for mol in $mols
    do
        mkdir -p $mol 
        cd $mol
                
        nswap=$(sed '2q;d' $geom_path/$mol.xyz | awk '{{print $3}}')
        heavyatom=$(sed '2q;d' $geom_path/$mol.xyz | awk '{{print $2}}')
        for baseIndex in "${{!regBases[@]}}"
        do
            reg=${{regBases[baseIndex]}}
            core=${{coreBases[baseIndex]}}
            mark=${{markList[baseIndex]}}
            zeta=$([ "$mark" = "None" ] && echo "$reg" || echo "$mark")
            mkdir -p $zeta
            cd $zeta

            cat > submit.sh << eof
#!/bin/sh
{partition}
#SBATCH -t $run_time
#SBATCH -o sbatch.out
#SBATCH -e sbatch.err
#SBATCH -c $nproc
#SBATCH -N $nodes
#SBATCH --mem=$mem

export OMP_NUM_THREADS=$nproc
export MKL_NUM_THREADS=$nproc
export OPENBLAS_NUM_THREADS=$nproc
export NUMEXPR_NUM_THREADS=$nproc

python3 $run_path/run_${{method}}.py atom=$heavyatom orbital=$nswap geom=$geom_path/$mol.xyz corebasis=$core regularbasis=$reg localize={localize_flag} {extraParams} > pyscf_output_$method.txt 
eof

            sbatch -J CEBE_${{mol}}_${{zeta}}{suffix} submit.sh
            
            echo $mol $method $zeta done
            cd .. # back to $mol
        done

        cd ..
    done
    cd ..
done
"""
            

#     def create_scripts(self, params):
#         """Creates a bash script with specified parameters

#         Args:
#             params (dict): a parsed dictionary of parameters created by load_parameters.py
#         """
#         atom = params['atom']
#         mols = params.get('mols', None)
#         molsLoc = params.get('mols(loc)', None)
#         memory = params.get('memory', configs.parameters.GENERAL_RAM_MEMORY) # if you don't specify memory, it'll be set to default value specified in configs.parameters
#         time = params.get('time', configs.parameters.GENERAL_TIME_LIMIT)
#         if "basespairs" in params:
#             regBases, coreBases, markList = self.process_base_pairs(params['basespairs']) 
#         else:
#             bases = params.get('bases', configs.parameters.GENERAL_BASES)
#             regBases, coreBases, markList = self.process_bases(bases) 
#         cores = params.get('cores', configs.parameters.GENERAL_CORES)
#         partition = params.get('partition', '')
#         memToCluster = configs.parameters.MEMORY_TO_CLUSTER
#         memToCores = configs.parameters.MEMORY_TO_CORES
        # extraParams = ''
        # for extraParam in {'diis_start', 'diis_space', 'max_scfcycles', 'level_shift'}:
        #     if extraParam in params:
        #         extraParams += f"{extraParam}={params[extraParam]}"
        # for extraParam in {'atoms_to_localize'}:
        #     if extraParam in params:
        #         extraParams += f"{extraParam}=\"{params[extraParam]}\""

        # cluster = memToCluster[max([mem for mem in memToCluster.keys() if int(memory) >= mem])] #select a cluster that corresponds to a largest key s.t. the key is smaller than requirement
#         # cores = memToCores[max([mem for mem in memToCores.keys() if int(memory) >= mem])]
#         memory = f'{memory}000' # convert GB to MB
#         if cluster == 'telemachus':
#             partition = f'\n#SBATCH -p {partition}' # for some reason telemachus requires you to specify a partition if you have a big request
#         else:
#             partition = ''

#         suffix = {True: '_loc', False: ''} # used for job name
        
#         for molList, doLoc in ((mols, False), (molsLoc, True)):
#             if not molList: continue
#             script = '''#!/bin/bash

# run_path=$PWD
# geom_path="$run_path/mols/exp/%s"
# out_path="$run_path/mom/%s/"

# run_time="%s:00:00"
# nproc=%s
# nodes=1
# mem=%s

# regBases=%s
# coreBases=%s
# markList=%s
# methods="mom"

# mols="%s"

# for method in $methods
# do
#     mkdir -p $out_path
#     cd $out_path
    
#         for mol in $mols
#         do
#                 mkdir -p $mol 
#                 cd $mol
                
#         #nelec=$(sed '2q;d' $geom_path/$mol.xyz | awk '{print $1}')
#         #core_atom=$(sed '2q;d' $geom_path/$mol.xyz | awk '{print $2}')
#         nswap=$(sed '2q;d' $geom_path/$mol.xyz | awk '{print $3}')
#         heavyatom=$(sed '2q;d' $geom_path/$mol.xyz | awk '{print $2}')
#         for baseIndex in "${!regBases[@]}"
#         do
#             reg=${regBases[baseIndex]}
#             core=${coreBases[baseIndex]}
#             mark=${markList[baseIndex]}
#         # done
#         # for zeta in $zetas
#         # do
#             # def_basis="cc-pV${zeta}Z"
#             # core_basis="cc-pCV${zeta}Z"
#             if [ $mark = "None" ]; then
#                 zeta=$reg
#             else
#                 zeta=$mark
#             fi
#             mkdir -p $zeta
#             cd $zeta

#                         cat > submit.sh << eof
# #!/bin/sh
# %s
# #SBATCH -t $run_time
# #SBATCH -o sbatch.out
# #SBATCH -e sbatch.err
# #SBATCH -c $nproc
# #SBATCH -N $nodes
# #SBATCH --mem=$mem

# export OMP_NUM_THREADS=$nproc
# export MKL_NUM_THREADS=$nproc
# export OPENBLAS_NUM_THREADS=$nproc
# export NUMEXPR_NUM_THREADS=$nproc

# python3 $run_path/run_${method}.py atom=$heavyatom orbital=$nswap geom=$geom_path/$mol.xyz corebasis=$core regularbasis=$reg localize=%s %s > pyscf_output_$method.txt 
# eof

#                         sbatch -J CEBE_${mol}_${zeta}%s submit.sh
            
#             echo $mol $method $zeta done	
#             cd .. # $mol
#         done

#                 cd ..
#         done
#     cd ..
# done
#     ''' % (atom, atom, time, cores, memory, regBases, coreBases, markList, molList, partition, doLoc, extraParams, suffix[doLoc])
#             with open(fname:=f'submitscripts{os.sep}{cluster}{os.sep}generated_{atom}_{self.cntr}{suffix[doLoc]}.sh', 'w') as f:
#                 f.write(script)
#                 self.fnames.add(fname)
#                 self.clusterToFnames.setdefault(cluster, set()).add(fname)
#                 self.cntr += 1

    def main(self):
        """Main workhorse. 
        """
        # First, ensure there are folders corresponding to each cluster
        for scriptDestination in 'ulysses telemachus'.split():
            if scriptDestination not in os.listdir(f'submitscripts{os.sep}'):
                os.mkdir(f'submitscripts{os.sep}{scriptDestination}/')

        # Then, for each block (dict) of parameters, create a launch script
        for params in self.launchParams:
            self.create_scripts(params)
        
        # for each cluster that will be used, create a "run_cluster.sh" that will specify all names of bash scripts
        for cluster, fnames in self.clusterToFnames.items():
            with open(f'run_{cluster}.sh', 'w') as f:
                f.write('#!/bin/bash\n')
                for fname in fnames:
                    f.write(f"bash {fname}\n")
                    # f.write(f"rm {fname}\n")



if __name__ == "__main__":
    sg = SubmitGenerator()
    sg.main()
    # o = sg.process_bases('D T Q STO-3G 4-31G')
    # print(o)