import os
import sys
sys.path.append('../core_excitations') # this allows us to import files from parent directories
import PRIVATE
import tools
from datetime import datetime
import numpy as np

class Timer:

    def __init__(self, methodToFolders):

        PATH = PRIVATE.BASE_PATH
        self.BASE = PATH
        self.SUBMODULE = self.BASE + ["analysis", "parser"]
        self.CALCULATIONS = self.BASE + ["calculations"]

        self.methodToFolders = methodToFolders
        self.molToBasToRun = {}

    def parse_timestamps(self):
        methodToFolders = self.methodToFolders
        molToBasToRun = self.molToBasToRun
        for algorithm, folders in methodToFolders.items():
            for folder in folders:
                mefolder = tools.join_path([algorithm, folder])

                molPath = tools.join_path(self.CALCULATIONS+[mefolder,], endWithSlash=True)
                molecules = [molecule for molecule in os.listdir(molPath) if os.path.isdir(tools.append_path(molPath, [molecule,], endWithSlash=True))]
                for molecule in molecules: 
                    if molecule not in molToBasToRun:
                        molToBasToRun[molecule] = {}

                    bases = [basis for basis in os.listdir(tools.append_path(molPath, [molecule,], endWithSlash=True)) if os.path.isdir(tools.append_path(molPath, [molecule, basis], endWithSlash=True))]
                    for basis in bases:
                        if basis not in molToBasToRun[molecule]:
                            molToBasToRun[molecule][basis] = {}
                        with open(tools.append_path(molPath, [molecule, basis, f"timestamps.txt"]), 'r') as f:
                            lines = f.readlines()
                            # print(molPath, molecule, basis, lines )
                        methods = set()
                        for i in range(0, len(lines)-2, 3):
                            rmethod, rstart, rend = lines[i], lines[i+1], lines[i+2]
                            method = rmethod.strip().split('- ')[1]
                            if method in methods:
                                method = f'U{method}'
                            methods.add(method)
                            start = datetime.strptime(rstart.strip(), "%Y-%m-%d %H:%M:%S")
                            end = datetime.strptime(rend.strip(), "%Y-%m-%d %H:%M:%S")
                            runtime = end-start
                            # print(method, start, end)
                            # print(molecule, basis, method)
                            molToBasToRun[molecule][basis][method] = runtime.seconds


                        #     if ':' in line:
                        #         key, val = rline.split("\n")[0].split(":")
                        #         molToData[molecule][basis][key] = val
                        #         continue
                        #     if "=" not in line: continue
                        #     rmethod, rvalue = line.split(' = ')
                        #     value = float(rvalue.split(' ')[0])
                        #     posthf_method = rmethod.split(' ')[1][1:-1]
                        #     molToData[molecule][basis][posthf_method] = value
    def export_errors(self):
        molToBasToRun = self.molToBasToRun
        # print(len(molToBasToRun))

        f = open(tools.join_path(self.BASE+ ["analysis", "extrapolation", "timer.md"]), 'w')
        f.write('# Runtimes\n\n')
        f.write("| Stat | Basis | RHF | UHF | MP2 | U-MP2 | MP3 | U-MP3 | CCSD | U-CCSD | (T) | U-(T) |\n")
        f.write("| - | - | - | - | - | - | - | - | - | - | - | - |\n")


        subtraction = {'MP3': 'MP2', 'UMP3': 'UMP2', 'CCSD(T)': 'CCSD', 'UCCSD(T)': 'UCCSD', 'U-MP3': 'U-MP2', 'U-CCSD(T)': 'U-CCSD'}
        newToOld = {'UMP2': 'U-MP2', 'UMP3': 'U-MP3', 'UCCSD': 'U-CCSD', 'UCCSD(T)': 'U-CCSD(T)'}
        basToMethods = {
            'D': set('HF UHF MP2 UMP2 MP3 UMP3 CCSD UCCSD CCSD(T) UCCSD(T)'.split()),
            'T': set('HF UHF MP2 UMP2 MP3 UMP3 CCSD UCCSD CCSD(T) UCCSD(T)'.split()),
            'Q': set('HF UHF MP2 UMP2 CCSD UCCSD CCSD(T) UCCSD(T)'.split()),
            '5': set('HF UHF MP2 UMP2'.split())
        }
        # print([molecule for molecule in molToBasToRun.keys() if 'CCSD' in molToBasToRun[molecule]['5']])

        for basis in 'D T Q 5'.split():
            avgRow, stdRow, medRow, maxRow = f'| Avg | {basis} |', f'| Std | {basis} |', f'| Med | {basis} |', f'| Max | {basis} |'
            if basis == '5':
                molecules = [mol for mol in molToBasToRun.keys() if 'UMP2' in molToBasToRun[mol][basis]]
            else:
                molecules = [mol for mol in molToBasToRun.keys() if 'UCCSD' in molToBasToRun[mol][basis]]

            for method in 'HF UHF MP2 UMP2 MP3 UMP3 CCSD UCCSD CCSD(T) UCCSD(T)'.split():
                runtimes = []
                for molecule in molecules:
                    if basis in molToBasToRun[molecule]:
                        if method not in basToMethods[basis]: continue
                        if method not in molToBasToRun[molecule][basis]:
                            if method not in newToOld or newToOld[method] not in molToBasToRun[molecule][basis]: continue
                            method = newToOld[method]
                        # if basis == '5' and method == 'UHF': print(molecule, method, molToBasToRun[molecule][basis][method])
                        runtime = molToBasToRun[molecule][basis][method]
                        if 'MP3' in method or 'CCSD(T)' in method:
                            if basis == '5': print(method, runtime)
                            runtime -= molToBasToRun[molecule][basis][subtraction[method]]
                        runtimes.append(runtime)
                if runtimes:
                    # print(basis, len(runtimes))
                    avg = round(np.mean(runtimes)/60, 1)
                    std = round(np.std(runtimes)/60, 1)
                    med = round(np.median(runtimes)/60, 1)
                    maxV = round(np.max(runtimes)/60, 1)
                else:
                    avg, std, med, maxV = 'ND', 'ND', 'ND', 'ND'
                avgRow += f" {avg} |"
                stdRow += f" {std} |"
                medRow += f" {med} |"
                maxRow += f" {maxV} |"
            avgRow += '\n'
            stdRow += '\n'
            medRow += '\n'
            maxRow += '\n'
            f.write(avgRow)
            f.write(stdRow)
            f.write(medRow)
            f.write(maxRow)

    def main(self):
        self.parse_timestamps()
        self.export_errors()

if __name__ == "__main__":
    methodToFolders = {'mom': ['o', 'f', 'c', 'cl', 'n', 's']}

    tobj = Timer(methodToFolders)
    tobj.main()