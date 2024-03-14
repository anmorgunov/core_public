import os
from datetime import datetime
from typing import Set, Dict
import numpy as np 

MethodTimeType = Dict[str, float]
BasisTimeType = Dict[str, MethodTimeType]
MoleculeTimeType = Dict[str, BasisTimeType]


class Timer:

    def __init__(self, calc_path:str, atoms:Set[str]):
        self.calc_path = calc_path
        self.atoms = atoms
        self.molToBasToRun:MoleculeTimeType = {}

    def parse_timestamps(self):
        molToBasToRun = self.molToBasToRun
        for atom in os.listdir(self.calc_path):
            if atom not in self.atoms: continue
            atomPath = os.path.join(self.calc_path, atom)
            molDirs = [molecule for molecule in os.listdir(atomPath) if os.path.isdir(os.path.join(atomPath, molecule)) and molecule != '.DS_Store']
            for molecule in molDirs:
                molPath = os.path.join(atomPath, molecule)
                basDirs = [basis for basis in os.listdir(molPath) if os.path.isdir(os.path.join(molPath, basis)) and basis != '.DS_Store']
                for basis in basDirs:
                    with open(os.path.join(molPath, basis, f"timestamps.txt"), 'r') as f:
                        lines = f.readlines()
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
                        molToBasToRun.setdefault(molecule, {}).setdefault(basis, {})[method] = runtime.seconds
        
    def export_errors(self, save_path:str):
        molToBasToRun = self.molToBasToRun
        # print(len(molToBasToRun))

        f = open(os.path.join(save_path, "timer.md"), 'w')
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
                molecules = [mol for mol in molToBasToRun.keys() if basis in molToBasToRun[mol] and 'UMP2' in molToBasToRun[mol][basis]]
            else:
                molecules = [mol for mol in molToBasToRun.keys() if basis in molToBasToRun[mol] and 'UCCSD' in molToBasToRun[mol][basis]]

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

    def main(self, save_path:str):
        self.parse_timestamps()
        self.export_errors(save_path)

if __name__ == "__main__":

    tobj = Timer()
    tobj.main()