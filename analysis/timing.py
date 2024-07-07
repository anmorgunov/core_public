import re
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Set

import numpy as np

from analysis import constants

MethodTimeType = Dict[str, float]
BasisTimeType = Dict[str, MethodTimeType]
MoleculeTimeType = Dict[str, BasisTimeType]


class Timer:
    subtraction = {
        "MP3": "MP2",
        "CCSD(T)": "CCSD",
        "UCCSD(T)": "UCCSD",
        "U-CCSD(T)": "U-CCSD",
    }  # timestamp for completion of MP3 combines MP2 and MP3
    # fmt: off
    newToOld = {"UMP2": "U-MP2", "UCCSD": "U-CCSD","UCCSD(T)": "U-CCSD(T)"}  # fmt:on
    basToMethods = {
        "D": set("HF UHF MP2 UMP2 CCSD UCCSD CCSD(T) UCCSD(T)".split()),
        "T": set("HF UHF MP2 UMP2 CCSD UCCSD CCSD(T) UCCSD(T)".split()),
        "Q": set("HF UHF MP2 UMP2 CCSD UCCSD CCSD(T) UCCSD(T)".split()),
        "5": set("HF UHF MP2 UMP2".split()),
    }
    methods = "HF UHF MP2 UMP2 CCSD UCCSD".split()

    def __init__(self, calc_path: Path, atom_mols: Dict[str, List[str]]) -> None:
        self.calc_path = calc_path
        self.atom_to_mols = atom_mols
        self.molToBasToRun: MoleculeTimeType = {}
        self.mol_bas_cores: Dict[str, Dict[str, int]] = {}
        self.used_cores: Set[int] = set()

        extra_bases = "pcX-1 pcX-2 pcX-3 ccX-DZ ccX-TZ ccX-QZ".split()
        for basis in extra_bases:
            self.basToMethods[basis] = self.basToMethods["Q"]
        extra_pen = "pcX-4 ccX-5Z".split()
        for basis in extra_pen:
            self.basToMethods[basis] = self.basToMethods["5"]

    def parse_timestamps(self) -> None:
        for atom, molecules in self.atom_to_mols.items():
            atomPath = self.calc_path / atom
            for molecule in molecules:
                molPath = atomPath / molecule
                bases = [f.name for f in molPath.glob("*/")]
                for basis in bases:
                    # if timestamps.txt does not exist, continue. Check using path object
                    ts_path = molPath / basis / "timestamps.txt"
                    if not ts_path.exists():
                        continue

                    submit_path = molPath / basis / "submit.sh"
                    with open(submit_path, "r") as f:
                        matches = re.search(r"#SBATCH -c (\d+)", f.read())
                        if matches:
                            self.mol_bas_cores.setdefault(molecule, {})[basis] = int(
                                matches.group(1)
                            )
                            self.used_cores.add(int(matches.group(1)))

                    with open(ts_path, "r") as f:
                        lines = f.readlines()
                    methods = set()
                    for i in range(0, len(lines) - 2, 3):
                        rmethod, rstart, rend = lines[i], lines[i + 1], lines[i + 2]
                        method = rmethod.strip().split("- ")[1]
                        if method in methods:
                            method = f"U{method}"
                        methods.add(method)
                        start = datetime.strptime(rstart.strip(), "%Y-%m-%d %H:%M:%S")
                        end = datetime.strptime(rend.strip(), "%Y-%m-%d %H:%M:%S")
                        runtime = end - start
                        self.molToBasToRun.setdefault(molecule, {}).setdefault(
                            basis, {}
                        )[method] = runtime.seconds

    def export_runtimes(self, save_path: Path, allow_cores: int = 8) -> None:
        body = "# Runtimes\n\n"
        body += "| Basis | RHF | UHF | MP2 | U-MP2 | CCSD | U-CCSD | N |\n"
        body += "| - | - | - | - | - | - | - | - |\n"

        bases = "D pcX-1 ccX-DZ T pcX-2 ccX-TZ Q pcX-3 ccX-QZ 5 pcX-4 ccX-5Z".split()
        for basis in bases:
            row = f"| {basis} |"
            if basis in {"5", "pcX-4", "ccX-5Z"}:
                molecules = [
                    mol
                    for mol in self.molToBasToRun.keys()
                    if basis in self.molToBasToRun[mol]
                    and "UMP2" in self.molToBasToRun[mol][basis]
                ]
            else:
                molecules = [
                    mol
                    for mol in self.molToBasToRun.keys()
                    if basis in self.molToBasToRun[mol]
                    and "UCCSD" in self.molToBasToRun[mol][basis]
                ]
            n_samples_set: Set[int] = set()
            for method in self.methods:
                runtimes_list = []
                for molecule in molecules:
                    # if self.mol_bas_cores[molecule][basis] != allow_cores:
                    #     continue
                    if basis in self.molToBasToRun[molecule]:
                        if method not in self.basToMethods[basis]:
                            continue
                        if method not in self.molToBasToRun[molecule][basis]:
                            if (
                                method not in self.newToOld
                                or self.newToOld[method]
                                not in self.molToBasToRun[molecule][basis]
                            ):
                                continue
                            method = self.newToOld[method]
                        # if basis == '5' and method == 'UHF': print(molecule, method, molToBasToRun[molecule][basis][method])
                        runtime = (
                            self.molToBasToRun[molecule][basis][method]
                            * self.mol_bas_cores[molecule][basis]
                        )
                        if "MP3" in method or "CCSD(T)" in method:
                            runtime -= self.molToBasToRun[molecule][basis][
                                self.subtraction[method]
                            ]
                        runtimes_list.append(runtime)
                if len(runtimes_list) > 0:
                    runtimes = np.array(runtimes_list) / 60
                    row += f" {np.mean(runtimes):.0f} ± {np.std(runtimes):.0f} |"
                    n_samples_set.add(len(runtimes))
                else:
                    row += " ND |"
            if len(n_samples_set) > 1:
                print(
                    f"Warning: {basis} has different number of samples, {n_samples_set}"
                )
            n_samples = n_samples_set.pop() if len(n_samples_set) > 0 else "N/D"
            row += f"{n_samples} |\n"
            body += row
        with open(save_path / "timing-cpu.md", "w") as file:
            file.write(body)

    def main(self, save_path: Path) -> None:
        self.parse_timestamps()
        # for n_cores in self.used_cores:
        # self.export_runtimes(save_path, n_cores)
        self.export_runtimes(save_path)


if __name__ == "__main__":
    data_path = Path(__file__).resolve().parent.parent / "Data"
    t = Timer(
        calc_path=data_path / "calculations" / "mom",
        atom_mols=constants.filtered_atom_to_mols,
    )
    t.main(data_path)
