from datetime import datetime
from pathlib import Path
from typing import Dict, Set

import numpy as np

MethodTimeType = Dict[str, float]
BasisTimeType = Dict[str, MethodTimeType]
MoleculeTimeType = Dict[str, BasisTimeType]


class Timer:
    subtraction = {
        "MP3": "MP2",
        "UMP3": "UMP2",
        "CCSD(T)": "CCSD",
        "UCCSD(T)": "UCCSD",
        "U-MP3": "U-MP2",
        "U-CCSD(T)": "U-CCSD",
    }  # timestamp for completion of MP3 combines MP2 and MP3
    # fmt: off
    newToOld = {"UMP2": "U-MP2","UMP3": "U-MP3","UCCSD": "U-CCSD","UCCSD(T)": "U-CCSD(T)"}  # fmt:on
    basToMethods = {
        "D": set("HF UHF MP2 UMP2 MP3 UMP3 CCSD UCCSD CCSD(T) UCCSD(T)".split()),
        "T": set("HF UHF MP2 UMP2 MP3 UMP3 CCSD UCCSD CCSD(T) UCCSD(T)".split()),
        "Q": set("HF UHF MP2 UMP2 CCSD UCCSD CCSD(T) UCCSD(T)".split()),
        "5": set("HF UHF MP2 UMP2".split()),
    }
    methods = "HF UHF MP2 UMP2 MP3 UMP3 CCSD UCCSD CCSD(T) UCCSD(T)".split()

    def __init__(self, calc_path: Path, atoms: Set[str]) -> None:
        self.calc_path = calc_path
        self.desired_atoms = atoms
        self.molToBasToRun: MoleculeTimeType = {}

        extra_bases = "pcX-1 pcX-2 pcX-3 ccX-DZ ccX-TZ ccX-QZ".split()
        for basis in extra_bases:
            self.basToMethods[basis] = self.basToMethods["Q"]
        extra_pen = "pcX-4 ccX-5Z".split()
        for basis in extra_pen:
            self.basToMethods[basis] = self.basToMethods["5"]

    def parse_timestamps(self) -> None:
        atoms = [f.name for f in self.calc_path.glob("*/")]
        for atom in atoms:
            if atom not in self.desired_atoms:
                continue
            atomPath = self.calc_path / atom
            molecules = [f.name for f in atomPath.glob("*/")]
            for molecule in molecules:
                molPath = atomPath / molecule
                bases = [f.name for f in molPath.glob("*/")]
                for basis in bases:
                    # if timestamps.txt does not exist, continue. Check using path object
                    ts_path = molPath / basis / "timestamps.txt"
                    if not ts_path.exists():
                        continue
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

    def export_errors(self, save_path: Path) -> None:
        body = "# Runtimes\n\n"
        body += "| Stat | Basis | RHF | UHF | MP2 | U-MP2 | MP3 | U-MP3 | CCSD | U-CCSD | (T) | U-(T) |\n"
        body += "| - | - | - | - | - | - | - | - | - | - | - | - |\n"

        bases = "D pcX-1 ccX-DZ T pcX-2 ccX-TZ Q pcX-3 ccX-QZ 5 pcX-4 ccX-5Z".split()
        for basis in bases:
            avgRow = f"| Avg | {basis} |"
            stdRow = f"| Std | {basis} |"
            medRow = f"| Med | {basis} |"
            maxRow = f"| Max | {basis} |"
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

            for method in self.methods:
                runtimes = []
                for molecule in molecules:
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
                        runtime = self.molToBasToRun[molecule][basis][method]
                        if "MP3" in method or "CCSD(T)" in method:
                            runtime -= self.molToBasToRun[molecule][basis][
                                self.subtraction[method]
                            ]
                        runtimes.append(runtime)
                if runtimes:
                    # print(basis, len(runtimes))
                    avg = f"{(np.mean(runtimes) / 60):.1f}"
                    std = f"{(np.std(runtimes) / 60):.1f}"
                    med = f"{(np.median(runtimes) / 60):.1f}"
                    maxV = f"{(np.max(runtimes) / 60):.1f}"
                else:
                    avg, std, med, maxV = "ND", "ND", "ND", "ND"
                avgRow += f" {avg} |"
                stdRow += f" {std} |"
                medRow += f" {med} |"
                maxRow += f" {maxV} |"
            avgRow += "\n"
            stdRow += "\n"
            medRow += "\n"
            maxRow += "\n"
            body += avgRow
            # body += stdRow
            # body += medRow
            body += maxRow
        with open(save_path / "timer.md", "w") as file:
            file.write(body)

    def main(self, save_path: Path) -> None:
        self.parse_timestamps()
        self.export_errors(save_path)


if __name__ == "__main__":
    pass
