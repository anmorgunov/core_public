from .Modules import LaTeX, Extrapolation, Parser
from .Modules.Parser import (
    ExperDataType,
    AtomDataType,
    BasisStatsType,
    AtomBasisStatsType,
)
from typing import Dict, Set, List
import os
from . import constants

BASE_PATH = os.path.dirname(__file__)


class SingleZetaResults:
    basisToHeader = {
        "5": ["Molecule", "$\Delta$HF", "$\Delta$MP2", "Experiment"],
        "Q": [
            "Molecule",
            "$\Delta$HF",
            "$\Delta$MP2",
            "$\Delta$CCSD",
            "$\Delta$CCSD(T)",
            "Experiment",
        ],
        "T": [
            "Molecule",
            "$\Delta$HF",
            "$\Delta$MP2",
            "$\Delta$MP3",
            "$\Delta$CCSD",
            "$\Delta$CCSD(T)",
            "Experiment",
        ],
        "D": [
            "Molecule",
            "$\Delta$HF",
            "$\Delta$MP2",
            "$\Delta$MP3",
            "$\Delta$CCSD",
            "$\Delta$CCSD(T)",
            "Experiment",
        ],
    }

    def __init__(self, atomToData: AtomDataType, molToExper: ExperDataType):
        self.atomToData = atomToData
        self.molToExper = molToExper
        self.f = lambda i: f"{i:.2f}"

    def _render_mol(self, mol: str):
        return constants.FNAME_TO_MOLS[mol]["latex"]

    def collect_series(self, atom: str, basis: str):
        body = []
        f = self.f  # f stands for formatter
        for molecule in sorted(self.atomToData[atom]):
            _data = self.atomToData[atom][molecule][basis]
            exp = self.molToExper[molecule]
            uhf, mp2 = _data["UHF"], _data["MP2"]
            row = ["\ch{%s}" % self._render_mol(molecule)]
            if basis == "5":
                energies = [uhf, mp2, exp]
            elif basis == "Q":
                ccsd, ccsdt = _data["CCSD"], _data["CCSD(T)"]
                energies = [uhf, mp2, ccsd, ccsdt, exp]
            else:
                ccsd, ccsdt = _data["CCSD"], _data["CCSD(T)"]
                mp3 = _data["MP3"]
                energies = [uhf, mp2, mp3, ccsd, ccsdt, exp]

            for e in energies:
                row.append(f(e))
            body.append(row)
        return body

    def series_table(self, atom: str, basis: str, save_path: str):
        body = self.collect_series(atom, basis)
        caption = f"K-Edge Ionization Energies (in eV) of {atom}-Series in cc-pV{basis}Z/cc-pCV{basis}Z"

        LaTeX.Table.export_table(
            caption=caption,
            label=f"{atom.lower()}-{basis.lower()}z",
            positioning="l " * len(self.basisToHeader[basis]),
            headers=self.basisToHeader[basis],
            body=body,
            name=save_path,
        )

    def all_results(self, save_folder: str):
        for atom in "C N O F".split():
            for basis in "D T Q 5".split():
                path = os.path.join(save_folder, f"{atom}-{basis}Z")
                self.series_table(atom.lower(), basis, path)


class MethodSummary:
    col_names = ["Basis", "Method", "MSE", "MAE", "MedAE", "MaxAE", "STD"]
    basisToMethods = {
        "D": "UHF MP2 MP2.5 MP3 CCSD CCSD(T)".split(),
        "T": "UHF MP2 MP2.5 MP3 CCSD CCSD(T)".split(),
        "Q": "UHF MP2 CCSD CCSD(T)".split(),
        "5": "UHF MP2".split(),
    }

    def __init__(
        self,
        basisToStats: BasisStatsType,
        atomToBasisStats: AtomBasisStatsType,
        save_folder: str,
        show_sample_size: bool = True,
        isPublication: bool = False,
    ):
        self.basisToStats = basisToStats
        self.atomToBasisStats = atomToBasisStats
        self.f = lambda i: f"{i:.2f}"
        self.save_folder = save_folder
        self.show_sample_size = show_sample_size
        self.isPublication = isPublication

    def collect_series(self, atom: str, bases: List[str]):
        _me = lambda x: f"$\Delta${x}" if x != "UHF" else f"$\Delta$HF"
        body = []
        f = self.f  # f stands for formatter
        for basis in bases:
            if self.isPublication:
                methods = self.basisToMethods[basis]
            else:
                methods = self.atomToBasisStats[atom][basis]
            for method in methods:
                stats = self.atomToBasisStats[atom][basis][method]
                row = [
                    basis,
                    _me(method),
                    f(stats["MSE"]),
                    f(stats["MAE"]),
                    f(stats["MedAE"]),
                    f(stats["MaxAE"]),
                    f(stats["STD(AE)"]),
                ]
                if self.show_sample_size:
                    row.append(stats["n"])
                body.append(row)
        return body

    def series_table(self, atom: str, bases: List[str], save_path: str):
        body = self.collect_series(atom, bases)

        if self.show_sample_size:
            self.col_names.append("Sample Size")
        header = ["\\textbf{%s}" % col_name for col_name in self.col_names]
        caption = f"Statistical analysis of accuracy of different methods at predicting K-Edge CEBEs (in eV) compared to experimental data for {atom}-series"
        LaTeX.Table.export_table(
            caption=caption,
            label=f"summary-{atom.lower()}",
            positioning="l " * len(header),
            headers=header,
            body=body,
            name=save_path,
        )

    def all_results(
        self,
    ):
        atoms = "C N O F".split()
        bases = "D T Q 5".split()
        for atom in atoms:
            path = os.path.join(self.save_folder, f"{atom}-summary")
            self.series_table(atom.lower(), bases, path)
