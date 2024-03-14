from .Modules import LaTeX, Extrapolation, Parser
from .Modules.Parser import MethodDataType, ExperDataType, AtomDataType
from typing import Dict, Set
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

    def __init__(
        self, atomToData: AtomDataType, molToExper: ExperDataType
    ):
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
        label = f"{atom.lower()}-{basis.lower()}z"

        LaTeX.Table.export_table(
            caption=caption,
            label=label,
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
