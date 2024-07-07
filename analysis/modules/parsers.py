from pathlib import Path
from typing import Dict, List, Literal, Optional, Set, TypedDict, cast

import numpy as np
import numpy.typing as npt
from openpyxl import Workbook, load_workbook
from openpyxl.styles import Font

# fmt:off
allowed_extra_keys = {"Frozen orbitals", "Swapped orbitals", "T1 for RHF", "T1 for UHF(a)", "T1 for UHF(b)", "error", "detailed"}
allowed_methods = {"HF", "UHF", "MP2", "MP3", "CCSD", "CCSD(T)", "CCSDT"}
LitKeyType = Literal["Frozen orbitals", "Swapped orbitals", "T1 for RHF", "T1 for UHF(a)", "T1 for UHF(b)", "UHF", "MP2", "MP3", "CCSD", "CCSD(T)", "CCSDT", "error", "detailed"]
# fmt:on
DataPointType = TypedDict(
    "DataPointType",
    {
        "Frozen orbitals": str,
        "Swapped orbitals": str,
        "UHF": float,
        "MP2": float,
        "MP2.5": float,
        "MP3": float,
        "CCSD": float,
        "CCSD(T)": float,
        "CCSDT": float,
        "T1 for RHF": str,
        "T1 for UHF(a)": str,
        "T1 for UHF(b)": str,
        "error": str,
        "detailed": str,
    },
    total=False,
)

ExperDataType = Dict[str, float]
# DataPointType = Dict[str, Union[str, float]]
BasisDataType = Dict[str, DataPointType]
MoleculeDataType = Dict[str, BasisDataType]
AtomDataType = Dict[str, MoleculeDataType]
AlgoDataType = Dict[str, AtomDataType]

MethodErrorType = Dict[str, float]
BasisErrorType = Dict[str, MethodErrorType]
MoleculeErrorType = Dict[str, BasisErrorType]
AtomErrorType = Dict[str, MoleculeErrorType]
AlgoErrorType = Dict[str, AtomErrorType]

MethodErrorsType = Dict[str, List[float]]
BasisErrorsType = Dict[str, MethodErrorsType]
AlgoErrorsType = Dict[str, BasisErrorsType]
AtomErrorsType = Dict[str, BasisErrorsType]
AlgoAtomErrorsType = Dict[str, AtomErrorsType]

# fmt:off
StatsKeyType = Literal["MSE", "MAE", "MedAE", "MaxAE", "STD(AE)", "n", "errors", "abs_errors"]
StatsType = TypedDict('StatsType', {"MSE":float, "MAE":float, "MedAE": float, "MaxAE": float, "STD(AE)":float, "n": int, "errors": npt.NDArray[np.float64], "abs_errors":npt.NDArray[np.float64]})
# fmt:on
# StatsType = Dict[str, Union[float, int]]
MethodStatsType = Dict[str, StatsType]
BasisStatsType = Dict[str, MethodStatsType]
AlgoStatsType = Dict[str, BasisStatsType]
AtomBasisStatsType = Dict[str, BasisStatsType]
AlgoAtomStatsType = Dict[str, AtomBasisStatsType]

# fmt:off
CALC_WB_COLS = {"B": "Molecule","C": "Method","D": "Basis","E": "Atom","F": "UHF","G": "MP2","H": "MP3","I": "CCSD","J": "CCSD(T)","K": "Exper.","L": "Frozen","M": "Swapped","N": "T1 for RHF","O": "T1 for UHF(a)","P": "T1 for UHF(b)","Q": "Error","R": "Comments","S": "Custom Notes"}
# fmt:on


def get_next_col(col: str, prefix: Optional[str] = None) -> str:
    if prefix is None:
        prefix = ""
    strings = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    if not col:
        return "A"
    if len(col) == 1:
        if col == "Z":
            return get_next_col(prefix, prefix=None) + get_next_col("", prefix=None)
        else:
            return prefix + strings[strings.index(col) + 1]
    else:
        return get_next_col(col[1:], prefix=prefix + col[0])


class ParsedData:
    """
    This object takes an excel file with parsed data and updates it with new results. The advantage of this approach is that any custom comments left in column L in CEBE_Data file are preserved.
    """

    # fmt:off
    specialBases = {"ccX-DZ","ccX-TZ","ccX-QZ","ccX-5Z","pcX-1","pcX-2","pcX-3","pcX-4",}  # fmt:on

    def __init__(
        self,
        experimental_wb: str,
        calculations_wb: str,
        algorithms: List[str],
        calculations_folder: Path,
        do_average_experimental: bool = True,
    ):
        """
        Creates a data object. Main reason for using a class is to avoid having to pass different data objects from function to funtion
        """
        self.experWB = load_workbook(experimental_wb, data_only=True)
        self.exper_wb_col = "L"
        self.readWB = load_workbook(calculations_wb)
        self.save_path = calculations_wb
        self.calc_path = calculations_folder

        self.saveWB = Workbook()

        self.molToExper: ExperDataType = {}
        self.doAverageExperimental = do_average_experimental

        self.algorithms = algorithms
        self.algoToData: AlgoDataType = {}

        self.debug = False

        self.font = Font(name="Helvetica", size=12)
        self.cols = "B C D E F G H I J K L M N O P Q R S"
        # algo -> atom -> molecule -> basis -> method -> data (str, float)

    def _parse_calculations(self) -> None:
        """
        Parse data from folders listed in self.algorithmToFolders and store it in a molToData dictionary

        molToData maps molecule (str) to basis (D/T, also str) to a dictionary:
            \{
                "atom": atom that is excited
                "error": error if any
                "UHF": result of UHF calculation
                "MP2": result of MP2 calculation
                "MP3": result of MP3 calculation
                "CCSD": result of CCSD calculation
            \}

        molToData is saved as a value of self.algoToData with relevant key

        P.S. This is already implemented to work with different algorithms (mom, sgm, coreless). Just use a different key in algorithmToFolders
        """
        algoToData = self.algoToData
        for algorithm in self.algorithms:
            algo_path = self.calc_path / algorithm
            atoms = [item.name for item in algo_path.glob("*/")]
            if self.debug:
                print(f"{atoms=}")
            for atom in atoms:
                molToData = algoToData.setdefault(algorithm, {}).setdefault(atom, {})
                if self.debug:
                    print(f"{self.algoToData=}, {molToData=}")
                atom_path = algo_path / atom
                molecules = [molecule.name for molecule in atom_path.glob("*/")]
                if self.debug:
                    print(f"{molecules=}")
                for molecule in molecules:
                    mol_path = atom_path / molecule
                    bases = [basis.name for basis in mol_path.glob("*/")]
                    if self.debug:
                        print(f"{bases=}")
                    for basis in bases:
                        basisToData = molToData.setdefault(molecule, {})
                        data_point = basisToData.setdefault(basis, {})
                        bas_path = mol_path / basis
                        result_file = bas_path / f"CEBE_{algorithm}.txt"

                        # Case 1. Calculation terminated normally
                        if result_file.exists():
                            lines = open(result_file, "r").readlines()
                            for rline in lines:
                                line = rline.split("\n")[0]
                                if ":" in line:
                                    _line_comps = rline.split("\n")[0].split(":")
                                    if _line_comps[0] not in allowed_extra_keys:
                                        raise KeyError(
                                            f"Received unrecognized entry {_line_comps[0]}, expected one of {allowed_extra_keys}"
                                        )
                                    key: LitKeyType = cast(LitKeyType, _line_comps[0])
                                    val: str = _line_comps[1]
                                    data_point[key] = val
                                    continue
                                if "=" not in line:
                                    continue
                                rmethod, rvalue = line.split(" = ")
                                value = float(rvalue.split(" ")[0])
                                posthf_method = cast(
                                    LitKeyType, rmethod.split(" ")[1][1:-1]
                                )
                                data_point[posthf_method] = value
                        # Case 2. Calculation was unsuccesful and did not terminate normally
                        else:
                            data_point["error"] = "No CEBE file"
                            err_file = bas_path / "sbatch.err"
                            if err_file.exists():
                                errFile = open(err_file, "r").readlines()
                                if errFile:
                                    data_point["detailed"] = errFile[-1].split("\n")[0]
                                else:
                                    data_point["detailed"] = "empty?"
                            else:
                                data_point["detailed"] = "empty?"
                                continue

    def _parse_experimental(self) -> None:
        """
        Parse data from experimental_wb file and create a dictionary mapping molecule name (which is used to denote file names) to experimental values
        """
        ws = self.experWB["Sheet1"]
        data = {}
        row = 2
        while True:
            if ws["A" + str(row)].value is None:
                # empty row reached = end of file
                break
            fname = ws["B" + str(row)].value
            data[fname] = ws[self.exper_wb_col + str(row)].value
            row += 1
        self.molToExper = data

    def _calculate_mp25(self) -> None:
        for atomData in self.algoToData.values():
            for molData in atomData.values():
                for basData in molData.values():
                    for data in basData.values():
                        if "MP2" in data and "MP3" in data:
                            if not isinstance(data["MP2"], float):
                                raise TypeError(
                                    f"Expected float for MP2, got {data['MP2']}"
                                )
                            if not isinstance(data["MP3"], float):
                                raise TypeError(
                                    f"Expected float for MP3, got {data['MP3']}"
                                )
                            data["MP2.5"] = (data["MP2"] + data["MP3"]) / 2

    def _update_calculation_wb(self) -> None:
        """
        Write the contents of self.algoToData from parse_results to a CEBE_Data file, preserving any comments in column L (currently not used)
        """
        readWS = self.readWB["Sheet"]
        ws = self.saveWB["Sheet"]

        r = 3
        for col, name in CALC_WB_COLS.items():
            ws[col + str(r)] = name

        for col in self.cols.split(" "):
            ws[col + str(r)].font = self.font
        r += 1
        for algorithm, atomData in self.algoToData.items():
            for atom, molData in atomData.items():
                for molecule, basToData in molData.items():
                    for basis, data in basToData.items():
                        ws["B" + str(r)] = molecule
                        ws["C" + str(r)] = algorithm
                        ws["D" + str(r)] = basis
                        ws["E" + str(r)] = atom
                        if molecule in self.molToExper:
                            ws["K" + str(r)] = self.molToExper[molecule]
                        if "error" in data:
                            ws["Q" + str(r)] = data["error"]
                            ws["R" + str(r)] = data["detailed"]
                        else:
                            if "UHF" in data:
                                ws["F" + str(r)] = data["UHF"]
                            if "MP2" in data:
                                ws["G" + str(r)] = data["MP2"]
                            if "MP3" in data:
                                ws["H" + str(r)] = data["MP3"]
                            if "CCSD" in data:
                                ws["I" + str(r)] = data["CCSD"]
                            if "CCSD(T)" in data:
                                ws["J" + str(r)] = data["CCSD(T)"]
                            if "CCSDT" in data:
                                ws["J" + str(r)] = data["CCSDT"]
                            if "Frozen orbitals" in data:
                                ws["L" + str(r)] = data["Frozen orbitals"]
                            if "Swapped orbitals" in data:
                                ws["M" + str(r)] = data["Swapped orbitals"]
                            if "T1 for RHF" in data:
                                ws["N" + str(r)] = round(float(data["T1 for RHF"]), 4)
                            if "T1 for UHF(a)" in data:
                                ws["O" + str(r)] = round(
                                    float(data["T1 for UHF(a)"]), 4
                                )
                            if "T1 for UHF(b)" in data:
                                ws["P" + str(r)] = round(
                                    float(data["T1 for UHF(b)"]), 4
                                )
                        ws["S" + str(r)] = readWS["S" + str(r)].value
                        for col in self.cols.split(" "):
                            ws[col + str(r)].font = self.font

                        r += 1

        self.saveWB.save(self.save_path)

    def extract_molecules(self, mols: Set[str], save_path: str) -> None:
        """ """
        new_wb = Workbook()
        new_ws = new_wb["Sheet"]

        wb = load_workbook(self.save_path)
        ws = wb["Sheet"]

        row = 3
        for col, name in CALC_WB_COLS.items():
            new_ws[col + str(row)] = name

        for col in self.cols.split(" "):
            new_ws[col + str(3)].font = self.font
        row, new_row = 4, 4
        while True:
            if (mol := ws["B" + str(row)].value) is None:
                break
            if mol in mols:
                new_ws["B" + str(new_row)] = mol
                for col in self.cols.split():
                    new_ws[col + str(new_row)] = ws[col + str(row)].value
                for col in self.cols.split(" "):
                    new_ws[col + str(new_row)].font = self.font
                new_row += 1

            row += 1
        new_wb.save(save_path)

    def filter_data_by_molecules(
        self, algoToData: AlgoDataType, atomToMols: Dict[str, Set[str]]
    ) -> AlgoDataType:
        filtered_algoToData: AlgoDataType = {}
        for algorithm, atomData in algoToData.items():
            for atom, molData in atomData.items():
                if atom not in atomToMols:
                    continue
                for molecule, basData in molData.items():
                    if molecule in atomToMols[atom]:
                        filtered_algoToData.setdefault(algorithm, {}).setdefault(
                            atom, {}
                        )[molecule] = basData
                        filtered_algoToData[algorithm][atom][molecule] = basData
        return filtered_algoToData

    def filter_by_presence_of_experimental(
        self, algoToData: AlgoDataType
    ) -> AlgoDataType:
        filtered_algoToData: AlgoDataType = {}
        for algorithm, atomData in algoToData.items():
            for atom, molData in atomData.items():
                for molecule, basData in molData.items():
                    if molecule in self.molToExper:
                        filtered_algoToData.setdefault(algorithm, {}).setdefault(
                            atom, {}
                        )[molecule] = basData
        return filtered_algoToData

    def calculate_errors(self, algoToData: AlgoDataType) -> None:
        valid_keys = {"UHF", "MP2", "MP2.5", "MP3", "CCSD", "CCSD(T)"}
        algoToError: AlgoErrorType = {}
        for algorithm, atomData in algoToData.items():
            for atom, molData in atomData.items():
                for molecule, basData in molData.items():
                    if molecule not in self.molToExper:
                        continue
                    for bas, methodData in basData.items():
                        for key, value in methodData.items():
                            if key not in valid_keys:
                                continue
                            if not isinstance(value, float):
                                raise TypeError(f"Expected float, got {value}")
                            error = value - self.molToExper[molecule]
                            if self.debug:
                                if error > 1 and bas in self.specialBases:
                                    print(
                                        f"Large error for {algorithm} {atom} {molecule} {bas} {key}: {error}"
                                    )
                            algoToError.setdefault(algorithm, {}).setdefault(
                                atom, {}
                            ).setdefault(molecule, {}).setdefault(bas, {})[key] = error
        self.algoToError = algoToError

    def calculate_series_statistics(self) -> None:
        if not hasattr(self, "algoToError"):
            raise AttributeError("You must call calculate_errors first")
        algoToAtomErrors: AlgoAtomErrorsType = {}
        for algorithm, atomData in self.algoToError.items():
            for atom, molData in atomData.items():
                for basData in molData.values():
                    for bas, methodData in basData.items():
                        for method, value in methodData.items():
                            algoToAtomErrors.setdefault(algorithm, {}).setdefault(
                                atom, {}
                            ).setdefault(bas, {}).setdefault(method, []).append(value)
        algoToAtomStats: AlgoAtomStatsType = {}
        for algorithm, atomErrData in algoToAtomErrors.items():
            for atom, basErrData in atomErrData.items():
                for bas, methodErrData in basErrData.items():
                    for method, error_list in methodErrData.items():
                        errors = np.array(error_list)
                        abs_errs = np.abs(errors)
                        algoToAtomStats.setdefault(algorithm, {}).setdefault(
                            atom, {}
                        ).setdefault(bas, {})[method] = {
                            "MSE": np.mean(errors),
                            "MAE": np.mean(abs_errs),
                            "MedAE": np.median(abs_errs),
                            "MaxAE": np.max(abs_errs),
                            "STD(AE)": np.std(abs_errs),
                            "n": len(errors),
                            "errors": errors,
                            "abs_errors": abs_errs,
                        }
        self.algoToAtomStats = algoToAtomStats

    def calculate_overall_statistics(self) -> None:
        if not hasattr(self, "algoToError"):
            raise AttributeError("You must call calculate_errors first")
        algoToErrors: AlgoErrorsType = {}
        for algorithm, atomData in self.algoToError.items():
            for molData in atomData.values():
                for basData in molData.values():
                    for bas, methodData in basData.items():
                        for method, value in methodData.items():
                            algoToErrors.setdefault(algorithm, {}).setdefault(
                                bas, {}
                            ).setdefault(method, []).append(value)
        algoToStats: AlgoStatsType = {}
        for algorithm, basErrData in algoToErrors.items():
            for bas, methodErrData in basErrData.items():
                for method, error_list in methodErrData.items():
                    errors = np.array(error_list)
                    abs_errs = np.abs(errors)
                    algoToStats.setdefault(algorithm, {}).setdefault(bas, {})[
                        method
                    ] = {
                        "MSE": np.mean(errors),
                        "MAE": np.mean(abs_errs),
                        "MedAE": np.median(abs_errs),
                        "MaxAE": np.max(abs_errs),
                        "STD(AE)": np.std(abs_errs),
                        "n": len(errors),
                        "errors": errors,
                        "abs_errors": abs_errs,
                    }
        self.algoToStats = algoToStats

    def process(self, save: bool = True) -> None:
        self._parse_experimental()
        self._parse_calculations()
        self._calculate_mp25()
        if save:
            self._update_calculation_wb()


if __name__ == "__main__":
    # perform_parsing()
    pass
