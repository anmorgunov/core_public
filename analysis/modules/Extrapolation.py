import plotly.graph_objects as go
import scipy.optimize
import numpy as np
from .Parser import AlgoDataType, BasisDataType, ExperDataType
from typing import List, Tuple, Union, Dict, Optional

Number = Union[float, int]


def rounder(dig):
    def rounder_to_dig(float):
        return format(np.around(float, dig), f".{dig}f")

    return rounder_to_dig


def _helgaker(zeta, a, b):
    """A standard helgaker extrapolation formula. Note: a, b are parameters which we seek. We pass them as parameters not because they're known,
    but because that's how scipy.optimize works (which actually finds those parameters).

    Args:
        zeta (int): an integer describing "zeta" value of basis set
        a (float): extrapolation parameter
        b (float): extrapolation parameter

    Returns:
        float: a fitted value at a given zeta level and given parameters
    """
    return a + b * (zeta ** (-3))


def extrapolate_energies(zetas, values):
    popt, pcov = scipy.optimize.curve_fit(_helgaker, zetas, values)
    return popt[0]


def parse_scheme(scheme: str) -> Tuple[str, List[str], str, bool]:
    if "+" in scheme:
        preargs, corr = scheme.split("+")
        args = preargs.replace("[", " ").replace("]", " ").split()
        method = args[0]
        bases = args[1:]
        corrBasis = corr.split("Dif")[-1]
        if "(T)" in corrBasis:
            corrBasis = corrBasis.split("(T)")[0]
            corrTriples = True
        else:
            corrTriples = False
    elif "[" in scheme:
        args = scheme.replace("[", " ").replace("]", " ").split()
        method = args[0]
        bases = args[1:]
        corrBasis = None
        corrTriples = None
    elif "-" in scheme:
        args = scheme.split("-")
        method = args[-1]
        bases = args[:-1]
        corrBasis = None
        corrTriples = None
    else:
        raise NameError(f"Scheme {scheme} is not recognized")
    return method, bases, corrBasis, corrTriples


def extrapolate_molecule_given_scheme(
    basisData: BasisDataType, scheme: str
) -> Optional[Dict[str, Optional[Number]]]:
    method, bases, corrBasis, corrTriples = parse_scheme(scheme)
    if method == "HF":
        method = "UHF"

    basToCoeff = {"D": 2, "T": 3, "Q": 4, "5": 5}
    if len(bases) == 1:
        assert (
            corrBasis is not None
        ), "You need to specify at least two basis sets for extrapolation"

        method_cbs = basisData[bases[0]][method]
    else:
        energies = [basisData[basis][method] for basis in bases]
        zetas = [basToCoeff[basis] for basis in bases]

        method_cbs = extrapolate_energies(zetas, energies)

    if corrBasis is not None:
        corrMethod = "CCSD(T)" if corrTriples else "CCSD"
        if (
            corrBasis not in basisData
            or "MP2" not in basisData[corrBasis]
            or corrMethod not in basisData[corrBasis]
        ):
            return None

        mp2 = basisData[corrBasis]["MP2"]
        cc = basisData[corrBasis][corrMethod]
        corr = cc - mp2
        return {"cbs": method_cbs, "cbs+corr": method_cbs + corr, "corr": corr}
    else:
        return {
            "cbs": method_cbs,
            "cbs+corr": None,
            "corr": None,
        }


class WholeDataset:
    """An object that is used to extrapolate results to CBS Limit for a given series.
    As of now, it assumes you have CCSD/CCSD(T) results for D,T,Q basis sets AND MP2 results for D,T,Q,5
    """

    def __init__(self):
        """
        Initialization
        """
        self.algoToCBS = {}
        self.smallBasisException = {}
        self.debug = False

    def extrapolate_all_data(self, algoData: AlgoDataType, schemes: List[str]):
        for algorithm, atomData in algoData.items():
            for atom, molData in atomData.items():
                for mol, basisData in molData.items():
                    for scheme in schemes:
                        result = extrapolate_molecule_given_scheme(
                            basisData, scheme
                        )
                        if result is None:
                            self.smallBasisException.setdefault(algorithm, {}).setdefault(
                                atom, {}
                            )[mol] = scheme
                        else:
                            self.algoToCBS.setdefault(algorithm, {}).setdefault(
                                atom, {}
                            ).setdefault(mol, {})[scheme] = result

    def calculate_errors(self, experimentalData: ExperDataType):
        algoToError = {}
        for algorithm, atomData in self.algoToCBS.items():
            for atom, molData in atomData.items():
                for mol, schemeData in molData.items():
                    for scheme, data in schemeData.items():
                        if data["cbs+corr"] is not None:
                            error = data["cbs+corr"] - experimentalData[mol]
                            algoToError.setdefault(algorithm, {}).setdefault(
                                atom, {}
                            ).setdefault(mol, {})[scheme] = error
                        else:
                            error = data["cbs"] - experimentalData[mol]
                            algoToError.setdefault(algorithm, {}).setdefault(
                                atom, {}
                            ).setdefault(mol, {})[scheme] = error
        self.algoToError = algoToError

    def calculate_series_statistics(self):
        assert hasattr(self, "algoToError"), "You need to call calculate_errors first"
        algoToAtomErrors = {}
        for algorithm, atomData in self.algoToError.items():
            for atom, molData in atomData.items():
                for schemeData in molData.values():
                    for scheme, error in schemeData.items():
                        algoToAtomErrors.setdefault(algorithm, {}).setdefault(
                            atom, {}
                        ).setdefault(scheme, []).append(error)
        algoToAtomStats = {}
        for algorithm, atomData in algoToAtomErrors.items():
            for atom, schemeData in atomData.items():
                for scheme, errors_list in schemeData.items():
                    errors = np.array(errors_list)
                    abs_errs = np.abs(errors)
                    algoToAtomStats.setdefault(algorithm, {}).setdefault(atom, {})[
                        scheme
                    ] = {
                        "MSE": np.mean(errors),
                        "MAE": np.mean(abs_errs),
                        "MedAE": np.median(abs_errs),
                        "MaxAE": np.max(abs_errs),
                        "STD(AE)": np.std(abs_errs),
                        "n": len(errors),
                    }
        self.algoToAtomStats = algoToAtomStats

    def calculate_overall_statistics(self):
        assert hasattr(self, "algoToError"), "You need to call calculate_errors first"
        algoToErrors = {}
        for algorithm, atomData in self.algoToError.items():
            for molData in atomData.values():
                for schemeData in molData.values():
                    for scheme, error in schemeData.items():
                        algoToErrors.setdefault(algorithm, {}).setdefault(
                            scheme, []
                        ).append(error)
        algoToStats = {}
        for algorithm, schemeData in algoToErrors.items():
            for scheme, errors in schemeData.items():
                errors = np.array(errors)
                abs_errs = np.abs(errors)
                algoToStats.setdefault(algorithm, {})[scheme] = {
                    "MSE": np.mean(errors),
                    "MAE": np.mean(abs_errs),
                    "MedAE": np.median(abs_errs),
                    "MaxAE": np.max(abs_errs),
                    "STD(AE)": np.std(abs_errs),
                    "n": len(errors),
                }
        self.algoToStats = algoToStats

    # def export_errors(self, names, atoms):
    #     """Exports errors as md table.

    #     The names used define what data will be printed in the table.
    #     """
    #     atomstr = "".join(atoms)
    #     self.fill_dictionaries(names, atoms)
    #     headers = [
    #         "Extrapolation Scheme",
    #         "MSE",
    #         "MAE",
    #         "MedAE",
    #         "MaxAE",
    #         "STD",
    #         "Sample Size",
    #     ]
    #     body = []
    #     for name in names:
    #         erD = self.nameToErrors[name]
    #         row = [
    #             name,
    #             f"{erD['MSE']}",
    #             f"{erD['MAE']}",
    #             f"{erD['MedAE']}",
    #             f"{erD['MaxAE']}",
    #             f"{erD['STD(AE)']}",
    #             f"{len(erD['raw'])}",
    #         ]
    #         body.append(row)
    #     return headers, body

    def _create_extrapolation_schemes(self):
        # CCSD schemes
        ccsd_schemes = []
        for method in "CCSD CCSD(T)".split():
            for bases in "D-T T-Q D-T-Q".split():
                ccsd_schemes.append(f"{bases}-{method}")
        if self.debug:
            print(f"{ccsd_schemes=}")

        hf_schemes = []
        for bases in "T-Q D-T-Q T-Q-5 D-T-Q-5".split():
            hf_schemes.append(f"{bases}-HF")
        if self.debug:
            print(f"{hf_schemes=}")

        mp2_schemes = []
        mp2_bases = "D | T | Q | 5 | D T | D T Q | T Q | Q 5 | T Q 5 | D T Q 5"
        for bases in mp2_bases.split(" | "):
            if len(bases) != 1:
                mp2_schemes.append(f"MP2[{bases}]")
            for (
                small_basis
            ) in "STO-3G STO-6G 3-21G 4-31G 6-31G def2svp def2svpd D".split():
                if bases == "D" and small_basis == "D":
                    continue
                mp2_schemes.append(f"MP2[{bases}]+Dif{small_basis}")
                mp2_schemes.append(f"MP2[{bases}]+Dif{small_basis}(T)")

        if self.debug:
            print(f"{mp2_schemes=}")

        self.schemes = {
            "CCSD": ccsd_schemes,
            "HF": hf_schemes,
            "MP2": mp2_schemes,
        }


if __name__ == "__main__":
    atomToNames = {
        "C": "D-T-Q-CCSD(T), MP2(D T Q)+DifD(T)".split(", "),
        "N": "D-T-Q-CCSD(T), MP2(D T Q 5)+DifD(T)".split(", "),
        "O": "T-Q-CCSD, MP2(Q)+DifD".split(", "),
        "F": "T-Q-CCSD(T), MP2(Q)+DifD".split(", "),
    }

    obj = WholeDataset()