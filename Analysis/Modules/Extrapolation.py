from typing import Dict, Iterable, List, Literal, Optional, Tuple, TypedDict, cast

import numpy as np
import numpy.typing as npt
import scipy.optimize

from Analysis.Modules.Parser import DataPointType, LitKeyType, allowed_methods

# fmt:off
StatsKeyType = Literal["MSE", "MAE", "MedAE", "MaxAE", "STD(AE)", "n", "errors", "abs_errors"]
StatsType = TypedDict('StatsType', {"MSE":float, "MAE":float, "MedAE": float, "MaxAE": float, "STD(AE)":float, "n": int, "errors": npt.NDArray[np.float64], "abs_errors":npt.NDArray[np.float64]})
# fmt:on
ExperDataType = Dict[str, float]
# DataPointType = Dict[str, Union[str, float]]
BasisDataType = Dict[str, DataPointType]
MoleculeDataType = Dict[str, BasisDataType]
AtomDataType = Dict[str, MoleculeDataType]
AlgoDataType = Dict[str, AtomDataType]

SchemeResultType = Dict[str, Optional[float]]
SchemesCBSType = Dict[str, SchemeResultType]
MolCBSType = Dict[str, SchemesCBSType]
AtomCBSType = Dict[str, MolCBSType]
AlgoCBSType = Dict[str, AtomCBSType]

SchemeErrType = Dict[str, float]
MolErrType = Dict[str, SchemeErrType]
AtomErrType = Dict[str, MolErrType]
AlgoErrType = Dict[str, AtomErrType]

SchemeErrsType = Dict[str, List[float]]
AlgoErrsType = Dict[str, SchemeErrsType]
AtomErrsType = Dict[str, SchemeErrsType]
AlgoAtomErrsType = Dict[str, AtomErrsType]

SchemeStatsType = Dict[str, StatsType]
AlgoStatsType = Dict[str, SchemeStatsType]
AtomStatsType = Dict[str, SchemeStatsType]
AlgoAtomStatsType = Dict[str, AtomStatsType]


def _helgaker(zeta: int, a: float, b: float) -> float:
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


def extrapolate_energies(zetas: List[int], values: List[float]) -> float:
    popt, pcov = scipy.optimize.curve_fit(_helgaker, zetas, values)
    return cast(float, popt[0])


def parse_scheme(scheme: str) -> Tuple[str, List[str], str, Optional[bool]]:
    if "+" in scheme:
        preargs, corr = scheme.split("+")
        args = preargs.replace("[", " ").replace("]", " ").split()
        if args[0] not in allowed_methods:
            raise KeyError(f"Received unrecognized method {args[0]}")
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
        if args[0] not in allowed_methods:
            raise KeyError(f"Received unrecognized method {args[0]}")
        method = args[0]
        bases = args[1:]
        corrBasis = None
        corrTriples = None
    elif "-" in scheme:
        args = scheme.split("-")
        if args[-1] not in allowed_methods:
            raise KeyError(f"Received unrecognized method {args[0]}")
        method = args[-1]
        bases = args[:-1]
        corrBasis = None
        corrTriples = None
    else:
        raise NameError(f"Scheme {scheme} is not recognized")
    return method, bases, corrBasis, corrTriples


def extrapolate_molecule_given_scheme(
    basisData: BasisDataType, scheme: str
) -> Optional[SchemeResultType]:
    method, bases, corrBasis, corrTriples = parse_scheme(scheme)
    if method == "HF":
        method = "UHF"
    method = cast(LitKeyType, method)

    # fmt:off
    basToCoeff = {"D": 2, "T": 3, "Q": 4, "5": 5, "pcX-1": 2, "pcX-2": 3, "pcX-3": 4, "pcX-4": 5, "ccX-DZ": 2, "ccX-TZ": 3, "ccX-QZ": 4, "ccX-5Z": 5,}  # fmt:on
    if len(bases) == 1:
        if corrBasis is None:
            raise TypeError(
                "You need to specify at least two basis sets for extrapolation"
            )
        method_cbs = cast(float, basisData[bases[0]][method])
    else:
        energies = [cast(float, basisData[basis][method]) for basis in bases]
        zetas = [basToCoeff[basis] for basis in bases]

        method_cbs = extrapolate_energies(zetas, energies)

    if corrBasis is None:
        return {
            "cbs": method_cbs,
            "cbs+corr": None,
            "corr": None,
        }
    else:
        corrMethod: LitKeyType = "CCSD(T)" if corrTriples else "CCSD"
        if (
            corrBasis not in basisData
            or "MP2" not in basisData[corrBasis]
            or corrMethod not in basisData[corrBasis]
        ):
            return None

        mp2 = basisData[corrBasis]["MP2"]
        cc = cast(float, basisData[corrBasis][corrMethod])
        corr = cc - mp2
        return {"cbs": method_cbs, "cbs+corr": method_cbs + corr, "corr": corr}


class WholeDataset:
    """An object that is used to extrapolate results to CBS Limit for a given series.
    As of now, it assumes you have CCSD/CCSD(T) results for D,T,Q basis sets AND MP2 results for D,T,Q,5
    """

    def __init__(self) -> None:
        """
        Initialization
        """
        self.algoToCBS: AlgoCBSType = {}
        self.smallBasisException: Dict[str, Dict[str, Dict[str, str]]] = {}
        self.debug = False

    def extrapolate_all_data(self, algoData: AlgoDataType, schemes: List[str]) -> None:
        for algorithm, atomData in algoData.items():
            for atom, molData in atomData.items():
                for mol, basisData in molData.items():
                    for scheme in schemes:
                        result = extrapolate_molecule_given_scheme(basisData, scheme)
                        if result is None:
                            self.smallBasisException.setdefault(
                                algorithm, {}
                            ).setdefault(atom, {})[mol] = scheme
                        else:
                            self.algoToCBS.setdefault(algorithm, {}).setdefault(
                                atom, {}
                            ).setdefault(mol, {})[scheme] = result

    def calculate_errors(self, experimentalData: ExperDataType) -> None:
        algoToError: AlgoErrType = {}
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
                            if data["cbs"] is None:
                                raise ValueError(
                                    f"Expected to find cbs value for {algorithm}->{atom}->{mol}->{scheme} but found None instead"
                                )
                            error = data["cbs"] - experimentalData[mol]
                            algoToError.setdefault(algorithm, {}).setdefault(
                                atom, {}
                            ).setdefault(mol, {})[scheme] = error
        self.algoToError = algoToError

    def calculate_series_statistics(self) -> None:
        assert hasattr(self, "algoToError"), "You need to call calculate_errors first"
        algoToAtomErrors: AlgoAtomErrsType = {}
        for algorithm, atomData in self.algoToError.items():
            for atom, molData in atomData.items():
                for schemeData in molData.values():
                    for scheme, error in schemeData.items():
                        algoToAtomErrors.setdefault(algorithm, {}).setdefault(
                            atom, {}
                        ).setdefault(scheme, []).append(error)

        algoToAtomStats: AlgoAtomStatsType = {}
        for algorithm, atomErrData in algoToAtomErrors.items():
            for atom, schemeErrData in atomErrData.items():
                for scheme, errors_list in schemeErrData.items():
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
                        "errors": errors,
                        "abs_errors": abs_errs,
                    }
        self.algoToAtomStats = algoToAtomStats

    def calculate_overall_statistics(self) -> None:
        assert hasattr(self, "algoToError"), "You need to call calculate_errors first"
        algoToErrors: AlgoErrsType = {}
        for algorithm, atomData in self.algoToError.items():
            for molData in atomData.values():
                for schemeData in molData.values():
                    for scheme, error in schemeData.items():
                        algoToErrors.setdefault(algorithm, {}).setdefault(
                            scheme, []
                        ).append(error)

        algoToStats: AlgoStatsType = {}
        for algorithm, schemeErrsData in algoToErrors.items():
            for scheme, errors_list in schemeErrsData.items():
                errors = np.array(errors_list)
                abs_errs = np.abs(errors)
                algoToStats.setdefault(algorithm, {})[scheme] = {
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

    def _create_extrapolation_schemes(self) -> None:
        # CCSD schemes
        ccsd_schemes = []
        for bases in "D-T T-Q D-T-Q".split():
            for method in "CCSD CCSD(T)".split():
                ccsd_schemes.append(f"{bases}-{method}")

        bases_combs = "pcX-1 pcX-2 | pcX-2 pcX-3 | pcX-1 pcX-2 pcX-3 | ccX-DZ ccX-TZ | ccX-TZ ccX-QZ | ccX-DZ ccX-TZ ccX-QZ".split(
            " | "
        )
        for bases in bases_combs:
            for method in "CCSD CCSD(T)".split():
                ccsd_schemes.append(f"{method}[{bases}]")
        if self.debug:
            print(f"{ccsd_schemes=}")

        hf_schemes = []
        for bases in "T-Q D-T-Q T-Q-5 D-T-Q-5".split():
            hf_schemes.append(f"{bases}-HF")
        if self.debug:
            print(f"{hf_schemes=}")

        mp2_schemes = []
        mp2_schemes_ext = []
        mp2_bases = "D | T | Q | 5 | D T | D T Q | T Q | Q 5 | T Q 5 | D T Q 5"
        for bases in mp2_bases.split(" | "):
            if len(bases) != 1:
                mp2_schemes.append(f"MP2[{bases}]")
            if bases != "D":
                mp2_schemes.append(f"MP2[{bases}]+DifD")
                mp2_schemes.append(f"MP2[{bases}]+DifD(T)")
            for (
                small_basis
            ) in "STO-3G STO-6G 3-21G 4-31G 6-31G def2svp def2svpd".split():
                mp2_schemes_ext.append(f"MP2[{bases}]+Dif{small_basis}")
                mp2_schemes_ext.append(f"MP2[{bases}]+Dif{small_basis}(T)")

        if self.debug:
            print(f"{mp2_schemes=}")

        self.schemeDict = {
            "HF": hf_schemes,
            "CCSD": ccsd_schemes,
            "MP2": mp2_schemes,
            "MP2_EXT": mp2_schemes_ext,
        }

        self._schemeIterKeys = ["HF", "CCSD", "MP2", "MP2_EXT"]

    def scheme_generator(self) -> Iterable[str]:
        assert hasattr(
            self, "schemeDict"
        ), "You need to call _create_extrapolation_schemes first"
        for key in self._schemeIterKeys:
            schemes = self.schemeDict[key]
            for scheme in schemes:
                yield scheme

    @property
    def schemes(self) -> Iterable[str]:
        print("call schemes")
        return self.scheme_generator()


if __name__ == "__main__":
    atomToNames = {
        "C": "D-T-Q-CCSD(T), MP2(D T Q)+DifD(T)".split(", "),
        "N": "D-T-Q-CCSD(T), MP2(D T Q 5)+DifD(T)".split(", "),
        "O": "T-Q-CCSD, MP2(Q)+DifD".split(", "),
        "F": "T-Q-CCSD(T), MP2(Q)+DifD".split(", "),
    }

    # obj = WholeDataset()
    o = parse_scheme("CCSD[pcX-1 pcX-2]")
    print(f"{o=}")
