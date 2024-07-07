from typing import Dict, List

# fmt:off
BASIS_TO_METHOD = {
    "D": "UHF MP2 CCSD CCSD(T)".split(),
    "T": "UHF MP2 CCSD CCSD(T)".split(),
    "Q": "UHF MP2 CCSD CCSD(T)".split(),
    "5": "UHF MP2".split(),
    "pcX-1": "UHF MP2 CCSD CCSD(T)".split(),
    "pcX-2": "UHF MP2 CCSD CCSD(T)".split(),
    "pcX-3": "UHF MP2 CCSD CCSD(T)".split(),
    "pcX-4": "UHF MP2".split(),
    "ccX-DZ": "UHF MP2 CCSD CCSD(T)".split(),
    "ccX-TZ": "UHF MP2 CCSD CCSD(T)".split(),
    "ccX-QZ": "UHF MP2 CCSD CCSD(T)".split(),
    "ccX-5Z": "UHF MP2".split(),
}

METHOD_TO_BASIS: Dict[str, List[str]] = {}
for basis, methods in BASIS_TO_METHOD.items():
    for method in methods:
        METHOD_TO_BASIS.setdefault(method, []).append(basis)
errBarrSummary_w5 = {
    "c": ["MP2[D T Q 5]", "MP2[D T Q 5]+DifSTO-3G(T)", "MP2[D T Q 5]+Dif3-21G(T)", "MP2[D T Q 5]+DifD(T)", "D-T-Q-CCSD(T)"],
    "n": ["MP2[T Q 5]", "MP2[T Q 5]+DifSTO-3G", "MP2[T Q 5]+Dif3-21G", "MP2[T Q 5]+DifD", "T-Q-CCSD"],
    "o": ["MP2[T Q 5]", "MP2[T Q 5]+DifSTO-3G", "MP2[T Q 5]+Dif3-21G", "MP2[T Q 5]+DifD", "T-Q-CCSD"],
    "f": ["MP2[T Q 5]", "MP2[T Q 5]+DifSTO-3G(T)", "MP2[T Q 5]+Dif3-21G(T)", "MP2[T Q 5]+DifD(T)", "T-Q-CCSD(T)"],
}

errBarSummary_no5 = {
    "c": ["MP2[D T Q]", "MP2[D T Q]+DifSTO-3G(T)", "MP2[D T Q]+Dif3-21G(T)", "MP2[D T Q]+DifD(T)", "D-T-Q-CCSD(T)"],
    "n": ["MP2[T Q]", "MP2[T Q]+DifSTO-3G", "MP2[T Q]+Dif3-21G", "MP2[T Q]+DifD", "T-Q-CCSD"],
    "o": ["MP2[T Q]", "MP2[T Q]+DifSTO-3G", "MP2[T Q]+Dif3-21G", "MP2[T Q]+DifD", "T-Q-CCSD"],
    "f": ["MP2[T Q]", "MP2[T Q]+DifSTO-3G(T)", "MP2[T Q]+Dif3-21G(T)", "MP2[T Q]+DifD(T)", "T-Q-CCSD(T)"],
}

extrapOrNot_w5 = {
    "c": ["MP2[D T Q 5]+DifD(T)", "D-T-Q-CCSD(T)", "CCSD(T)", "CCSD(T)", "CCSD(T)"],
    "n": ["MP2[T Q 5]+DifD", "T-Q-CCSD", "CCSD", "CCSD", "CCSD"],
    "o": ["MP2[T Q 5]+DifD", "T-Q-CCSD", "CCSD", "CCSD", "CCSD"],
    "f": ["MP2[T Q 5]+DifD(T)", "T-Q-CCSD(T)", "CCSD", "CCSD", "CCSD"],
}

extrapOrNot_no5 = {
    "c": ["MP2[D T Q]+DifD(T)", "D-T-Q-CCSD(T)", "CCSD(T)", "CCSD(T)", "CCSD(T)"],
    "n": ["MP2[T Q]+DifD", "T-Q-CCSD", "CCSD", "CCSD", "CCSD"],
    "o": ["MP2[T Q]+DifD", "T-Q-CCSD", "CCSD", "CCSD", "CCSD"],
    "f": ["MP2[T Q]+DifD(T)", "T-Q-CCSD(T)", "CCSD", "CCSD", "CCSD"],
}

smallBasisStudy = {
    "c": ["MP2[D T Q]+DifSTO-3G(T)", "MP2[D T Q]+DifSTO-6G(T)", "MP2[D T Q]+Dif3-21G(T)", "MP2[D T Q]+Dif4-31G(T)", "MP2[D T Q]+Dif6-31G(T)", "MP2[D T Q]+DifD(T)"],
    "n": ["MP2[T Q]+DifSTO-3G", "MP2[T Q]+DifSTO-6G", "MP2[T Q]+Dif3-21G", "MP2[T Q]+Dif4-31G", "MP2[T Q]+Dif6-31G", "MP2[T Q]+DifD"],
    "o": ["MP2[T Q]+DifSTO-3G", "MP2[T Q]+DifSTO-6G", "MP2[T Q]+Dif3-21G", "MP2[T Q]+Dif4-31G", "MP2[T Q]+Dif6-31G", "MP2[T Q]+DifD"],
    "f": ["MP2[T Q]+DifSTO-3G(T)", "MP2[T Q]+DifSTO-6G(T)", "MP2[T Q]+Dif3-21G(T)", "MP2[T Q]+Dif4-31G(T)", "MP2[T Q]+Dif6-31G(T)", "MP2[T Q]+DifD(T)"],
}

bigBasisStudy = {
    "c": ["MP2[D T Q]+DifD(T)", "MP2[D T Q 5]+DifD(T)", "MP2[T Q]+DifD(T)", "MP2[T Q 5]+DifD(T)", "MP2[Q 5]+DifD(T)", "MP2[5]+DifD(T)"],
    "n": ["MP2[D T Q]+DifD", "MP2[D T Q 5]+DifD", "MP2[T Q]+DifD", "MP2[T Q 5]+DifD", "MP2[Q 5]+DifD", "MP2[5]+DifD"],
    "o": ["MP2[D T Q]+DifD", "MP2[D T Q 5]+DifD", "MP2[T Q]+DifD", "MP2[T Q 5]+DifD", "MP2[Q 5]+DifD", "MP2[5]+DifD",],
    "f": ["MP2[D T Q]+DifD(T)", "MP2[D T Q 5]+DifD(T)", "MP2[T Q]+DifD(T)", "MP2[T Q 5]+DifD(T)", "MP2[Q 5]+DifD(T)", "MP2[5]+DifD(T)"],
}

corrSummary = {
    "do5": {
        "c": [["MP2[D T Q 5]+DifD", "D-T-Q-CCSD"], ["MP2[D T Q 5]+Dif3-21G", "D-T-Q-CCSD"], ["MP2[D T Q 5]+DifD(T)", "D-T-Q-CCSD(T)"], ["MP2[D T Q 5]+Dif3-21G(T)", "D-T-Q-CCSD(T)"]],
        "n": [["MP2[T Q 5]+DifD", "T-Q-CCSD"], ["MP2[T Q 5]+Dif3-21G", "T-Q-CCSD"], ["MP2[T Q 5]+DifD(T)", "T-Q-CCSD(T)"], ["MP2[T Q 5]+Dif3-21G(T)", "T-Q-CCSD(T)"]],
        "o": [["MP2[T Q 5]+DifD", "T-Q-CCSD"], ["MP2[T Q 5]+Dif3-21G", "T-Q-CCSD"], ["MP2[T Q 5]+DifD(T)", "T-Q-CCSD(T)"], ["MP2[T Q 5]+Dif3-21G(T)", "T-Q-CCSD(T)"]],
        "f": [["MP2[T Q 5]+DifD", "T-Q-CCSD"], ["MP2[T Q 5]+Dif3-21G", "T-Q-CCSD"], ["MP2[T Q 5]+DifD(T)", "T-Q-CCSD(T)"], ["MP2[T Q 5]+Dif3-21G(T)", "T-Q-CCSD(T)"]],
    },
    "no5": {
        "c": [["MP2[D T Q]+DifD", "D-T-Q-CCSD"], ["MP2[D T Q]+Dif3-21G", "D-T-Q-CCSD"], ["MP2[D T Q]+DifD(T)", "D-T-Q-CCSD(T)"], ["MP2[D T Q]+Dif3-21G(T)", "D-T-Q-CCSD(T)"]],
        "n": [["MP2[D T Q]+DifD", "D-T-Q-CCSD"],["MP2[D T Q]+Dif3-21G", "D-T-Q-CCSD"],["MP2[D T Q]+DifD(T)", "D-T-Q-CCSD(T)"],["MP2[D T Q]+Dif3-21G(T)", "D-T-Q-CCSD(T)"]],
        "o": [["MP2[T Q]+DifD", "T-Q-CCSD"], ["MP2[T Q]+Dif3-21G", "T-Q-CCSD"], ["MP2[T Q]+DifD(T)", "T-Q-CCSD(T)"], ["MP2[T Q]+Dif3-21G(T)", "T-Q-CCSD(T)"]],
        "f": [["MP2[T Q]+DifD", "T-Q-CCSD"], ["MP2[T Q]+Dif3-21G", "T-Q-CCSD"], ["MP2[T Q]+DifD(T)", "T-Q-CCSD(T)"], ["MP2[T Q]+Dif3-21G(T)", "T-Q-CCSD(T)"]],
    },
}

# molecules with convergence issues in 3-21G basis
corrSmBasisException = {"hccf", "me2-n-cho", "h2n-cho", "ch3-c-ho", "hcooch3", "mnh2pyridi-n-e", "hconh2", "cf3co-o-h", "c-h3cch", "hno2", "ch3-c-n", "onh2pyridi-n-e", "hcooh", "c-h2o", "c-o", "ch3cof", "hco-n-hch3"}
