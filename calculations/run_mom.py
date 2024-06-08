import pprint
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict

import configs
import numpy as np
from BasisSets import (
    pyfiles as nwchem,  # see basis_sets module to understand where it comes from
)
from dmp2 import sgscf
from pyscf import gto
from pyscf.tools import cubegen

pp = pprint.PrettyPrinter(indent=4)

basisDictType = Dict[str, str]


def separate_txt_by_element(path: str) -> Dict[str, str]:
    with open(path, "r") as file:
        lines = file.readlines()

    basis_sets: Dict[str, str] = {}
    withinBasis = False
    curElement = None
    for i, rline in enumerate(lines):
        line = rline.strip()
        if line.startswith("#BASIS SET"):
            withinBasis = True
            curElement = lines[i + 1].strip().split()[0]
            continue
        if line.startswith("END"):
            break
        if not withinBasis:
            continue
        if curElement not in basis_sets:
            assert isinstance(curElement, str), "Element not instantiated"
            basis_sets[curElement] = ""
        basis_sets[curElement] += rline
    return basis_sets


def find_basis_files() -> Dict[str, basisDictType]:
    basis = Path(__file__).resolve().parent / "BasisSets"
    # find all files that end with .txt
    basis_files = list(basis.glob("*.txt"))
    keyToBasis: Dict[str, basisDictType] = {}
    for file in basis_files:
        fname = file.stem
        atomToBasis = separate_txt_by_element(file)
        atomToBasis["H"] = "cc-pV"
        keyToBasis[fname] = atomToBasis
    return keyToBasis


def zeta_to_h_basis(zeta: int) -> str:
    return f"cc-pV{zeta}Z"


# find_basis_files()
class CEBECalculator:

    def __init__(self, PARAMS):
        # heavyatom, nswap, geom, core_basis, full_basis, doLoc, testMode
        self.HEAVYATOM = PARAMS["atom"]
        self.NSWAP = PARAMS["orbital"]
        self.GEOM = PARAMS["geom"]
        self.CORE_BASIS = PARAMS["corebasis"]
        self.FULL_BASIS = PARAMS["regularbasis"]
        self.DOLOC = PARAMS.get("localize", False)
        self.CONFIGKEY = PARAMS.get("config", self.FULL_BASIS)
        self.DIIS_START = PARAMS.get("diis_start", 0)
        self.DIIS_SPACE = PARAMS.get("diis_space", 50)
        self.MAX_SCF_CYCLES = PARAMS.get("scf_maxcycles", 1000)
        self.LEVEL_SHIFT = PARAMS.get("level_shift", 0)
        self.ATOMS_TO_LOCALIZE = PARAMS.get("atoms_to_localize", self.HEAVYATOM).split()

        self.externalBasisLibrary = {
            # 'cc-pVTZ': nwchem.deconccPVTZ.BASIS_SET,
            # '4-31G': nwchem.g431.BASIS_SET,
            "6-311+G(3df)": nwchem.g6311plus3df.BASIS_SET,
        } | find_basis_files()

        self.keyToConfig = {
            "cc-pVQZ": configs.setups.Q_CONFIG,
            "pc-3": configs.setups.Q_CONFIG,
            "cc-pV5Z": configs.setups.PENTUPLE_CONFIG,
            "pc-4": configs.setups.PENTUPLE_CONFIG,
            "test": configs.setups.TEST_CONFIG,
            "mp3": configs.setups.MP3_CONFIG,
            "local": configs.setups.LOCAL_CONFIG,
        }

        self.GEN_CONFIG = configs.setups.GEN_CONFIG

        self.config = self.keyToConfig.get(self.CONFIGKEY, self.GEN_CONFIG)

        self.doTime = configs.parameters.doTimestamps
        self.DO_T1_DIAGN = configs.parameters.DO_T1_DIAGN
        self.VERBOSITY = configs.parameters.VERBOSITY

        if self.doTime:
            self.timefile = open("timestamps.txt", "w")

        self.doMP2 = self.config.get("doMP2", self.GEN_CONFIG["doMP2"])
        self.doMP3 = self.config.get("doMP3", self.GEN_CONFIG["doMP3"])
        self.doCCSD = self.config.get("doCCSD", self.GEN_CONFIG["doCCSD"])
        self.doTriples = self.config.get("doTriples", self.GEN_CONFIG["doTriples"])
        self.doSFX = self.config.get("doSFX", self.GEN_CONFIG["doSFX"])
        self.toFreeze = self.config.get("toFreeze", self.GEN_CONFIG["toFreeze"])
        self.doChengBasis = self.config.get(
            "doChengBasis", self.GEN_CONFIG["doChengBasis"]
        )
        self.doSpecialBasis = PARAMS.get(
            "doSpecialBasis", self.GEN_CONFIG["doSpecialBasis"]
        )
        self.basisLibKey = self.config.get(
            "basisLibKey", self.GEN_CONFIG["basisLibKey"]
        )
        self.toPrintDensity = self.config.get(
            "toPrintDensity", self.GEN_CONFIG["toPrintDensity"]
        )

    def get_time_now(self):
        return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    def core_basis_dict(self, mol):
        """
            This function makes a dictionary giving core_basis to the heaviest
            element and full_basis to the other elements

        Args:
            mol (obj): a gto.M object
            heavyatomsymbol (str): element symbol
            core_basis (str): a basis for core atom
            full_basis (str): a basis for noncore atoms

        Returns:
            dict: mapping atom symbols to basis sets
        """
        # First make a set of all atom symbols.
        atom_symbols = set()
        for a in range(mol.natm):
            atom_symbols.add(mol.atom_symbol(a))

        # Make the basis dict
        basis_dict = {}
        heavy_atom = self.HEAVYATOM
        core_basis = self.CORE_BASIS
        full_basis = self.FULL_BASIS
        doChengBasis = self.doChengBasis

        if doChengBasis:
            zeta = full_basis.split("V")[1][0]
            if zeta == "T":
                basisDic = nwchem.chengT.BASIS_SET
            elif zeta == "D":
                basisDic = nwchem.chengD.BASIS_SET
            basis_dict[heavy_atom] = gto.basis.parse(basisDic[heavy_atom]["core"])
            for a in atom_symbols:
                if a != heavy_atom:
                    basis_dict[a] = gto.basis.parse(basisDic[a]["reg"])
        elif self.doSpecialBasis:
            basisDic = self.externalBasisLibrary[core_basis]
            for a in atom_symbols:
                basis_dict[a] = full_basis
            basis_dict[heavy_atom] = gto.basis.parse(basisDic[heavy_atom])
        else:
            basis_dict[heavy_atom] = core_basis
            for a in atom_symbols:
                if a != heavy_atom:
                    basis_dict[a] = full_basis
        return basis_dict

    def _print_welcome_message(self):
        print("\t*********************************************** ")
        print("\t***** Beginning CEBE Calculation with MOM ***** ")
        print("\t*****       Calculation Details:          ***** ")
        print("\t      Geometry:\t", self.GEOM)
        print("\t      Core Basis:\t", self.CORE_BASIS)
        print("\t      Full Basis:\t", self.FULL_BASIS)
        print("\t      DIIS (start, space):\t", self.DIIS_START, self.DIIS_SPACE)
        print("\t      Level shift:\t", self.LEVEL_SHIFT)
        print("\t      Do localization?:\t", self.DOLOC)
        if self.DOLOC:
            print("\t      Atoms that were localized?:\t", self.ATOMS_TO_LOCALIZE)
        print("\t      Do SFX2C1E?:\t", self.doSFX)
        print("\t      Freeze orbitals?:\t", self.toFreeze)
        print("\t      Use Cheng orbitals?:\t", self.doChengBasis)
        print("\t      Do Special Basis?:\t", self.doSpecialBasis)
        if self.doSpecialBasis:
            print("\t      Do Special Basis?:\t", self.CORE_BASIS)
        print("\t*********************************************** ")
        print("\t*********************************************** ")

    def _build_mol(self):
        mol = gto.M()
        mol.fromfile(self.GEOM)
        mol.basis = self.core_basis_dict(mol)
        mol.build()
        self.mol = mol

    def _do_rhf(self):
        # RHF for ground state
        mf_gs = sgscf.RHF(self.mol)
        if self.doSFX:
            mf_gs.sfx2c1e()
        mf_gs.kernel()

        if self.toPrintDensity:
            self.P_GS = mf_gs.make_rdm1()  # store density matrix

        self.C0 = mf_gs.mo_coeff  # Ground State MOs
        self.mf_gs = mf_gs

    def _do_mp2(self):
        self.mf_gs.do_mp2()

    def _do_mp3(self):
        self.mf_gs.do_mp3()

    def _do_ccsd(self, triples=False):
        self.mf_gs.do_ccsd(
            do_ccsd_t=triples, do_t1_diagn=self.DO_T1_DIAGN, verbose=self.VERBOSITY
        )

    # def _do_ccsd_t(self):
    #     self.mf_gs.do_ccsd(do_ccsd_t=True, do_t1_diagn=self.DO_T1_DIAGN, verbose=self.VERBOSITY)

    def _do_localization(self):
        mol = self.mol

        C_old = self.C0.copy()
        C_pyloc = self.C0.copy()
        from pyscf.lo.boys import Boys

        mo = Boys(mol, C_old[:, : mol.nelec[1]])
        mo.init_guess = None
        mo = mo.kernel()
        C_pyloc[:, : mol.nelec[1]] = mo
        self.C_new = C_pyloc

    def _prepare_mol_for_uhf(self):
        # UHF for core excited state
        self.mol.charge = 1
        self.mol.spin = -1  # Take away one alpha electron
        self.mol.build()
        # Prepare list of orbitals to freeze
        na, nb = self.mol.nelec
        orb_swap = [[[self.NSWAP, na]], []]
        self.orb_swap = orb_swap

        if self.toFreeze:
            self.freeze_orbs = [
                [
                    na,
                ],
                [],
            ]
        else:
            self.freeze_orbs = False

        umf_core = sgscf.UHF(self.mol)
        # umf_core = sgscf.SGUHF(self.mol)

        umf_core.method = "diis"
        umf_core.mom_start = 0  # Means do MOM
        umf_core.orb_swap = orb_swap  # Swaps core orbital with LUMO
        umf_core.max_cycle = self.MAX_SCF_CYCLES

        sh = self.LEVEL_SHIFT
        umf_core.level_shift_eta = [sh, sh]
        if self.doSFX:
            umf_core.sfx2c1e()
        self.umf_core = umf_core
        self.diis = sgscf.sgopt.DIIS(self.umf_core)
        self.diis.adiis.min_space = self.DIIS_START
        self.diis.adiis.space = self.DIIS_SPACE

    def _do_uhf(self):
        if self.DOLOC:
            self._do_localization()
            self.umf_core.kernel(
                mo_coeff0=np.asarray([self.C_new, self.C_new]), optimizer=self.diis
            )
            self.umf_core.analyze(self.C_new)
        else:
            self.umf_core.kernel(
                mo_coeff0=np.asarray(
                    [
                        self.C0,
                    ]
                    * 2
                ),
                optimizer=self.diis,
            )

        # umf_core.analyze(C0 = C0) # Analyze orbital excitation
        if self.toPrintDensity:
            P_ES = self.umf_core.make_rdm1()
            cubegen.density(self.mol, "pdiff.cube", self.P_GS - P_ES[0] - P_ES[1])

    def _do_ump2(self):
        self.umf_core.do_mp2(frozen=self.freeze_orbs)

    def _do_ump3(self):
        self.umf_core.do_mp3(frozen=self.freeze_orbs)

    def _do_uccsd(self, triples=False):
        self.umf_core.do_ccsd(
            do_ccsd_t=triples,
            do_t1_diagn=self.DO_T1_DIAGN,
            verbose=self.VERBOSITY,
            frozen=self.freeze_orbs,
        )

    def run(self):
        self._print_welcome_message()

        if self.doTime:
            self.timefile.write("- HF\n")
            self.timefile.write(f"{self.get_time_now()}\n")
        self._build_mol()
        # --------- RUN RHF --------
        self._do_rhf()
        # --------------------------
        if self.doTime:
            self.timefile.write(f"{self.get_time_now()}\n")

        if self.doMP2:
            if self.doTime:
                self.timefile.write("- MP2\n")
                self.timefile.write(f"{self.get_time_now()}\n")
            # --------- RUN MP2 --------
            self._do_mp2()
            # --------------------------
            if self.doTime:
                self.timefile.write(f"{self.get_time_now()}\n")

        if self.doMP3:
            if self.doTime:
                self.timefile.write("- MP3\n")
                self.timefile.write(f"{self.get_time_now()}\n")
            # --------- RUN MP3 --------
            self._do_mp3()
            # --------------------------
            if self.doTime:
                self.timefile.write(f"{self.get_time_now()}\n")

        if self.doCCSD:
            if self.doTime:
                self.timefile.write("- CCSD\n")
                self.timefile.write(f"{self.get_time_now()}\n")
            # --------- RUN CCSD --------
            self._do_ccsd(triples=self.doTriples)
            # --------------------------
            if self.doTime:
                self.timefile.write(f"{self.get_time_now()}\n")

        self._prepare_mol_for_uhf()
        if self.doTime:
            self.timefile.write("- UHF\n")
            self.timefile.write(f"{self.get_time_now()}\n")
        # --------- RUN UHF --------
        self._do_uhf()
        # --------------------------
        if self.doTime:
            self.timefile.write(f"{self.get_time_now()}\n")

        def get_de(x, y):
            return (y - x) * 27.211399

        cebe_uhf = get_de(self.mf_gs.e_tot, self.umf_core.e_tot)
        f_out = open("CEBE_mom.txt", "w+")
        f_out.write("Frozen orbitals: %s\n" % self.freeze_orbs)
        f_out.write("Swapped orbitals: %s\n" % self.orb_swap)
        f_out.write("CEBE (UHF)    = %.6f eV\n" % cebe_uhf)

        if self.doMP2:
            if self.doTime:
                self.timefile.write("- UMP2\n")
                self.timefile.write(f"{self.get_time_now()}\n")
            # --------- RUN MP2 --------
            self._do_ump2()
            cebe_mp2 = get_de(self.mf_gs.e_mp2, self.umf_core.e_mp2)
            f_out.write("CEBE (MP2)    = %.6f eV\n" % cebe_mp2)
            # --------------------------
            if self.doTime:
                self.timefile.write(f"{self.get_time_now()}\n")

        if self.doMP3:
            if self.doTime:
                self.timefile.write("- UMP3\n")
                self.timefile.write(f"{self.get_time_now()}\n")
            # --------- RUN MP3 --------
            self._do_ump3()
            cebe_mp3 = get_de(self.mf_gs.e_mp3, self.umf_core.e_mp3)
            f_out.write("CEBE (MP3)    = %.6f eV\n" % cebe_mp3)
            # --------------------------
            if self.doTime:
                self.timefile.write(f"{self.get_time_now()}\n")

        if self.doCCSD:
            if self.doTime:
                self.timefile.write("- UCCSD\n")
                self.timefile.write(f"{self.get_time_now()}\n")
            # --------- RUN CCSD --------
            self._do_uccsd(triples=self.doTriples)
            cebe_ccsd = get_de(self.mf_gs.e_ccsd, self.umf_core.e_ccsd)
            f_out.write("CEBE (CCSD)   = %.6f eV\n" % cebe_ccsd)
            f_out.write("T1 for RHF: %s\n" % self.mf_gs.t1_diagn)
            f_out.write("T1 for UHF(a): %s\n" % self.umf_core.t1_diagn_a)
            f_out.write("T1 for UHF(b): %s\n" % self.umf_core.t1_diagn_b)
            # --------------------------
            if self.doTime:
                self.timefile.write(f"{self.get_time_now()}\n")

        if self.doTriples:
            cebe_ccsd_t = get_de(self.mf_gs.e_ccsd_t, self.umf_core.e_ccsd_t)
            f_out.write("CEBE (CCSD(T))   = %.6f eV\n" % cebe_ccsd_t)

        # print("PySCF PyANALYZE")
        # self.mf_gs.pyanalyze(verbose = self.VERBOSITY)
        f_out.close()


if __name__ == "__main__":
    """
    Command line arguments are:
        geom: (string) Path to .xyz file
        core_basis (string): Basis for core atom
        full_basis (string): Basis for other atoms
    """
    strToBool = {"false": False, "true": True}
    parameters = sys.argv
    assert parameters[0].endswith(
        ".py"
    ), "The first parameter should be path to the .py launch file"
    paramDic = {}
    case_sensitive = {
        "corebasis",
        "regularbasis",
        "atom",
        "geom",
        "atoms_to_localize",
        "dospecialbasis",
    }
    for i, param in enumerate(parameters):
        if i == 0:
            continue
        key, val = param.split("=")
        if key.lower() in case_sensitive:
            paramDic[key] = val
        else:
            paramDic[key.lower()] = val.lower()
    assert "geom" in paramDic, "You must provide path to geometry file (use geom=/)"
    assert (
        "atom" in paramDic
    ), "You must provide symbol of an element to which core basis will be assigned to (use atom=)"
    assert (
        "orbital" in paramDic
    ), "You must provide index (0-based indexing) of an orbital to be excited (use orbital=)"
    assert (
        "regularbasis" in paramDic
    ), "You must specify the basis set (use regularbasis=)"
    assert (
        "corebasis" in paramDic
    ), "You must specify the basis set for atom that will be excited (use corebasis=)"
    if "localize" in paramDic:
        paramDic["localize"] = strToBool[paramDic["localize"]]
    if "doSpecialBasis" in paramDic:
        paramDic["doSpecialBasis"] = strToBool[paramDic["doSpecialBasis"]]

    toInt = "orbital scf_maxcycles diis_start diis_space".split()
    for key in toInt:
        if key in paramDic:
            paramDic[key] = int(paramDic[key])
    toInt = "level_shift".split()
    for key in toInt:
        if key in paramDic:
            paramDic[key] = float(paramDic[key])
    predictor = CEBECalculator(paramDic)
    # print(paramDic)
    predictor.run()
    # CebeCalc = CEBECalculator({'atom': "O", "orbital": 0, "geom": "geom.xyz", "corebasis": "cc-pVDZ", "regularbasis": "cc-pVDZ", })
