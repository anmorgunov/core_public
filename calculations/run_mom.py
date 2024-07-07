import argparse
import pprint
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Set, Union

import configs
import numpy as np
from BasisSets import pyfiles as nwchem
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
        atomToBasis = separate_txt_by_element(str(file))
        atomToBasis["H"] = "cc-pV"
        keyToBasis[fname] = atomToBasis
    return keyToBasis


def zeta_to_h_basis(zeta: int) -> str:
    return f"cc-pV{zeta}Z"


# find_basis_files()
class CEBECalculator:
    keyToConfig: Dict[str, configs.setups.ConfigType] = {
        "quadruple": configs.setups.Q_CONFIG,
        "pentuple": configs.setups.PENTUPLE_CONFIG,
        "test": configs.setups.TEST_CONFIG,
        "mp3": configs.setups.MP3_CONFIG,
        "local": configs.setups.LOCAL_CONFIG,
    }

    def __init__(
        self,
        element_to_ionize: str,
        orbital_to_ionize: int,
        path_to_geometry: str,
        core_basis: str,
        regular_basis: str,
        do_special_basis: bool = False,
        localize: bool = False,
        config_key: str = "quadruple",
        diis_start: int = 0,
        diis_space: int = 50,
        scf_maxcycles: int = 1000,
        level_shift: float = 0,
        atoms_to_localize: str = "",
    ) -> None:
        # heavyatom, nswap, geom, core_basis, full_basis, doLoc, testMode
        self.element_to_ionize = element_to_ionize
        self.n_swap = orbital_to_ionize
        self.path_to_geometry = path_to_geometry
        self.core_basis = core_basis
        self.regular_basis = regular_basis
        self.do_special_basis = do_special_basis
        self.localize = localize
        self.config_key = config_key
        self.diis_start = diis_start
        self.diis_space = diis_space
        self.scf_maxcycles = scf_maxcycles
        self.level_shift = level_shift
        self.atoms_to_localize = atoms_to_localize

        self.externalBasisLibrary: Dict[str, basisDictType] = {
            "6-311+G(3df)": nwchem.g6311plus3df.BASIS_SET,
        } | find_basis_files()

        self.config: configs.setups.ConfigType = self.keyToConfig.get(
            self.config_key, configs.setups.GEN_CONFIG
        )

        self.save_time = configs.parameters.doTimestamps
        self.do_t1_diagnostics = configs.parameters.DO_T1_DIAGN
        self.verbosity = configs.parameters.VERBOSITY

        if self.save_time:
            self.timefile = open("timestamps.txt", "w")

    def get_time_now(self) -> str:
        return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    def core_basis_dict(self, element_symbols: Set[str]) -> basisDictType:
        """
        This function makes a dictionary giving core_basis to the heaviest
        element and full_basis to the other elements
        """

        basis_dict: basisDictType = {}

        if self.config["doChengBasis"]:
            zeta = self.regular_basis.split("V")[1][0]
            if zeta == "T":
                basisDic = nwchem.chengT.BASIS_SET
            elif zeta == "D":
                basisDic = nwchem.chengD.BASIS_SET
            basis_dict[self.element_to_ionize] = gto.basis.parse(
                basisDic[self.element_to_ionize]["core"]
            )
            for element in element_symbols:
                if element != self.element_to_ionize:
                    basis_dict[element] = gto.basis.parse(basisDic[element]["reg"])
        elif self.do_special_basis:
            specBasisDic = self.externalBasisLibrary[self.core_basis]
            for element in element_symbols:
                basis_dict[element] = self.regular_basis
            basis_dict[self.element_to_ionize] = gto.basis.parse(
                specBasisDic[self.element_to_ionize]
            )
        else:
            basis_dict[self.element_to_ionize] = self.core_basis
            for element in element_symbols:
                if element != self.element_to_ionize:
                    basis_dict[element] = self.regular_basis
        return basis_dict

    def _print_welcome_message(self) -> None:
        print("\t*********************************************** ")
        print("\t***** Beginning CEBE Calculation with MOM ***** ")
        print("\t*****       Calculation Details:          ***** ")
        print("\t      Geometry:\t", self.path_to_geometry)
        print("\t      Core Basis:\t", self.core_basis)
        print("\t      Full Basis:\t", self.regular_basis)
        print("\t      DIIS (start, space):\t", self.diis_start, self.diis_space)
        print("\t      Level shift:\t", self.level_shift)
        print("\t      Do localization?:\t", self.localize)
        if self.localize:
            print("\t      Atoms that were localized?:\t", self.atoms_to_localize)
        print("\t      Do SFX2C1E?:\t", self.config["doSFX"])
        print("\t      Freeze orbitals?:\t", self.config["toFreeze"])
        print("\t      Use Cheng orbitals?:\t", self.config["doChengBasis"])
        print("\t      Do Special Basis?:\t", self.do_special_basis)
        if self.do_special_basis:
            print("\t      Special Basis:\t", self.core_basis)
        print("\t*********************************************** ")
        print("\t*********************************************** ")

    def _build_mol(self) -> None:
        mol = gto.M()
        mol.fromfile(self.path_to_geometry)
        atom_symbols = set()
        for a in range(mol.natm):
            atom_symbols.add(mol.atom_symbol(a))
        mol.basis = self.core_basis_dict(atom_symbols)
        mol.build()
        self.mol = mol

    def _do_rhf(self) -> None:
        # RHF for ground state
        mf_gs = sgscf.RHF(self.mol)
        if self.config["doSFX"]:
            mf_gs.sfx2c1e()
        mf_gs.kernel()

        if self.config["toPrintDensity"]:
            self.P_GS = mf_gs.make_rdm1()  # store density matrix

        self.C0 = mf_gs.mo_coeff  # Ground State MOs
        self.mf_gs = mf_gs

    def _do_mp2(self) -> None:
        self.mf_gs.do_mp2()

    def _do_mp3(self) -> None:
        self.mf_gs.do_mp3()

    def _do_ccsd(self, triples: bool = False) -> None:
        self.mf_gs.do_ccsd(
            do_ccsd_t=triples,
            do_t1_diagn=self.do_t1_diagnostics,
            verbose=self.verbosity,
        )

    def _do_localization(self) -> None:
        mol = self.mol

        C_old = self.C0.copy()
        C_pyloc = self.C0.copy()
        from pyscf.lo.boys import Boys

        mo = Boys(mol, C_old[:, : mol.nelec[1]])
        mo.init_guess = None
        mo = mo.kernel()
        C_pyloc[:, : mol.nelec[1]] = mo
        self.C_new = C_pyloc

    def _prepare_mol_for_uhf(self) -> None:
        # UHF for core excited state
        self.mol.charge = 1
        self.mol.spin = -1  # Take away one alpha electron
        self.mol.build()
        # Prepare list of orbitals to freeze
        na, nb = self.mol.nelec
        orb_swap = [[[self.n_swap, na]], []]
        self.orb_swap = orb_swap

        if self.config["toFreeze"]:
            self.freeze_orbs: Union[List[List[int]], bool] = [[na], []]
        else:
            self.freeze_orbs = False

        umf_core = sgscf.UHF(self.mol)
        # umf_core = sgscf.SGUHF(self.mol)

        umf_core.method = "diis"
        umf_core.mom_start = 0  # Means do MOM
        umf_core.orb_swap = orb_swap  # Swaps core orbital with LUMO
        umf_core.max_cycle = self.scf_maxcycles

        sh = self.level_shift
        umf_core.level_shift_eta = [sh, sh]
        if self.config["doSFX"]:
            umf_core.sfx2c1e()
        self.umf_core = umf_core
        self.diis = sgscf.sgopt.DIIS(self.umf_core)
        self.diis.adiis.min_space = self.diis_start
        self.diis.adiis.space = self.diis_space

    def _do_uhf(self) -> None:
        if self.localize:
            self._do_localization()
            self.umf_core.kernel(
                mo_coeff0=np.asarray([self.C_new, self.C_new]), optimizer=self.diis
            )
            self.umf_core.analyze(self.C_new)
        else:
            self.umf_core.kernel(
                mo_coeff0=np.asarray([self.C0] * 2),
                optimizer=self.diis,
            )

        # umf_core.analyze(C0 = C0) # Analyze orbital excitation
        if self.config["toPrintDensity"]:
            P_ES = self.umf_core.make_rdm1()
            cubegen.density(self.mol, "pdiff.cube", self.P_GS - P_ES[0] - P_ES[1])

    def _do_ump2(self) -> None:
        self.umf_core.do_mp2(frozen=self.freeze_orbs)

    def _do_ump3(self) -> None:
        self.umf_core.do_mp3(frozen=self.freeze_orbs)

    def _do_uccsd(self, triples: bool = False) -> None:
        self.umf_core.do_ccsd(
            do_ccsd_t=triples,
            do_t1_diagn=self.do_t1_diagnostics,
            verbose=self.verbosity,
            frozen=self.freeze_orbs,
        )

    def run(self) -> None:
        self._print_welcome_message()

        if self.save_time:
            self.timefile.write("- HF\n")
            self.timefile.write(f"{self.get_time_now()}\n")
        self._build_mol()
        # --------- RUN RHF --------
        self._do_rhf()
        # --------------------------
        if self.save_time:
            self.timefile.write(f"{self.get_time_now()}\n")

        if self.config["doMP2"]:
            if self.save_time:
                self.timefile.write("- MP2\n")
                self.timefile.write(f"{self.get_time_now()}\n")
            # --------- RUN MP2 --------
            self._do_mp2()
            # --------------------------
            if self.save_time:
                self.timefile.write(f"{self.get_time_now()}\n")

        if self.config["doMP3"]:
            if self.save_time:
                self.timefile.write("- MP3\n")
                self.timefile.write(f"{self.get_time_now()}\n")
            # --------- RUN MP3 --------
            self._do_mp3()
            # --------------------------
            if self.save_time:
                self.timefile.write(f"{self.get_time_now()}\n")

        if self.config["doCCSD"]:
            if self.save_time:
                self.timefile.write("- CCSD\n")
                self.timefile.write(f"{self.get_time_now()}\n")
            # --------- RUN CCSD --------
            self._do_ccsd(triples=self.config["doTriples"])
            # --------------------------
            if self.save_time:
                self.timefile.write(f"{self.get_time_now()}\n")

        self._prepare_mol_for_uhf()
        if self.save_time:
            self.timefile.write("- UHF\n")
            self.timefile.write(f"{self.get_time_now()}\n")
        # --------- RUN UHF --------
        self._do_uhf()
        # --------------------------
        if self.save_time:
            self.timefile.write(f"{self.get_time_now()}\n")

        def get_de(x: float, y: float) -> float:
            return (y - x) * 27.211399

        cebe_uhf = get_de(self.mf_gs.e_tot, self.umf_core.e_tot)
        f_out = open("CEBE_mom.txt", "w+")
        f_out.write("Frozen orbitals: %s\n" % self.freeze_orbs)
        f_out.write("Swapped orbitals: %s\n" % self.orb_swap)
        f_out.write("CEBE (UHF)    = %.6f eV\n" % cebe_uhf)

        if self.config["doMP2"]:
            if self.save_time:
                self.timefile.write("- UMP2\n")
                self.timefile.write(f"{self.get_time_now()}\n")
            # --------- RUN MP2 --------
            self._do_ump2()
            cebe_mp2 = get_de(self.mf_gs.e_mp2, self.umf_core.e_mp2)
            f_out.write("CEBE (MP2)    = %.6f eV\n" % cebe_mp2)
            # --------------------------
            if self.save_time:
                self.timefile.write(f"{self.get_time_now()}\n")

        if self.config["doMP3"]:
            if self.save_time:
                self.timefile.write("- UMP3\n")
                self.timefile.write(f"{self.get_time_now()}\n")
            # --------- RUN MP3 --------
            self._do_ump3()
            cebe_mp3 = get_de(self.mf_gs.e_mp3, self.umf_core.e_mp3)
            f_out.write("CEBE (MP3)    = %.6f eV\n" % cebe_mp3)
            # --------------------------
            if self.save_time:
                self.timefile.write(f"{self.get_time_now()}\n")

        if self.config["doCCSD"]:
            if self.save_time:
                self.timefile.write("- UCCSD\n")
                self.timefile.write(f"{self.get_time_now()}\n")
            # --------- RUN CCSD --------
            self._do_uccsd(triples=self.config["doTriples"])
            cebe_ccsd = get_de(self.mf_gs.e_ccsd, self.umf_core.e_ccsd)
            f_out.write("CEBE (CCSD)   = %.6f eV\n" % cebe_ccsd)
            f_out.write("T1 for RHF: %s\n" % self.mf_gs.t1_diagn)
            f_out.write("T1 for UHF(a): %s\n" % self.umf_core.t1_diagn_a)
            f_out.write("T1 for UHF(b): %s\n" % self.umf_core.t1_diagn_b)
            # --------------------------
            if self.save_time:
                self.timefile.write(f"{self.get_time_now()}\n")

        if self.config["doTriples"]:
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
    parser = argparse.ArgumentParser(description="CEBE Calculator")
    parser.add_argument("--geom", type=str, required=True, help="Path to .xyz file")
    parser.add_argument(
        "--element_to_ionize",
        type=str,
        required=True,
        help="Symbol of the element to which core basis will be assigned",
    )
    parser.add_argument(
        "--orbital",
        type=int,
        required=True,
        help="Index (0-based) of the orbital to be ionized",
    )
    parser.add_argument(
        "--regularbasis",
        type=str,
        required=True,
        help="Basis set for other atoms",
    )
    parser.add_argument(
        "--corebasis",
        type=str,
        required=True,
        help="Basis set for the atom that will be excited",
    )
    parser.add_argument("--localize", action="store_true", help="Enable localization")
    parser.add_argument(
        "--dospecialbasis", action="store_true", help="Enable special basis"
    )
    parser.add_argument(
        "--config", type=str, default="quadruple", help="Configuration key"
    )
    parser.add_argument("--diis_start", type=int, default=0, help="DIIS start")
    parser.add_argument("--diis_space", type=int, default=50, help="DIIS space")
    parser.add_argument(
        "--scf_maxcycles", type=int, default=1000, help="Maximum number of SCF cycles"
    )
    parser.add_argument("--level_shift", type=float, default=0.0, help="Level shift")
    parser.add_argument(
        "--atoms_to_localize", type=str, default=None, help="Atoms to localize"
    )

    args = parser.parse_args()

    predictor = CEBECalculator(
        element_to_ionize=args.element_to_ionize,
        orbital_to_ionize=args.orbital,
        path_to_geometry=args.geom,
        core_basis=args.corebasis,
        regular_basis=args.regularbasis,
        do_special_basis=args.dospecialbasis,
        localize=args.localize,
        config_key=args.config,
        diis_start=args.diis_start,
        diis_space=args.diis_space,
        scf_maxcycles=args.scf_maxcycles,
        level_shift=args.level_shift,
        atoms_to_localize=args.atoms_to_localize,
    )
    # print(paramDic)
    predictor.run()
    # CebeCalc = CEBECalculator({'atom': "O", "orbital": 0, "geom": "geom.xyz", "corebasis": "cc-pVDZ", "regularbasis": "cc-pVDZ", })
