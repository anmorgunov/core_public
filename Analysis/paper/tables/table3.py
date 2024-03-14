import sys, Analysis.constants as constants

sys.path.append(
    "../../core_excitations"
)  # this allows us to import files from parent directories
import PRIVATE
import tools
import Analysis.constants as constants
import analysis.modules.table
import analysis.regression.evaluate
import analysis.extrapolation.extrapolations
from typing import Optional


class Tables:
    def __init__(self):
        self.latT = analysis.modules.table.LaTeXTable()
        self.evalObj = analysis.regression.evaluate.export()
        self.extrapObj = analysis.extrapolation.extrapolations.export()
        self.molToExper = self.evalObj.molToExper
        # self.atomToMols = constants.ATOM_TO_MOLS
        self.molToAtom = self.evalObj.molWithExpToAtom
        self.data = self.evalObj.parserObj.methodToData
        self.MOL = constants.FNAME_TO_MOLS

        self.atomToMols = {}
        for mol, atom in self.molToAtom.items():
            self.atomToMols.setdefault(atom, set()).add(mol)

    def series_table(
        self, header, body, positioning, atoms, suffix: Optional[str] = None
    ):
        if len(atoms) == 1:
            atom = atoms[0]
        else:
            atom = "all"
        caption = (
            f"Statistical analysis of accuracy of different extrapolation schemes at predicting K-Edge CEBEs (in eV) compared to experimental data for {atom}-series"+ suffix)
        label = f"summary-{atom.lower()}"
        path = tools.join_path(
            PRIVATE.BASE_PATH
            + [
                "analysis",
                "paper",
                "output",
                "tables",
                "extrapschemes",
                f"{atom}-summary",
            ]
        )
        self.latT.export_table(caption, label, positioning, header, body, path)

    def summary(self):
        atomCombos = (("C",), ("N",), ("O"), ("F"), ("C", "N", "O", "F"))
        for atoms in atomCombos:
            header, body = self.extrapObj.export_errors(self.extrapObj.all, atoms)
            header = [
                "\\textbf{%s}" % colname
                for i, colname in enumerate(header)
                if i != len(header) - 1
            ]
            n_samples = set([row[-1] for row in body])
            assert (
                len(n_samples) == 1
            ), "Number of samples must be the same for all methods"
            n_samples = n_samples.pop()
            suffix = f" ({n_samples} molecules)"
            new_body = []
            for i, row in enumerate(body):
                newRow = row[:-1]
                name = row[0]
                if "MP2" in name:
                    comps = name.split('+')
                    print(comps)
                    if len(comps) == 1:
                        big = comps[0]
                        corr = ""
                    else:
                        big, corr = comps
                        corr = "+" + corr
                    big = big.replace("(", "[").replace(")", "]")
                    name = f"{big}{corr}"
                    newRow[0] = name
                new_body.append(newRow)
            positioning = "l " * len(header)
            self.series_table(header, new_body, positioning, atoms, suffix)

    def separate_in_two(self):
        atomCombos = (("C",), ("N",), ("O"), ("F"), ("C", "N", "O", "F"))
        atomToBody = {}
        for atoms in atomCombos:
            header, body = self.extrapObj.export_errors(self.extrapObj.all, atoms)
            if len(atoms) != 1:
                continue
            atomToBody[atoms[0]] = body

        # ['Extrapolation Scheme', 'MSE', 'MAE', 'MedAE', 'MaxAE', 'STD', 'Sample Size']
        topheader = ["", "\multicolumn{3}{c}{C-Series}", "\multicolumn{3}{c}{N-Series}"]
        lowheader = [
            "Extrapolation Scheme",
            "MAE",
            "MaxAE",
            "STD",
            "MAE",
            "MaxAE",
            "STD",
        ]
        positioning = "l " * len(lowheader)
        newBody = []
        for i, row in enumerate(atomToBody["C"]):
            n, o, f = atomToBody["N"], atomToBody["O"], atomToBody["F"]
            newRow = row[:1] + row[2:3] + row[4:6] + n[i][2:3] + n[i][4:6]
            newBody.append(newRow)
        caption = f"Errors of K-Edge Ionization Energies (in eV) in the CBS limit compared to experimental data"
        label = f"cn-series-wide"
        path = tools.join_path(
            PRIVATE.BASE_PATH
            + ["analysis", "paper", "output", "tables", "extrapschemes", f"cn-summary"]
        )
        header = [topheader, lowheader]
        self.latT.export_table(caption, label, positioning, header, newBody, path)

        topheader = ["", "\multicolumn{3}{c}{O-Series}", "\multicolumn{3}{c}{F-Series}"]
        newBody = []
        for i, row in enumerate(atomToBody["C"]):
            n, o, f = atomToBody["N"], atomToBody["O"], atomToBody["F"]
            newRow = row[:1] + o[i][2:3] + o[i][4:6] + f[i][2:3] + f[i][4:6]
            newBody.append(newRow)
        caption = f"Errors of K-Edge Ionization Energies (in eV) in the CBS limit compared to experimental data"
        label = f"of-series-wide"
        path = tools.join_path(
            PRIVATE.BASE_PATH
            + ["analysis", "paper", "output", "tables", "extrapschemes", f"of-summary"]
        )
        header = [topheader, lowheader]
        self.latT.export_table(caption, label, positioning, header, newBody, path)

    def convergence_issues(self):
        header, body = self.extrapObj.export_small_basis_convergence_issues()
        caption = f"Number of molecule with non-converging UHF in small basis sets"
        label = f"conv-issues"
        positioning = "l " + "c " * (len(header) - 1)
        path = tools.join_path(
            PRIVATE.BASE_PATH
            + ["analysis", "paper", "output", "tables", "extrapschemes", f"conv-issues"]
        )
        self.latT.export_table(caption, label, positioning, header, body, path)


if __name__ == "__main__":
    tables = Tables()
    tables.summary()
    # tables.separate_in_two()
    # tables.convergence_issues()
