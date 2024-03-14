import sys, Analysis.constants as constants
sys.path.append('../../core_excitations') # this allows us to import files from parent directories
import PRIVATE
import tools
import Analysis.constants as constants
import analysis.modules.table 
import analysis.regression.evaluate
from typing import Optional

class Tables:

    def __init__(self):
        self.latT = analysis.modules.table.LaTeXTable()
        self.evalObj = analysis.regression.evaluate.export()
        self.molToExper = self.evalObj.molToExper
        # self.atomToMols = constants.ATOM_TO_MOLS
        self.molToAtom = self.evalObj.molWithExpToAtom
        self.data = self.evalObj.parserObj.methodToData
        self.MOL = constants.FNAME_TO_MOLS

        self.atomToMols = {}
        for mol, atom in self.molToAtom.items():
            self.atomToMols.setdefault(atom, set()).add(mol)
    
    def series_table(self, header, body, positioning, atoms, suffix:Optional[str]=None):
        if len(atoms) == 1: atom = atoms.pop()
        else: atom = 'all'
        caption = f'Statistical analysis of accuracy of different methods at predicting K-Edge CEBEs (in eV) compared to experimental data for {atom}-series' + suffix if suffix else ''
        label = f'summary-{atom.lower()}'
        path = tools.join_path(PRIVATE.BASE_PATH+['analysis', 'paper', 'output', 'tables', 'methodsummary', f'{atom}-summary'])
        self.latT.export_table(caption, label, positioning, header, body, path)

    def summary(self):
        atomCombos = ({'C',}, {'N',}, {'O'}, {'F'}, {'C', 'N', 'O', 'F'})
        header, atomsBodies = self.evalObj.summary_of_deviations(atomCombos)
        header = ["\\textbf{%s}" % colname for i, colname in enumerate(header) if i != len(header)-1]
        positioning = 'l '*len(header)
        for atoms, body in atomsBodies:
            n_samples = set([row[-1] for row in body])
            assert len(n_samples) == 1, 'Number of samples must be the same for all methods'
            n_samples = n_samples.pop()
            suffix = f' ({n_samples} molecules)'
            body = [row[:-1] for row in body]
            self.series_table(header, body, positioning, atoms, suffix)
    
    def combined_table(self):
        atomCombos = (('C',), ('N',), ('O'), ('F'), ('C', 'N', 'O', 'F'))
        header, atomsBodies = self.evalObj.summary_of_deviations(atomCombos)
        atomToBody = {}
        for atoms, body in atomsBodies:
            if len(atoms) != 1: continue
            atomToBody[atoms[0]] = body

        topheader = ['Basis', 'Method', '\multicolumn{2}{c}{C-Series}', '\multicolumn{2}{c}{N-Series}', '\multicolumn{2}{c}{O-Series}', '\multicolumn{2}{c}{F-series}']
        lowheader = ['', '', 'MAE', 'STD', 'MAE', 'STD', 'MAE', 'STD', 'MAE', 'STD']
        newBody = []
        for i, row in enumerate(atomToBody['C']):
            n, o, f = atomToBody['N'], atomToBody['O'], atomToBody['F']
            newRow = row[:2] + row[3:4] + row[6:7] + n[i][3:4] + n[i][6:7] + o[i][3:4] + o[i][6:7] + f[i][3:4] + f[i][6:7]
            newBody.append(newRow)
        caption = f'Errors of K-Edge Ionization Energies (in eV) compared to experimental data'
        label = f'summary-series-wide'
        path = tools.join_path(PRIVATE.BASE_PATH+['analysis', 'paper', 'output', 'tables', 'methodsummary', f'series-summary'])
        header = [topheader, lowheader]
        positioning = 'l ' * len(lowheader)
        self.latT.export_table(caption, label, positioning, header, newBody, path)

    def separate_in_two(self):
        atomCombos = (('C',), ('N',), ('O'), ('F'), ('C', 'N', 'O', 'F'))
        header, atomsBodies = self.evalObj.summary_of_deviations(atomCombos)
        atomToBody = {}
        for atoms, body in atomsBodies:
            if len(atoms) != 1: continue
            atomToBody[atoms[0]] = body

        topheader = ['Basis', 'Method', '\multicolumn{3}{c}{C-Series}', '\multicolumn{3}{c}{N-Series}']
        lowheader = ['', '', 'MAE', 'MaxAE', 'STD', 'MAE', 'MaxAE', 'STD']
        positioning = 'l ' * len(lowheader)
        newBody = []
        for i, row in enumerate(atomToBody['C']):
            n, o, f = atomToBody['N'], atomToBody['O'], atomToBody['F']
            newRow = row[:2] + row[3:4] + row[5:7] + n[i][3:4] + n[i][5:7]
            newBody.append(newRow)
        caption = f'Errors of K-Edge Ionization Energies (in eV) compared to experimental data'
        label = f'cn-series-wide'
        path = tools.join_path(PRIVATE.BASE_PATH+['analysis', 'paper', 'output', 'tables', 'methodsummary', f'cn-summary'])
        header = [topheader, lowheader]
        self.latT.export_table(caption, label, positioning, header, newBody, path)

        topheader = ['Basis', 'Method', '\multicolumn{3}{c}{O-Series}', '\multicolumn{3}{c}{F-Series}']
        newBody = []
        for i, row in enumerate(atomToBody['C']):
            n, o, f = atomToBody['N'], atomToBody['O'], atomToBody['F']
            newRow = row[:2] + o[i][3:4] + o[i][5:7] + f[i][3:4] + f[i][5:7]
            newBody.append(newRow)
        caption = f'Errors of K-Edge Ionization Energies (in eV) compared to experimental data'
        label = f'of-series-wide'
        path = tools.join_path(PRIVATE.BASE_PATH+['analysis', 'paper', 'output', 'tables', 'methodsummary', f'of-summary'])
        header = [topheader, lowheader]
        self.latT.export_table(caption, label, positioning, header, newBody, path)


if __name__ == "__main__":
    tables = Tables()
    tables.summary()
    # tables.combined_table()
    # tables.separate_in_two()



    