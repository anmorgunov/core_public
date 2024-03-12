import sys, constants
sys.path.append('../../core_excitations') # this allows us to import files from parent directories
import PRIVATE
import tools
import constants
import analysis.modules.table 
import analysis.regression.evaluate
import analysis.extrapolation.extrapolations

class Tables:

    def __init__(self):
        self.latT = analysis.modules.table.LaTeXTable()
        self.evalObj = analysis.regression.evaluate.export()
        self.molToExper = self.evalObj.molToExper
        # self.atomToMols = constants.ATOM_TO_MOLS
        self.molToAtom = self.evalObj.molWithExpToAtom
        self.data = self.evalObj.parserObj.methodToData
        self.MOL = constants.FNAME_TO_MOLS
        self.basisToHeader = {
            '5': ['Molecule', '$\Delta$HF', '$\Delta$MP2', 'Experiment'],
            'Q': ['Molecule', '$\Delta$HF', '$\Delta$MP2', '$\Delta$CCSD', '$\Delta$CCSD(T)', 'Experiment'],
            'T': ['Molecule', '$\Delta$HF', '$\Delta$MP2', '$\Delta$MP3', '$\Delta$CCSD', '$\Delta$CCSD(T)', 'Experiment'],
            'D': ['Molecule', '$\Delta$HF', '$\Delta$MP2', '$\Delta$MP3', '$\Delta$CCSD', '$\Delta$CCSD(T)', 'Experiment']
        }
        self.rounder = tools.rounder
        # self.positioning = 'l l l l l l'

        self.atomToMols = {}
        for mol, atom in self.molToAtom.items():
            self.atomToMols.setdefault(atom, set()).add(mol)

    def collect_series(self, atom, basis):
        # general shape 
        # Molecule DHF DMP2 DCCSD DCCSD(T) Exp
        body = []
        for molecule in self.atomToMols[atom]:
            _data = self.data[f'mom/{atom.lower()}'][molecule][basis]
            mol = self.MOL[molecule]['latex']
            exp = self.molToExper[molecule]
            uhf, mp2 = _data['UHF'], _data['MP2']
            f = lambda x: f"${round(x, 2)}$ $({round(x-exp,2)})$" #f stands for format
            f2 = lambda x: f"{self.rounder(2)(x)}"
            f = f2
            if basis == '5':
                row = ["\ch{%s}" % mol, f(uhf), f(mp2), f2(exp)]
            elif basis == 'Q':
                ccsd, ccsdt = _data['CCSD'], _data['CCSD(T)']
                row = ["\ch{%s}" % mol, f(uhf), f(mp2), f(ccsd), f(ccsdt), f2(exp)]
            else:
                ccsd, ccsdt = _data['CCSD'], _data['CCSD(T)']
                mp3 = _data['MP3']
                row = ["\ch{%s}" % mol, f(uhf), f(mp2), f(mp3), f(ccsd), f(ccsdt), f2(exp)]
            body.append(row)
        return body
    
    def series_table(self, atom, basis):
        body = self.collect_series(atom, basis)
        caption = f'K-Edge Ionization Energies (in eV) of {atom}-Series in cc-pV{basis}Z/cc-pCV{basis}Z'
        label = f'{atom.lower()}-{basis.lower()}z'
        path = tools.join_path(PRIVATE.BASE_PATH+['analysis', 'paper', 'output', 'tables', 'singlezeta', f'{atom}-{basis}Z'])
        self.latT.export_table(caption, label, 'l '*len(self.basisToHeader[basis]), self.basisToHeader[basis], body, path)

    def all_results(self):
        for atom in 'C N O F'.split():
            for basis in 'D T Q 5'.split():
                self.series_table(atom, basis)


if __name__ == "__main__":
    tables = Tables()
    tables.all_results()



    