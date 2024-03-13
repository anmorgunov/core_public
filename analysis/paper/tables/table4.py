import sys, Analysis.constants as constants
sys.path.append('../../core_excitations') # this allows us to import files from parent directories
import PRIVATE
import tools
import Analysis.constants as constants
import analysis.modules.table 
import analysis.regression.evaluate
import analysis.regression.regress
import analysis.extrapolation.extrapolations

class Tables:

    def __init__(self):
        self.latT = analysis.modules.table.LaTeXTable()
        self.regressObj = analysis.regression.regress.export()
        self.rounder = tools.rounder
        self.methodToBases = {
            'UHF': 'D T Q 5'.split(),
            'MP2': 'D T Q 5'.split(),
            'MP3': 'D T'.split(),
            'CCSD': 'D T Q'.split(),
            'CCSD(T)': 'D T Q'.split(),
        }
        self.method1ToBases = {
            'MP2': 'T Q 5'.split(),
            'MP3': 'T'.split(),
        }
        self.method2ToBases = {
            'CCSD': 'T Q'.split(),
            'CCSD(T)': 'T Q'.split(),
        }
        self.methodCBSToBases = {
            'MP2': 'T Q 5'.split(),
            'MP3': 'T'.split(),
            'CCSD': 'T Q'.split(),
            'CCSD(T)': 'T Q'.split(),
        }
        self.CBSnames = ['MP2(T Q 5)+DifD', 'T-Q-CCSD', 'MP2(T Q 5)+DifD(T)', 'T-Q-CCSD(T)', 'MP2(D T Q 5)+DifD', 'D-T-Q-CCSD', 'MP2(D T Q 5)+DifD(T)', 'D-T-Q-CCSD(T)']
    
    def regression_quality(self, caption, label, header, body, positioning, atom, suffix):       
        path = tools.join_path(PRIVATE.BASE_PATH+['analysis', 'paper', 'output', 'tables', 'regression', f'{atom}-regress-{suffix}'])
        self.latT.export_table(caption, label, positioning, header, body, path)
    
    def prepare_exp_regression(self, atoms):
        if len(atoms) == 1: atom = atoms[0]
        else: atom = 'all'
        caption = f'Regression Parameters for comparison with experimental values for {atom}-series'
        label = f'regression-exp-{atom.lower()}'
        f = self.rounder(4)
        body = []
        header = ['Method', 'Basis', 'Slope', 'Intercept', 'R^2']
        positioning = 'l ' * len(header)
        for method, bases in self.methodToBases.items():
            for basis in bases:
                res = self.regressObj.experimental(atoms, method, basis)[-1]
                row = [method, basis, f(res.slope), f(res.intercept), f(res.rvalue**2)]
                body.append(row)
        self.regression_quality(caption, label, header, body, positioning, atom, 'exp')

    def prepare_two_methods_regression(self, atoms):
        if len(atoms) == 1: atom = atoms[0]
        else: atom = 'all'
        caption = f'Regression Parameters for two methods comparison for {atom}-series'
        label = f'regression-two-{atom.lower()}'
        f = self.rounder(4)
        body = []
        header = ['Method(x)', 'Basis(x)', 'Basis(y)', 'Method(y)', 'Slope', 'Intercept', 'R^2']
        positioning = 'l ' * len(header)
        for method1, xbases in self.method1ToBases.items():
            for xbasis in xbases:
                for method2, ybases in self.method2ToBases.items():
                    for ybasis in ybases:
                        res = self.regressObj.two_methods(atoms, method1, method2, xbasis, ybasis)[-1]
                        row = [method1, xbasis, method2, ybasis, f(res.slope), f(res.intercept), f(res.rvalue**2)]
                        body.append(row)
        self.regression_quality(caption, label, header, body, positioning, atom, 'two')

    def prepare_cbs_regression(self, atoms):
        if len(atoms) == 1: atom = atoms[0]
        else: atom = 'all'
        caption = f'Regression Parameters for CBS comparison for {atom}-series'
        label = f'regression-two-{atom.lower()}'
        f = self.rounder(4)
        body = []
        header = ['Method', 'Basis', 'Extrapolation Scheme', 'Slope', 'Intercept', 'R^2']
        positioning = 'l ' * len(header)
        for method, bases in self.methodCBSToBases.items():
            for basis in bases:
                for name in self.CBSnames:
                    res = self.regressObj.method_to_cbs(atoms, method, name, basis)[-1]
                    row = [method, basis, name, f(res.slope), f(res.intercept), f(res.rvalue**2)]
                    body.append(row)
        self.regression_quality(caption, label, header, body, positioning, atom, 'cbs')

    def summary(self):
        atomCombos = (('C',), ('N',), ('O'), ('F'), ('C', 'N', 'O', 'F'))
        for atoms in atomCombos:
            # self.prepare_exp_regression(atoms)
            # self.prepare_two_methods_regression(atoms)
            self.prepare_cbs_regression(atoms)



if __name__ == "__main__":
    tables = Tables()
    tables.summary()
