import scipy.stats
import plotly.graph_objects as go
import numpy as np
# import my_constants as constants
import constants
import os
import sys
sys.path.append('../../core_excitations') # this allows us to import files from parent directories
import PRIVATE
import tools
import analysis.parser.parse
import analysis.extrapolation.extrapolations
import analysis.regression.evaluate
#constants.PATH look for

class Regress:

    def __init__(self):
        self.parserObj = analysis.parser.parse.perform_parsing()
        self.extrapObj = analysis.extrapolation.extrapolations.export()
        self.evalObj = analysis.regression.evaluate.export()
        self.molToExper = self.parserObj.molToExper
        self.algorithmToData = self.parserObj.methodToData # method -> molecule -> basis -> data
        self.atomToMols = constants.ATOM_TO_MOLS
        self.molToAtom = {mol:atom for (atom, molecules) in self.atomToMols.items() for mol in molecules.split()}

    def do_linear_regression(self, xvals, yvals):
        xs, ys = np.array(xvals), np.array(yvals)
        res = scipy.stats.linregress(x=xs, y=ys)
        return res

    def prepare_exp_comparison(self, method, atoms, basis):
        xVals, yVals, labels = [], [], []
        for molecule, exp in self.molToExper.items():
            if molecule not in self.molToAtom: continue
            if (atom:=self.molToAtom[molecule]) not in atoms: continue
            yVals.append(exp)
            x = self.algorithmToData[f"mom/{atom.lower()}"][molecule][basis][method]
            xVals.append(x)
            labels.append(molecule)
        return xVals, yVals, labels

    def prepare_two_method_comparison(self, method1, method2, atoms, basis, basis2=None):
        if basis2 is None: basis2 = basis
        xVals, yVals, labels = [], [], []
        for molecule in self.molToExper.keys():
            if molecule not in self.molToAtom: continue
            if (atom:=self.molToAtom[molecule]) not in atoms: continue
            val1 = self.algorithmToData[f"mom/{atom.lower()}"][molecule][basis][method1]
            val2 = self.algorithmToData[f"mom/{atom.lower()}"][molecule][basis2][method2]
            xVals.append(val1)
            yVals.append(val2)
            labels.append(molecule)
        return xVals, yVals, labels

    def prepare_method_cbs_comparison(self, method, name, atoms, basis):
        xVals, yVals, labels = [], [], []
        for molecule in self.molToExper.keys():
            if molecule not in self.molToAtom: continue
            if (atom:=self.molToAtom[molecule]) not in atoms: continue
            val1 = self.algorithmToData[f"mom/{atom.lower()}"][molecule][basis][method]            
            val2 = self.extrapObj.atomToMolToCBS[atom][molecule][name]
            xVals.append(val1)
            yVals.append(val2)
            labels.append(molecule)
        return xVals, yVals, labels

    def do_linear_regression(self, xvals, yvals):
        """A function that determines the slope and interception of linear regression model. Launched immediately after parse_results
        """
        xs, ys = np.array(xvals), np.array(yvals)
        res = scipy.stats.linregress(x=xs, y=ys)
        return res
        # res.slope, res.intercept, res.rvalue

    def experimental(self, atoms, method, basis):
        xVals, yVals, labels = self.prepare_exp_comparison(method, atoms, basis)
        res = self.do_linear_regression(xVals, yVals)
        return xVals, yVals, labels, res

    def two_methods(self, atoms, method1, method2, basis, basis2=None):
        xVals, yVals, labels = self.prepare_two_method_comparison(method1, method2, atoms, basis, basis2)
        res = self.do_linear_regression(xVals, yVals)
        return xVals, yVals, labels, res

    def method_to_cbs(self, atoms, method, name, basis):
        xVals, yVals, labels = self.prepare_method_cbs_comparison(method, name, atoms, basis)
        res = self.do_linear_regression(xVals, yVals)
        return xVals, yVals, labels, res
    
def export():
    obj = Regress()
    return obj    

if __name__ == "__main__":
    reg = Regress()
    reg.experimental(('O',), 'MP2', 'T')
    reg.two_methods(('O',), 'MP2', 'CCSD', 'T')
    reg.method_to_cbs(('O',), 'MP2', 'MP2(T Q 5)+DifD', 'T')