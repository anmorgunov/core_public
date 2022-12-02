import plotly.graph_objects as go 
import sys 
sys.path.append("../../core_excitations")
from analysis.parser import parse
import scipy.optimize
import numpy as np
import constants
import tools
import PRIVATE

class CBSLimit:
    """An object that is used to extrapolate results to CBS Limit for a given series.
    As of now, it assumes you have CCSD/CCSD(T) results for D,T,Q basis sets AND MP2 results for D,T,Q,5
    """
    
    def __init__(self):
        """Initialization

        Args:
            mols (str): a string of molecule names separated by spaces
            method (str): essentially an address to calculation results: /calculations/{method}/. In my case it's usually mom/{atom}
        """
        self.atomToMols = constants.ATOM_TO_MOLS
        self.METHOD = lambda x: f'mom/{x.lower()}'
        self.parsed = parse.perform_parsing(fillMol=False) # parse all calculation results
        self.atomToMolToCBS = {} 
        self.atomToMolToMP2 = {} 
        self.basToCoeff = {'D': 2, 'T': 3, 'Q': 4, '5': 5} # used for Helgaker-type extrapolation
        self.EXTRAPOLATION = PRIVATE.BASE_PATH + ['analysis', 'extrapolation']


    def _helgaker(self, zeta, a, b):
        """A standard helgaker extrapolation formula. Note: a, b are parameters which we seek. We pass them as parameters not because they're known,
        but because that's how scipy.optimize works (which actually finds those parameters).

        Args:
            zeta (int): an integer describing "zeta" value of basis set
            a (float): extrapolation parameter
            b (float): extrapolation parameter

        Returns:
            float: a fitted value at a given zeta level and given parameters
        """
        return a + b * (zeta**(-3))

    def extrapolate_two_basis(self, bases, method):
        """Despite the name, this function handles more than two bases.

        Args:
            bases (str): a string of "zetas" separated by space. These zetas are used for extrapolation
            method (str): a name of the method for which we do the extrapolation
        """
        name = '-'.join(bases.split()) + f'-{method}' # naming scheme: list zetas w/ hyphens and specify method
        for atom, mols in self.atomToMols.items():
            if atom not in self.atomToMolToCBS:
                self.atomToMolToCBS[atom] = {}
            molToCBS = self.atomToMolToCBS[atom]
            for mol in mols.split():
            # for mol in self.mols:
                basToData = self.parsed.methodToData[self.METHOD(atom)][mol] # get data for a molecule
                if len(bases) == 1: # trivial case when no extrapolation is needed
                    extrapolated = basToData[bases[0]][method]
                else:
                    for bas in bases.split():
                        if method not in basToData[bas]:
                            print('this should not happen', mol, bas)
                    values = [basToData[bas][method] for bas in bases.split()]
                    zetas = [self.basToCoeff[bas] for bas in bases.split()]
                    popt, pcov = scipy.optimize.curve_fit(self._helgaker, zetas, values)
                    #popt is a tuple of optimal parameters a, b defined by _helgaker.
                    extrapolated = popt[0]
                molToCBS.setdefault(mol, {})[name] = extrapolated
            # self.atomToMolToCBS[atom] = molToCBS
        
            
    def fill_dictionaries(self, names, atoms):
        """This method fills nameToErrors and molToNameToData dictionaries

        Names - the list of extrapolation methods that should be included in the dictionary.

        Args:
            names (list): a list of extrapolation methods that should be analyzed
        """
        round2 = self.rounder(2)
        nameToErr = {}
        molToNameToData = {}
        #redacted
        #end redacted
        nameToErr[name] = {}
        nameToErr[name]['raw'] = errors
        nameToErr[name]['MSE'] = round2(np.mean(errors)) #MSE - Mean Signed Error
        nameToErr[name]['MAE'] = round2(np.mean([abs(err) for err in errors])) # Mean Absolute Error
        nameToErr[name]['MaxAE'] = round2(max([abs(err) for err in errors])) # Max Absolute Error
        nameToErr[name]['MedAE'] = round2(np.median([abs(err) for err in errors])) # Median Absolute Error
        nameToErr[name]['STD(AE)'] = round2(np.std([abs(err) for err in errors])) # STD of Absolute Errors

        self.nameToErrors = nameToErr
        self.molToNameToData = molToNameToData
    
    def export_mp2_extrapolated_values(self, atoms, name):
        #redacted
        pass


    def plot_mp2_study(self):
        """A helper function used to study the effect of difference in small basis set (whether it's really improving results).
        Short answer: yes, it does improve.
        """
        #redacted
        pass

    def export_errors(self, names, atoms):
        """Exports errors as md table.

        The names used define what data will be printed in the table.
        """
        self.fill_dictionaries(names, atoms)
        headers = ['Extrapolation Scheme', 'MSE', 'MAE', 'MedAE', 'MaxAE', 'STD', 'Sample Size']
        body = []
        for name in names:
            erD = self.nameToErrors[name]
            row = [name, f"{erD['MSE']}", f"{erD['MAE']}", f"{erD['MedAE']}", f"{erD['MaxAE']}", f"{erD['STD(AE)']}", f"{len(erD['raw'])}"]
            body.append(row)
        return headers, body

    def initialize(self):
        """This function extrapolates calculation results to CBS limit 
        """
        # redacted
        pass

    def save_fig(self, figure, fname, html=True, jpg=False):
        if html: figure.write_html(tools.join_path(self.EXTRAPOLATION+['output', f"{fname}.html"]), include_plotlyjs='cdn')
        if jpg: figure.write_image(tools.join_path(self.EXTRAPOLATION+['output', f"{fname}.jpg"]), scale=4.0)

    def triples_study(self, names):
        round4 = self.rounder(4)
        self.fill_dictionaries(names)
        #self.molToCBS[mol][name]
        methods = 'D-T-Q-CCSD D-T-Q-CCSD(T)'.split()
        t_corrs = []
        t_corrs_rel = []
        t_corrs_rel_el = []
        t_corrs_el = []
        for mol in self.mols:
            t_corr = self.molToCBS[mol][methods[1]] - self.molToCBS[mol][methods[0]]
            t_corrs.append(round4(t_corr))
            # t_corrs_rel.append(round2(t_corr/self.molToCBS[mol][methods[0]]*100, 2))
            t_corrs_rel.append(t_corr/self.molToCBS[mol][methods[0]]*100)
            nel = sum(self.parsed.molToData[mol]['nelec'])
            # print(nel)
            # t_corrs_el.append(round2(t_corr/(self.molToCBS[mol][methods[0]]*nel)*100*100, 2))
            t_corrs_rel_el.append(t_corr/(self.molToCBS[mol][methods[0]]*nel)*100*50)
            t_corrs_el.append(t_corr/nel)


        fig = go.Figure()
        fig.add_trace(self.create_bar_trace(t_corrs, self.mols, 'Absolute correction'))
        fig.add_trace(self.create_bar_trace(t_corrs_rel, self.mols, 'Relative correction'))
        fig.add_trace(self.create_bar_trace(t_corrs_rel_el, self.mols, 'Relative correction (per 50 el)'))
        fig.add_trace(self.create_bar_trace(t_corrs_el, self.mols, 'Absolute correction (per el)'))
        self._update_fig(fig)
        self.save_fig(fig, f'triples-{self.ATOM}-DTQ')
        pass

    def export_extrapolations(self, names, molecules, fname):
        for molecule in molecules:
            #print(f"exp for {molecule} is {self.parsed.molToExper[molecule]}")
            for name in names:
                data = self.atomToMolToCBS['N'][molecule][name]
                print(f"{molecule} in {name}: {self.rounder(2)(data)}")


def export():
    obj = CBSLimit()
    obj.initialize()
    return obj

if __name__ == "__main__":
    obj = CBSLimit()
    obj.initialize()