import sys, Analysis.constants as constants
sys.path.append('../../core_excitations') # this allows us to import files from parent directories
import PRIVATE
import tools
import Analysis.constants as constants
import analysis.modules.table 
from analysis.modules.graph import Graph
import analysis.regression.evaluate
import analysis.extrapolation.extrapolations
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
from scipy import stats

class Figures(Graph):

    def __init__(self):
        super().__init__()
        self.PAPER = PRIVATE.BASE_PATH + ['analysis', 'paper']
        self.latT = analysis.modules.table.LaTeXTable()
        self.evalObj = analysis.regression.evaluate.export()
        self.extrapObj = analysis.extrapolation.extrapolations.export()
        self.molToExper = self.evalObj.molToExper
        # self.atomToMols = constants.ATOM_TO_MOLS
        self.molToAtom = self.evalObj.molWithExpToAtom
        self.data = self.evalObj.parserObj.methodToData
        self.MOL = constants.FNAME_TO_MOLS
        self.colorscale = ['#03045E', '#00B4D8', '#5A189A', '#E0AAFF', '#C77DFF']
        # self.spectrum = ['#F72585', '#B5179E', '#7209B7', '#560BAD', '#480CA8', '#3A0CA3', '#3F37C9', '#4361EE', '#4895EF', '#4CC9F0']
        self.spectrum = ['#F72585', '#B5179E', '#7209B7', '#3A0CA3', '#4361EE', '#4895EF', '#4CC9F0']
        self.rounder = tools.rounder
        self.atomToMols = {}
        for mol, atom in self.molToAtom.items():
            self.atomToMols.setdefault(atom, set()).add(mol)
    
    def create_histogram_trace(self, data, name, color, showlegend, step=0.1, opacity=0.75, normalization=''):
        return go.Histogram(x=data, name=name, xbins=dict(start=min(data)-1, end=max(data)+1, size=step),autobinx=False,
                             marker_color=color, showlegend=showlegend, opacity=opacity, histnorm=normalization) 
    
    def create_scatter_trace(self, x, y, color, thick=3, showlegend=False): # label, mode, texts, color):
        return go.Scatter(x=x, y=y, line=dict(color=color, width=thick), showlegend=showlegend)

    def find_max_bin_count(self, trace):
        x, xbins = trace['x'], trace['xbins']
        start, end, step = xbins['start'], xbins['end'], xbins['size']
        pass

    def error_bars_extrap_general(self, names, suffix, interact=False, colorscale=None):
        if colorscale is None: colorscale = self.colorscale
        atomCombos = ({'C',}, {'N',}, {'O'}, {'F'}, {'C', 'N', 'O', 'F'})
        atoms =  {'C', 'N', 'O', 'F'} 
        self.extrapObj.fill_dictionaries(self.extrapObj.all, atoms)
        nameToErr = self.extrapObj.nameToErrors
        fig = go.Figure()
        for i, name in enumerate(names):
            errdic = nameToErr[name]
            trace = self.create_histogram_trace(errdic['raw'], name, colorscale[i%len(colorscale)], showlegend=True, step=0.10, opacity=0.7)
            fig.add_trace(trace)
            mean, std = np.mean(errdic['raw']), np.std(errdic['raw'])
            xvals = np.arange(-1, 1, 0.0025)
            # self.find_max_bin_count(trace)
            normal = stats.norm(loc=mean, scale=std).pdf(xvals)
            yvals = normal * 12
            fig.add_trace(self.create_scatter_trace(xvals, yvals, color=colorscale[i%len(colorscale)]))
        fig.update_xaxes(range=[-2, 2])
        fig.update_layout(barmode='overlay', width=720, height=500)
        self._update_axes(fig, xtitle='Deviation from exp. (eV)', xdtick=0.2)
        self._update_fig(fig, )
        if interact:
            fig.show()
        else:
            self._save_fig(fig, 'extrap_histogram_all'+suffix)

    def error_bars_extrap_series(self, names, suffix, interact=False):
        atomCombos = (('C',), ('N',), ('O',), ('F',))
        self.TICK_SIZE = 12
        fig = make_subplots(rows=2, cols=2, subplot_titles=("C-Series", "N-Series", "O-Series", "F-Series"), vertical_spacing=0.15)
        scale = [4, 6, 5, 3]
        f = self.rounder(2)
        for i, atoms in enumerate(atomCombos):
            showLegend = True if i == 0 else False
            self.extrapObj.fill_dictionaries(self.extrapObj.all, atoms)
            nameToErr = self.extrapObj.nameToErrors
            for j, name in enumerate(names):
                errdic = nameToErr[name]
                # absErrs = [abs(val) for val in errdic['raw']]
                # mean, std = f(np.mean(absErrs)), f(np.std(absErrs))
                # print(mean, std)
                trace = self.create_histogram_trace(errdic['raw'], name, self.colorscale[j], showLegend, step=0.1, opacity=0.7)
                fig.add_trace(trace, row=1+i//2,col = 1+i%2)
                # fig.add_annotation(dict(xref=f'x{i+1}', yref=f'y{i+1}', x=0.8, y=0.8,
                #               xanchor='left', yanchor='top',
                #               text=f"MAE: {f(np.mean(absErrs))}\n STD: {f(np.std(absErrs))}", #:.{num_digs}f
                #               font=dict(family=self.FONT,
                #                         size=self.ANNOTATION_SIZE,
                #                         color=self.BLACK),
                #               showarrow=False))
                mean, std = np.mean(errdic['raw']), np.std(errdic['raw'])
                xvals = np.arange(-1, 1, 0.0025)
                normal = stats.norm(loc=mean, scale=std).pdf(xvals)
                yvals = normal * scale[i]
                fig.add_trace(self.create_scatter_trace(xvals, yvals, color=self.colorscale[j], thick=2), row=1+i//2,col = 1+i%2)
        fig.update_xaxes(range=[-1, 1])
        fig.update_layout(barmode='overlay', width=700, height=500, showlegend=True)
        self._update_fig(fig,) #xtitle='Deviation from exp. (eV)'
        self._update_axes(fig, xdtick=0.3)
        fig.update_layout(
            margin=dict(l=20, r=20, t=30, b=0),
        )
        if interact: fig.show()
        else: self._save_fig(fig, 'extrap_histogram_series'+suffix)

    def error_bars_extrap_series_specific(self, interact=False, suffix=''):
        from scipy.stats import gaussian_kde
        atomToNames = {
            'C': ["MP2(D T Q)+DifD(T)", "D-T-Q-CCSD(T)"],
            'N': ["MP2(T Q)+DifD", "T-Q-CCSD"],
            'O': ["MP2(T Q)+DifD", "T-Q-CCSD"],
            'F': ["MP2(T Q)+DifD(T)", "T-Q-CCSD(T)"]
        }
        self.TICK_SIZE = 12
        # fig = make_subplots(rows=2, cols=2, subplot_titles=("C-Series", "N-Series", "O-Series", "F-Series"), vertical_spacing=0.15)
        fig = go.Figure()
        errVals = [[], []]
        scale = [4, 6, 5, 3]
        f = self.rounder(2)
        labels = []
        for i, (atom, names) in enumerate(atomToNames.items()):
            showLegend = True if i == 0 else False
            self.extrapObj.fill_dictionaries(self.extrapObj.all, (atom,))
            nameToErr = self.extrapObj.nameToErrors
            for j, name in enumerate(names):
                errdic = nameToErr[name]
                errVals[j].extend(errdic['raw'])
                if i == 0:
                    labels.append(name)

        combined_data = errVals[0] + errVals[1]
        min_val, max_val = min(combined_data), max(combined_data)
        bins = list(range(int(min_val), int(max_val) + 2))
        n_bins = 10

        colors = [
            ["#5a189a", "#9d4edd"], #MP2
            ["#0077b6", "#00b4d8"] # CCSD
        ]
        for j in range(2):
            color_id = 1 if j == 0 else 4
            data = errVals[j]
            density = gaussian_kde(data)
            xs = np.linspace(np.min(data), np.max(data), 200)
            fig.add_trace(go.Histogram(x=data, marker_color=colors[j][1], xbins=dict(start=min_val, end=max_val, size=(max_val - min_val) / n_bins), opacity=0.7, histnorm="probability", name=labels[j]))
            # fig.add_trace(self.create_histogram_trace(data, name, self.colorscale[j], showLegend, step=0.1, opacity=0.7, normalization="probability"))
            fig.add_trace(self.create_scatter_trace(xs, 0.13*density(xs), color=colors[j][0], thick=4))
        fig.update_xaxes(range=[-1, 1])
        fig.update_layout(barmode='overlay', width=700, height=500, showlegend=True)
        self._update_fig(fig,) #xtitle='Deviation from exp. (eV)'
        self._update_axes(fig, xdtick=0.3, ydtick=0.1)
        fig.update_layout(
            margin=dict(l=20, r=20, t=30, b=0),
        )
        if interact: fig.show()
        else: self._save_fig(fig, 'extrap_histogram_element_specific'+suffix)

    def summary(self):
        names = 'T-Q-CCSD, MP2(T Q)+DifD, D-T-Q-CCSD, MP2(D T Q)+DifD'.split(', ')
        self.error_bars_extrap_general(names, '')
        self.error_bars_extrap_series(names, '')
        names = 'T-Q-CCSD(T), MP2(T Q)+DifD(T), D-T-Q-CCSD(T), MP2(D T Q)+DifD(T)'.split(', ')
        self.error_bars_extrap_general(names, '-triples')
        self.error_bars_extrap_series(names, '-triples')
        
    def interactive(self):
        names = self.extrapObj.all
        self.error_bars_extrap_general(names, '-combi', False, self.spectrum)




if __name__ == "__main__":
    figures = Figures()
    # figures.summary()

    figures.error_bars_extrap_series_specific()
    # figures.interactive()
