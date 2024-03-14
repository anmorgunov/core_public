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
import plotly.io as pio   
# pio.kaleido.scope.mathjax = None
import plotly.express as px
import time
from decimal import Decimal

class Figures(Graph):

    def __init__(self):
        super().__init__()
        self.PAPER = PRIVATE.BASE_PATH + ['analysis', 'paper']
        self.evalObj = analysis.regression.evaluate.export()
        self.extrapObj = analysis.extrapolation.extrapolations.export()
        self.MOL = constants.FNAME_TO_MOLS
        self.methodscale = ['#b2df8a', '#E0AAFF', '#48cae4', '#0077b6'] #darker green 33a02c
        self.rounder = tools.rounder
        
        self.basToMethods = {
            'D': 'UHF MP2 CCSD CCSD(T)',
            'T': 'UHF MP2 CCSD CCSD(T)',
            'Q': 'UHF MP2 CCSD CCSD(T)',
            '5': 'UHF MP2'
        }
        self.methodToLabel = {'UHF': 'ΔHF', 'MP2': 'ΔMP2', 'CCSD': 'ΔCCSD', 'CCSD(T)': 'ΔCCSD(T)'}
        self.methodToBases = {}
        for basis, methods in self.basToMethods.items():
            for method in methods.split():
                self.methodToBases.setdefault(method, []).append(basis)

    def create_bar_trace(self, xVals, yVals, label, errors, color, showlegend=True):
        return go.Bar(name=label, x=xVals, y=yVals, marker_color=color, error_y=dict(type='data', array=errors, color=self.BLACK), showlegend=showlegend)

    def add_p_value(self, figure, pval, x, y, dy=0.05, num_digs=4):
        if 'E' in pval:
            text = "$p = %s\mathrm{e}{%s}$" % (pval[:-4], pval[-3:])
        else:
            text = "$p = %s $" % pval
        figure.add_annotation(dict(xref='paper', yref='paper', x=x, y=y+dy,
                              xanchor='left', yanchor='top',
                              text=text, #:.{num_digs}f {self.rounder(num_digs)(
                              font=dict(family=self.FONT,
                                        size=self.ANNOTATION_SIZE,
                                        color=self.BLACK),
                              showarrow=False))

    def number_in_scientific(self, number, threshold=3):
        sci = '%.2E' % Decimal(str(number))
        if abs(int(sci[-3:])) <= threshold:
            return self.rounder(threshold+2)(number)
        return sci

    def extrap_err_bars_general(self, interact=False, suffix=''):
        fig = go.Figure()
        atoms = ['C', 'N', 'O', 'F']
        names = ["MP2(D T Q 5)+DifD", "D-T-Q-CCSD", "MP2(T Q 5)+DifD", "T-Q-CCSD", "MP2(D T Q 5)+DifD(T)", "D-T-Q-CCSD(T)", "MP2(T Q 5)+DifD(T)", "T-Q-CCSD(T)"]
        extrapscale = ['#E0AAFF', '#48cae4',] # '#E0AAFF', '#48cae4', '#9d4edd', '#4895ef', '#9d4edd', '#4895ef']
        
        self.extrapObj.fill_dictionaries(self.extrapObj.all, atoms)
        nameToErr = self.extrapObj.nameToErrors
        yVals, yErrs = [], []
        for name in names:
            yVals.append(float(nameToErr[name]['MAE']))
            yErrs.append(float(nameToErr[name]['STD(AE)']))
        fig.add_trace(self.create_bar_trace(names, yVals, 'general', yErrs, extrapscale*4))

        fig.update_yaxes(range=[0, 0.6]) 
        self._update_axes(fig, ytitle='Error compared to exp. (eV)', ydtick=0.1)#ydtick=0.1
        self._update_fig(fig)
        fig.update_layout(
            margin=dict(l=20, r=20, t=30, b=0),
        )
        if interact: fig.show()
        else: self._save_fig(fig, 'extrap_bars_general'+suffix)

    def extrap_err_bars_for_toc(self, interact=False, suffix=''):
        fig = go.Figure()
        atoms = ['C', 'N', 'O', 'F']
        self.TICK_SIZE = 20
        self.AXIS_TITLE_SIZE = 22

        atomToNames = {
            'C': ["MP2(D T Q)", "MP2(D T Q)+DifD(T)", "D-T-Q-CCSD(T)"],
            'N': ["MP2(T Q)", "MP2(T Q)+DifD", "T-Q-CCSD"],
            'O': ["MP2(T Q)", "MP2(T Q)+DifD", "T-Q-CCSD"],
            'F': ["MP2(T Q)", "MP2(T Q)+DifD(T)", "T-Q-CCSD(T)"]
        }
        # names = ["MP2(T Q 5)+DifSTO-3G", "MP2(T Q 5)+DifSTO-3G(T)", "MP2(T Q 5)+Dif3-21G", "MP2(T Q 5)+Dif3-21G(T)", "MP2(T Q 5)+DifD", "MP2(T Q 5)+DifD(T)", "T-Q-CCSD", "T-Q-CCSD(T)"]
        labels = ['MP2[∞]', 'MP2[∞]<br>+δ(D)', 'CCSD[∞]']
        extrapscale = ['#c77dff'] + ['#9d4edd']*1 + ['#48cae4']*1 # '#E0AAFF', '#48cae4', '#9d4edd', '#4895ef', '#9d4edd', '#4895ef']
        
        devs = [[] for _ in range(len(atomToNames['C']))]
        for atom, names in atomToNames.items():
            self.extrapObj.fill_dictionaries(names, [atom])
            nameToErr = self.extrapObj.nameToErrors
            for i, name in enumerate(names):
                devs[i].extend([abs(err) for err in nameToErr[name]['raw']])
                # yVals.append(sum([float(nameToErr[name]['MAE']) for name in names])/len(names))
                # yErrs.append(sum([float(nameToErr[name]['STD(AE)']) for name in names])/len(names))

        yVals, yErrs = [], []
        for devList in devs:
            yVals.append(np.mean(devList))
            yErrs.append(np.std(devList))
        
        fig.add_trace(self.create_bar_trace(labels, yVals, 'general', yErrs, extrapscale, showlegend=False))

        fig.update_yaxes(range=[0, 0.59]) 
        self._update_axes(fig, ytitle='MAE compared with<br>experimental CEBEs (eV)', ydtick=0.1)#ydtick=0.1
        self._update_fig(fig)
        fig.update_layout(
            width=500, height=500,
            margin=dict(l=20, r=20, t=30, b=0),
        )
        if interact: fig.show()
        else: self._save_fig(fig, 'bars_toc'+suffix)

    

    def extrap_err_bars_dstudy(self, dEffect, interact=False, suffix='', extrapscale = ['#64b5f6', '#1976d2']):
        self.TICK_SIZE = 12
        self.AXIS_TITLE_SIZE = 14
        fig = make_subplots(rows=2, cols=2, subplot_titles=("C-Series", "N-Series", "O-Series", "F-Series"), vertical_spacing=0.10)
        atoms = ['C', 'N', 'O', 'F']
        # dEffect = [["MP2(D T Q 5)+DifD", "D-T-Q-CCSD", "MP2(D T Q 5)+DifD(T)", "D-T-Q-CCSD(T)"], ["MP2(T Q 5)+DifD", "T-Q-CCSD", "MP2(T Q 5)+DifD(T)", "T-Q-CCSD(T)"]]
        labels = ["with D", "without D"]
        for j, atom in enumerate(atoms):
            self.extrapObj.fill_dictionaries(self.extrapObj.all, [atom,])
            nameToErr = self.extrapObj.nameToErrors
            yVals, yErrs = [], []
            # for name in names:
            #     yVals.append(float(nameToErr[name]['MAE']))
            #     yErrs.append(float(nameToErr[name]['STD(AE)']))
            for names in dEffect:
                yVals.append(sum([float(nameToErr[name]['MAE']) for name in names])/len(names))
                yErrs.append(sum([float(nameToErr[name]['STD(AE)']) for name in names])/len(names))
            a, b = nameToErr[dEffect[0][0]]['raw'], nameToErr[dEffect[1][0]]['raw']
            o = stats.mannwhitneyu(a, b)
            if j//2 == 0:
                fig.add_trace(self.create_bar_trace(labels, yVals, 'general', yErrs, extrapscale, showlegend=False), row=1+j//2,col = 1+j%2)
            else:
                fig.add_trace(self.create_bar_trace(labels, yVals, 'general', yErrs, extrapscale, showlegend=False), row=1+j//2,col = 1+j%2)
            # print('about to add', atom, o.pvalue, 0.01+0.54*(j%2), 0.8+0.1*(1-j//2))
            # print(f"d study {atom} size {len(a)} stat {o} {suffix}")
            self.add_p_value(fig, self.number_in_scientific(o.pvalue, threshold=0), 0.23+0.55*(j%2), 0.39+0.55*(1-j//2))
        fig.update_yaxes(range=[0, 0.69]) 
        fig.update_layout(
            margin=dict(l=10, r=20, t=30, b=0),
            width=500, height=500
        )
        self._update_axes(fig, ytitle='Error compared to exp. (eV)', ydtick=0.1)#ydtick=0.1
        self._update_fig(fig)
        fig.update_xaxes()
        fig['layout']['xaxis1'].update(showticklabels = False)
        fig['layout']['xaxis2'].update(showticklabels = False)
        fig['layout']['yaxis2'].update(title='')
        fig['layout']['yaxis4'].update(title='')
        if interact: fig.show()
        else: self._save_fig(fig, 'bars_dstudy'+suffix)

    def extrap_err_bars_tstudy(self, includeD=False, interact=False, suffix=''):
        self.TICK_SIZE = 12
        self.AXIS_TITLE_SIZE = 14
        fig = make_subplots(rows=2, cols=2, subplot_titles=("C-Series", "N-Series", "O-Series", "F-Series"), vertical_spacing=0.10)
        atoms = ['C', 'N', 'O', 'F']
        names = ["MP2(D T Q 5)+DifD", "D-T-Q-CCSD", "MP2(T Q 5)+DifD", "T-Q-CCSD", "MP2(D T Q 5)+DifD(T)", "D-T-Q-CCSD(T)", "MP2(T Q 5)+DifD(T)", "T-Q-CCSD(T)"]
        # tEffect = [["MP2(D T Q 5)+DifD(T)", "D-T-Q-CCSD(T)", "MP2(T Q 5)+DifD(T)", "T-Q-CCSD(T)"], ["MP2(D T Q 5)+DifD", "D-T-Q-CCSD", "MP2(T Q 5)+DifD", "T-Q-CCSD"]]
        nameLibrary = {'noD':{
            'C': [["T-Q-CCSD"], ["T-Q-CCSD(T)",]],
            'N': [["T-Q-CCSD"], ["T-Q-CCSD(T)",]],
            'O': [["T-Q-CCSD"], ["T-Q-CCSD(T)",]],
            'F': [["T-Q-CCSD"], ["T-Q-CCSD(T)",]],
            },
            'doD': {
            'C': [["D-T-Q-CCSD"], ["D-T-Q-CCSD(T)",]],
            'N': [["D-T-Q-CCSD"], ["D-T-Q-CCSD(T)",]],
            'O': [["D-T-Q-CCSD"], ["D-T-Q-CCSD(T)",]],
            'F': [["D-T-Q-CCSD"], ["D-T-Q-CCSD(T)",]],
            }}
        if includeD: atomTotEffect = nameLibrary['doD']
        else: atomTotEffect = nameLibrary['noD']
        labels = ["without (T)", "with (T)"]
        extrapscale = ['#64b5f6', '#1976d2']
        for j, atom in enumerate(atoms):
            self.extrapObj.fill_dictionaries(self.extrapObj.all, [atom,])
            nameToErr = self.extrapObj.nameToErrors
            yVals, yErrs = [], []
            # for name in names:
            #     yVals.append(float(nameToErr[name]['MAE']))
            #     yErrs.append(float(nameToErr[name]['STD(AE)']))
            for names in atomTotEffect[atom]:
                yVals.append(sum([float(nameToErr[name]['MAE']) for name in names])/len(names))
                yErrs.append(sum([float(nameToErr[name]['STD(AE)']) for name in names])/len(names))
            if j//2 == 0:
                fig.add_trace(self.create_bar_trace(labels, yVals, 'general', yErrs, extrapscale, showlegend=False), row=1+j//2,col = 1+j%2)
            else:
                fig.add_trace(self.create_bar_trace(labels, yVals, 'general', yErrs, extrapscale, showlegend=False), row=1+j//2,col = 1+j%2)
            tEffect = atomTotEffect[atom]
            a, b = nameToErr[tEffect[0][0]]['raw'], nameToErr[tEffect[1][0]]['raw']
            o = stats.mannwhitneyu(a, b)
            # print(f"t study {atom} size {len(a)} stat {o} d? {includeD}")
            self.add_p_value(fig, self.number_in_scientific(o.pvalue, threshold=3), 0.23+0.55*(j%2), 0.39+0.55*(1-j//2))
        fig.update_yaxes(range=[0, 0.69]) 
        fig.update_layout(
            margin=dict(l=10, r=20, t=30, b=0),
            width=500, height=500
        )
        self._update_axes(fig, ytitle='Error compared to exp. (eV)', ydtick=0.1)#ydtick=0.1
        self._update_fig(fig)
        fig['layout']['xaxis1'].update(showticklabels = False)
        fig['layout']['xaxis2'].update(showticklabels = False)
        fig['layout']['yaxis2'].update(title='')
        fig['layout']['yaxis4'].update(title='')
        if interact: fig.show()
        else: self._save_fig(fig, 'bars_tstudy'+suffix)


    def extrap_err_bars_series(self, interact=False, suffix=''):
        self.TICK_SIZE = 12
        self.AXIS_TITLE_SIZE = 10
        fig = make_subplots(rows=2, cols=2, subplot_titles=("C-Series", "N-Series", "O-Series", "F-Series"), vertical_spacing=0.20)
        atoms = ['C', 'N', 'O', 'F']
        names = ["MP2(D T Q 5)+DifD", "D-T-Q-CCSD", "MP2(T Q 5)+DifD", "T-Q-CCSD", "MP2(D T Q 5)+DifD(T)", "D-T-Q-CCSD(T)", "MP2(T Q 5)+DifD(T)", "T-Q-CCSD(T)"]
        labels = ["MP2(DTQ5)+DifD", "CCSD(DTQ)", "MP2(TQ5)+DifD", "CCSD(TQ)", "MP2(DTQ5)+DifD(T)", "CCSD(T)(DTQ)", "MP2(TQ5)+DifD(T)", "CCSD(T)(TQ)"]
        extrapscale = ['#E0AAFF', '#48cae4',] # '#E0AAFF', '#48cae4', '#9d4edd', '#4895ef', '#9d4edd', '#4895ef']
        for j, atom in enumerate(atoms):
            self.extrapObj.fill_dictionaries(self.extrapObj.all, [atom,])
            nameToErr = self.extrapObj.nameToErrors
            yVals, yErrs = [], []
            for name in names:
                yVals.append(float(nameToErr[name]['MAE']))
                yErrs.append(float(nameToErr[name]['STD(AE)']))
            if j//2 == 0:
                fig.add_trace(self.create_bar_trace(labels, yVals, 'general', yErrs, extrapscale*4), row=1+j//2,col = 1+j%2)
            else:
                fig.add_trace(self.create_bar_trace(labels, yVals, 'general', yErrs, extrapscale*4), row=1+j//2,col = 1+j%2)

        fig.update_yaxes(range=[0, 0.7]) 
        self._update_axes(fig, ytitle='Error compared to exp. (eV)', ydtick=0.1)#ydtick=0.1
        self._update_fig(fig, 'Accuracy of methods')
        fig['layout']['xaxis1'].update(showticklabels = False)
        fig['layout']['xaxis2'].update(showticklabels = False)
        fig['layout']['yaxis2'].update(title='')
        fig['layout']['yaxis4'].update(title='')
        if interact: fig.show()
        else: self._save_fig(fig, 'extrap_bars_series'+suffix)

    def effect_of_small_basis_difference(self, atoms, formatter, yrange, interact=False, suffix='', title=None):
        fig = go.Figure()
        if title is None: title = f"{atoms[0]}-series"
        names = ['MP2(T Q 5)+DifD', 'MP2(T Q 5)+DifSTO-3G', 'MP2(T Q 5)+DifSTO-6G', 'MP2(T Q 5)+Dif3-21G', 'MP2(T Q 5)+Dif4-31G', 'MP2(T Q 5)+Dif6-31G', 'MP2(T Q 5)+Difdef2svp', 'MP2(T Q 5)+Difdef2svpd']
        colorscale = ['#9d4edd']*(len(names)+1) #, '#E0AAFF', '#48cae4', '#0077b6']#['#03045E', '#00B4D8', '#5A189A', '#E0AAFF', '#C77DFF']
        self.extrapObj.fill_dictionaries(names, atoms)
        fig = go.Figure()
        xVals, yVals, yErrs = [], [], []
        for i, name in enumerate(names):
            data = self.extrapObj.export_mp2_extrapolated_values(atoms, name)
            bigErr, extrErr = data['errorOfBig'], data['extrapErr']
            f = formatter
            if i == 0:
                absErrs = [f(err) for err in bigErr]
                xVals.append(name.split("+")[0])
                yVals.append(np.mean(absErrs))
                yErrs.append(np.std(absErrs))
            absExtr = [f(extr) for extr in extrErr]
            xVals.append("+"+name.split("+")[1])
            yVals.append(np.mean(absExtr))
            yErrs.append(np.std(absExtr))
        fig.add_trace(self.create_bar_trace(xVals, yVals, "oxygen", yErrs, colorscale, showlegend=False))

        fig.update_yaxes(range=yrange)
        # fig.update_layout(width=720, height=500)
        fig.update_layout(
            margin=dict(l=0, r=20, t=30, b=20),
            # legend=dict(yanchor="bottom", y=-0.15, xanchor="left", x=0.24, orientation='h', bgcolor='rgba(0,0,0,0)')
        )
        self._update_axes(fig, ytitle='Error compared to exp. (eV)', ydtick=0.1)
        self._update_fig(fig, title=title)
        if interact: fig.show()
        else: self._save_fig(fig, 'mp2_small'+suffix)

    def small_basis_subplots(self, formatter, interact=False, suffix='', title=''):
        self.TICK_SIZE = 12
        self.AXIS_TITLE_SIZE = 14
        fig = make_subplots(rows=2, cols=2, subplot_titles=[f"<b>{t}</b>" for t in ("C-Series", "N-Series", "O-Series", "F-Series")], vertical_spacing=0.08, horizontal_spacing=0.07)
        atoms = ['C', 'N', 'O', 'F']
        atomToNames = {
            'C': ['MP2(D T Q)+DifSTO-3G(T)', 'MP2(D T Q)+DifSTO-6G(T)', 'MP2(D T Q)+Dif3-21G(T)', 'MP2(D T Q)+Dif4-31G(T)', 'MP2(D T Q)+Dif6-31G(T)', 'MP2(D T Q)+DifD(T)'],
            'N': ['MP2(T Q)+DifSTO-3G', 'MP2(T Q)+DifSTO-6G', 'MP2(T Q)+Dif3-21G', 'MP2(T Q)+Dif4-31G', 'MP2(T Q)+Dif6-31G', 'MP2(T Q)+DifD'],
            'O': ['MP2(T Q)+DifSTO-3G', 'MP2(T Q)+DifSTO-6G', 'MP2(T Q)+Dif3-21G', 'MP2(T Q)+Dif4-31G', 'MP2(T Q)+Dif6-31G', 'MP2(T Q)+DifD' ],
            'F': ['MP2(T Q)+DifSTO-3G(T)', 'MP2(T Q)+DifSTO-6G(T)', 'MP2(T Q)+Dif3-21G(T)', 'MP2(T Q)+Dif4-31G(T)', 'MP2(T Q)+Dif6-31G(T)', 'MP2(T Q)+DifD(T)'],
            }
        labels = ['MP2[∞]', '+STO-3G', '+STO-6G', '+3-21G', '+4-31G', '+6-31G', '+D']
        colorscale = ['#9d4edd']*(len(atomToNames['C'])+1) #, '#E0AAFF', '#48cae4', '#0077b6']#['#03045E', '#00B4D8', '#5A189A', '#E0AAFF', '#C77DFF']
        for j, atom in enumerate(atoms):
            self.extrapObj.fill_dictionaries(atomToNames[atom], [atom])
            xVals, yVals, yErrs = [], [], []
            for i, name in enumerate(atomToNames[atom]):
                data = self.extrapObj.export_mp2_extrapolated_values([atom,], name)
                bigErr, extrErr = data['errorOfBig'], data['extrapErr']
                f = formatter
                if i == 0:
                    absErrs = [f(err) for err in bigErr]
                    xVals.append(name.split("+")[0])
                    yVals.append(np.mean(absErrs))
                    yErrs.append(np.std(absErrs))
                absExtr = [f(extr) for extr in extrErr]
                xVals.append("+"+name.split("+")[1][3: ])
                yVals.append(np.mean(absExtr))
                yErrs.append(np.std(absExtr))
            fig.add_trace(self.create_bar_trace(labels, yVals, "oxygen", yErrs, colorscale, showlegend=False), row=1+j//2, col=1+j%2)

        
        # fig.update_layout(width=720, height=500)
        fig.update_layout(
            margin=dict(l=0, r=20, t=30, b=20),
            # legend=dict(yanchor="bottom", y=-0.15, xanchor="left", x=0.24, orientation='h', bgcolor='rgba(0,0,0,0)')
        )
        self._update_axes(fig, ytitle='MAE compared with<br>experimental CEBEs (eV)', ydtick=0.1)
        self._update_fig(fig, title=title)
        fig['layout']['xaxis1'].update(showticklabels = False, tickvals=[])
        fig['layout']['xaxis2'].update(showticklabels = False, tickvals=[])
        fig['layout']['yaxis1'].update(range=[0, 0.65], )
        fig['layout']['yaxis2'].update(range=[0, 0.65], title='')
        fig['layout']['yaxis3'].update(range=[0, 0.65], )
        fig['layout']['yaxis4'].update(range=[0, 0.65], title='')
        if interact: fig.show()
        else: self._save_fig(fig, 'mp2_small'+suffix)

    def big_basis_subplots(self, formatter, interact=False, suffix='', title=''):
        self.TICK_SIZE = 12
        self.AXIS_TITLE_SIZE = 14
        fig = make_subplots(rows=2, cols=2, subplot_titles=[f"<b>{t}</b>" for t in ("C-Series", "N-Series", "O-Series", "F-Series")], vertical_spacing=0.08, horizontal_spacing=0.07)
        atoms = ['C', 'N', 'O', 'F']
        atomToNames = {
            'C': ['MP2(D T Q)+DifD(T)', 'MP2(D T Q 5)+DifD(T)', 'MP2(T Q)+DifD(T)', 'MP2(T Q 5)+DifD(T)', 'MP2(Q 5)+DifD(T)', 'MP2(5)+DifD(T)'],
            'N': ['MP2(D T Q)+DifD', 'MP2(D T Q 5)+DifD', 'MP2(T Q)+DifD', 'MP2(T Q 5)+DifD', 'MP2(Q 5)+DifD', 'MP2(5)+DifD'],
            'O': ['MP2(D T Q)+DifD', 'MP2(D T Q 5)+DifD', 'MP2(T Q)+DifD', 'MP2(T Q 5)+DifD', 'MP2(Q 5)+DifD', 'MP2(5)+DifD'],
            'F': ['MP2(D T Q)+DifD(T)', 'MP2(D T Q 5)+DifD(T)', 'MP2(T Q)+DifD(T)', 'MP2(T Q 5)+DifD(T)', 'MP2(Q 5)+DifD(T)', 'MP2(5)+DifD(T)'],
            }
        colorscale = ['#9d4edd']*(len(atomToNames['C'])+1) #, '#E0AAFF', '#48cae4', '#0077b6']#['#03045E', '#00B4D8', '#5A189A', '#E0AAFF', '#C77DFF']
        # labels = ['D-T-Q-MP2', 'D-T-Q-5-MP2', 'T-Q-MP2', 'T-Q-5-MP2', 'Q-5-MP2', '5-MP2']
        labels = ['MP2[D T Q]', 'MP2[D T Q 5]', 'MP2[T Q]', 'MP2[T Q 5]', 'MP2[Q 5]', 'MP2[5]']
        for j, atom in enumerate(atoms):
            self.extrapObj.fill_dictionaries(atomToNames[atom], [atom])
            xVals, yVals, yErrs = [], [], []
            for i, name in enumerate(atomToNames[atom]):
                data = self.extrapObj.export_mp2_extrapolated_values([atom,], name)
                bigErr, extrErr = data['errorOfBig'], data['extrapErr']
                f = formatter
                # if i == 0:
                #     absErrs = [f(err) for err in bigErr]
                #     xVals.append(name.split("+")[0])
                #     yVals.append(np.mean(absErrs))
                #     yErrs.append(np.std(absErrs))
                absExtr = [f(extr) for extr in extrErr]
                xVals.append(name.split("+")[0])
                yVals.append(np.mean(absExtr))
                yErrs.append(np.std(absExtr))
            fig.add_trace(self.create_bar_trace(labels, yVals, "oxygen", yErrs, colorscale, showlegend=False), row=1+j//2, col=1+j%2)

        
        # fig.update_layout(width=720, height=500)
        fig.update_layout(
            margin=dict(l=0, r=20, t=30, b=20),
            # legend=dict(yanchor="bottom", y=-0.15, xanchor="left", x=0.24, orientation='h', bgcolor='rgba(0,0,0,0)')
        )
        self._update_axes(fig, ytitle='MAE compared with<br>experimental CEBEs (eV)', ydtick=0.1)
        self._update_fig(fig, title=title)
        fig['layout']['xaxis1'].update(showticklabels = False, tickvals=[])
        fig['layout']['xaxis2'].update(showticklabels = False, tickvals=[])
        fig['layout']['yaxis1'].update(range=[0, 0.65], )
        fig['layout']['yaxis2'].update(range=[0, 0.65], title='')
        fig['layout']['yaxis3'].update(range=[0, 0.65], )
        fig['layout']['yaxis4'].update(range=[0, 0.65], title='')
        if interact: fig.show()
        else: self._save_fig(fig, 'mp2_big'+suffix)


    def summary(self):
        figure="test.pdf"
        fig=px.scatter(x=[0, 1, 2, 3, 4], y=[0, 1, 4, 9, 16])
        fig.write_image(figure, format="pdf")
        time.sleep(2)


        # for atom, yrangeAbs, yrange in zip("C N O F".split(), [[0, 1.1], [0, 0.6], [0, 0.5], [0, 0.5]], [[-0.1, 1.1], [-0.4, 0.6], [-0.4, 0.5], [-0.4, 0.5]]):
            # self.effect_of_small_basis_difference([atom], lambda x: abs(x), yrangeAbs, interact=False, suffix=f'{atom}-series')
        # #     self.effect_of_small_basis_difference([atom], lambda x: x, yrange, interact=False, suffix=f'{atom}-series-signed')

        # self.small_basis_subplots(lambda x: abs(x))
        # self.big_basis_subplots(lambda x: abs(x))


        # self.extrap_err_bars_with_bad(interact=False, suffix='-with5', do5=True)
        # self.extrap_err_bars_with_bad(interact=False, suffix='-no5', do5=False)
        # self.extrap_err_bars_for_toc(interact=False, suffix='-no5')

        # # self.extrap_err_bars_general(interact=False)
        # self.extrap_err_bars_series(interact=False)



        # atoms = ['C', 'N', 'O', 'F']
        # self.effect_of_small_basis_difference(atoms, lambda x: abs(x), [0, 0.7], interact=False, suffix=f'series', title='')
        # self.effect_of_small_basis_difference(atoms, lambda x: x, [-0.3, 0.7], interact=False, suffix=f'series-signed', title='')
        pass



if __name__ == "__main__":
    figures = Figures()
    figures.summary()
    # figures.interactive()