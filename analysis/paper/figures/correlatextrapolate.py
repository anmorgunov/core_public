import sys, constants
sys.path.append('../../core_excitations') # this allows us to import files from parent directories
import PRIVATE
import tools
import constants
import analysis.modules.table 
from analysis.modules.graph import Graph
import analysis.regression.evaluate
import analysis.extrapolation.extrapolations
import plotly.graph_objects as go
import plotly.express as px
import time
from plotly.subplots import make_subplots
import numpy as np
from scipy import stats
import plotly.io as pio   
from decimal import Decimal
# pio.kaleido.scope.mathjax = None

class Figures(Graph):

    def __init__(self):
        super().__init__()
        self.PAPER = PRIVATE.BASE_PATH + ['analysis', 'paper']
        self.extrapObj = analysis.extrapolation.extrapolations.export()
        self.evalObj = analysis.regression.evaluate.export()

        self.molToAtom = self.evalObj.molWithExpToAtom
        self.atomToMols = {}
        for mol, atom in self.molToAtom.items():
            self.atomToMols.setdefault(atom, set()).add(mol)
        # self.atomToMols = constants.ATOM_TO_MOLS
        self.correlate = ['#03045e', '#00b4d8', '#f72585'] #darker green 33a02c
        self.rounder = tools.rounder
        self.ccscale = ['#90caf9', '#64b5f6', '#2196f3', '#1976d2']
        self.mpscale = ['#e0aaff', '#c77dff', '#9d4edd', '#7b2cbf']

        
        self.basToMethods = {
            'D': 'UHF MP2 CCSD CCSD(T)',
            'T': 'UHF MP2 CCSD CCSD(T)',
            'Q': 'UHF MP2 CCSD CCSD(T)',
            '5': 'UHF MP2'
        }
        self.methodToBases = {}
        for basis, methods in self.basToMethods.items():
            for method in methods.split():
                self.methodToBases.setdefault(method, []).append(basis)

    def number_in_scientific(self, number, threshold=3, sigdigs=5):
        sci = '%.2E' % Decimal(str(number))
        if abs(int(sci[-3:])) <= threshold:
            return self.rounder(sigdigs)(number)
        return sci

    def calculate_r2(self, xs, ys):
        y_bar = np.mean(ys)
        residual_squares = np.sum([(y-x)**2 for x, y in zip(xs, ys)])
        total_squares = np.sum([(y-y_bar)**2 for y in ys])
        return 1 - residual_squares/total_squares

    def calculate_statistics(self, xs, ys):
        rmsd = np.sqrt(np.mean([(y-x)**2 for x, y in zip(xs, ys)]))
        mae = np.mean([abs(y-x) for x, y in zip(xs, ys)])
        return rmsd, mae

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


    def add_r_equation(self, figure, r_sq, x, y, dy=0.05, num_digs=4):
        figure.add_annotation(dict(xref='paper', yref='paper', x=x, y=y+dy,
                              xanchor='left', yanchor='top',
                              text=f"$R^2 = {self.rounder(num_digs)(r_sq)}$", #:.{num_digs}f
                              font=dict(family=self.FONT,
                                        size=self.ANNOTATION_SIZE,
                                        color=self.BLACK),
                              showarrow=False))

    def create_scatter_trace(self, xvals, yvals, mode, label, color, thick=4, showlegend=True, markerSize=10, markerBorW=1, text=''): 
        return go.Scatter(x=xvals, y=yvals, mode=mode, name=label, line=dict(color=color, width=thick, dash='dot'), 
            marker=dict(color=color, size=markerSize, line=dict(color='rgba(51, 51, 51, 1)', width=markerBorW)), text=text, connectgaps=False, showlegend=showlegend)

    def correlate_extrapolation_sandbox(self, interact=False, suffix=''):
        fig = go.Figure()
        bases = ['D', 'T', 'Q', '5']
        # atoms = ['C', 'N', 'O', 'F']
        atoms = ['O']
        names = ['MP2(T Q 5)+DifD', 'T-Q-CCSD']
        xVals, yVals = [], []
        for atom in atoms:
            for molecule in self.atomToMols[atom]:
                data = self.extrapObj.atomToMolToCBS[atom][molecule]
                xVals.append(data[names[0]])
                yVals.append(data[names[1]])
        ref = [min(xVals), max(xVals)]
        fig.add_trace(self.create_scatter_trace(ref, ref, 'lines', label='y=x', color=self.correlate[0]))
        fig.add_trace(self.create_scatter_trace(xVals, yVals, 'markers', label=' vs. '.join(names), color=self.correlate[1]))

        self._update_axes(fig, xtitle='Deviation from exp. (eV)')
        fig.update_xaxes(showgrid=True, gridcolor=self.GREY)
        self._update_fig(fig, 'Accuracy of methods')
        fig.update_layout(width=720, height=500)
        if interact:
            fig.show()
        else:
            self._save_fig(fig, 'methods_bars_all'+suffix)

    def correlate_extrapolation_within_series(self, interact=False, suffix=''):
        self.TITLE_SIZE = 12
        atoms = ['C', 'N', 'O', 'F']
        names = [["MP2(D T Q 5)+DifD", "D-T-Q-CCSD"], ["MP2(T Q 5)+DifD", "T-Q-CCSD"], ["MP2(D T Q 5)+DifD(T)", "D-T-Q-CCSD(T)"], ["MP2(T Q 5)+DifD(T)", "T-Q-CCSD(T)"]]
        # names = [["MP2(D T Q 5)+Dif3-21G", "D-T-Q-CCSD"], ["MP2(T Q 5)+Dif3-21G", "T-Q-CCSD"], ["MP2(D T Q 5)+Dif3-21G(T)", "D-T-Q-CCSD(T)"], ["MP2(T Q 5)+Dif3-21G(T)", "T-Q-CCSD(T)"]]
        f = lambda x: ' vs. '.join(x)
        fig = make_subplots(rows=2, cols=2, subplot_titles=(f(names[0]), f(names[1]), f(names[2]), f(names[3])), vertical_spacing=0.12)
        # bases = ['D', 'T', 'Q', '5']
        for j, extrapNames in enumerate(names):
            xVals, yVals = [], []
            for atom in atoms:
                atomX, atomY, mols = [], [], []
                xmin, ymin = float('inf'), float('inf')
                for molecule in self.atomToMols[atom]:
                    if molecule in self.extrapObj.smallBasisException: continue
                    data = self.extrapObj.atomToMolToCBS[atom][molecule]
                    atomX.append(x:=data[extrapNames[0]])
                    atomY.append(y:=data[extrapNames[1]])
                    mols.append(molecule)
                    if x < xmin: 
                        xmin = x
                        ymin = y
                shifter = lambda x, ref: x - ref+0.5
                shiftedX = [shifter(x, xmin) for x in atomX]
                shiftedY = [shifter(y, ymin) for y in atomY]
                xVals.extend(shiftedX)
                yVals.extend(shiftedY)
                
            ref = [min(xVals), max(xVals)]
            fig.add_trace(self.create_scatter_trace(ref, ref, 'lines', label='y=x', color=self.correlate[0], showlegend=False), row=1+j//2,col = 1+j%2)
            fig.add_trace(self.create_scatter_trace(xVals, yVals, 'markers', label=f(extrapNames), color=self.correlate[1], markerBorW=1, markerSize=8, showlegend=False, text=mols), row=1+j//2,col = 1+j%2)
            # res = stats.linregress(xVals, yVals)
            # r_sq = res.rvalue
            r_sq = self.calculate_r2(xVals, yVals)
            self.add_r_equation(fig, r_sq, (j%2)*0.55+0.28, -0.01+0.557*(1-j//2))

        self._update_axes(fig, xtitle='')
        fig.update_xaxes(showgrid=True, gridcolor=self.GREY)
        fig.update_yaxes(range=[0, 13]) 
        fig.update_xaxes(range=[0, 13]) 
        fig.update_annotations(font_size=14)
        self._update_fig(fig)
        fig.update_layout(width=600, height=570, margin=dict(l=20, r=20, t=30, b=0))
        if interact:
            fig.show()
        else:
            self._save_fig(fig, f'correlate_{atoms[0]}'+suffix)

    def correlate_extrapolation_summary(self, names, labels, nrows, ncols, width=600, height=570, interact=False, suffix='', do5=True, fRsqx=None, fRsqy=None):
        self.TITLE_SIZE = 12
        atoms = ['C', 'N', 'O', 'F']
        if do5:
            atomToNames = names['do5']
        else:
            atomToNames = names['no5']
        
        fig = make_subplots(rows=nrows, cols=ncols, subplot_titles=labels, vertical_spacing=0.12)
        for j in range(nrows*ncols):
            xVals, yVals = [], []
            for atom in atoms:
                atomX, atomY, mols = [], [], []
                xmin, ymin = float('inf'), float('inf')
                for molecule in self.atomToMols[atom]:
                    if molecule in self.extrapObj.smallBasisException: continue
                    data = self.extrapObj.atomToMolToCBS[atom][molecule]
                    atomX.append(x:=data[atomToNames[atom][j][0]])
                    atomY.append(y:=data[atomToNames[atom][j][1]])
                    mols.append(molecule)
                    if x < xmin: 
                        xmin = x
                        ymin = y
                shifter = lambda x, ref: x - ref+0.5
                shiftedX = [shifter(x, xmin) for x in atomX]
                shiftedY = [shifter(y, ymin) for y in atomY]
                xVals.extend(shiftedX)
                yVals.extend(shiftedY)
                
            ref = [min(xVals), max(xVals)]
            fig.add_trace(self.create_scatter_trace(ref, ref, 'lines', label='y=x', color=self.correlate[0], showlegend=False), row=1+j//2,col = 1+j%2)
            fig.add_trace(self.create_scatter_trace(xVals, yVals, 'markers', label='bestdata', color=self.correlate[1], markerBorW=1, markerSize=8, showlegend=False, text=mols), row=1+j//2,col = 1+j%2)
            # res = stats.linregress(xVals, yVals)
            # r_sq = res.rvalue
            r_sq = self.calculate_r2(xVals, yVals)
            rmsd, mae = self.calculate_statistics(xVals, yVals)
            f = tools.rounder(2)
            print(j, f(rmsd), f(mae))
            if fRsqx is None:
                fRsqx = lambda j : (j%2)*0.55+0.29
            if fRsqy is None:
                fRsqy = lambda j : -0.005+0.557*(1-j//2)
            self.add_r_equation(fig, r_sq, fRsqx(j), fRsqy(j))

        self._update_axes(fig, xtitle='')
        fig.update_xaxes(showgrid=False, gridcolor=self.GREY)
        fig.update_yaxes(range=[0, 13], showgrid=False) 
        fig.update_xaxes(range=[0, 13]) 
        fig.update_annotations(font_size=14)
        self._update_fig(fig)
        fig.update_layout(width=width, height=height, margin=dict(l=20, r=20, t=30, b=0))
        # axtitle = 'CEBE in eV shifted from lowest value '
        # fig['layout']['xaxis3'].update(title=axtitle)
        # fig['layout']['xaxis4'].update(title=axtitle)
        # fig['layout']['yaxis1'].update(title=axtitle)
        # fig['layout']['yaxis3'].update(title=axtitle)
        if interact:
            fig.show()
        else:
            self._save_fig(fig, f'correlate_summary'+suffix)

            
    def correlate_methods(self, atoms, interact=False, suffix=''):
        self.TITLE_SIZE = 12
        names = [["MP2(Q)", "Q-CCSD"], ["MP2(Q)", "Q-CCSD(T)"]]
        f = lambda x: ' vs. '.join(x)
        fig = make_subplots(rows=1, cols=2, subplot_titles=(f(names[0]), f(names[1])))
        for j, extrapNames in enumerate(names):
            xVals, yVals, mols = [], [], []
            for atom in atoms:
                atomX, atomY = [], []
                xmin, ymin = float('inf'), float('inf')
                for molecule in self.atomToMols[atom]:
                    data = self.extrapObj.atomToMolToCBS[atom][molecule]
                    atomX.append(x:=data[extrapNames[0]])
                    atomY.append(y:=data[extrapNames[1]])
                    mols.append(molecule)
                    if x < xmin: 
                        xmin = x
                        ymin = y
                shifter = lambda x, ref: x - ref+0.5
                shiftedX = [shifter(x, xmin) for x in atomX]
                shiftedY = [shifter(y, ymin) for y in atomY]
                xVals.extend(shiftedX)
                yVals.extend(shiftedY)
                
            ref = [min(xVals), max(xVals)]
            fig.add_trace(self.create_scatter_trace(ref, ref, 'lines', label='y=x', color=self.correlate[0], showlegend=False), row=1+j//2,col = 1+j%2)
            fig.add_trace(self.create_scatter_trace(xVals, yVals, 'markers', label=f(extrapNames), color=self.correlate[1], markerBorW=1, markerSize=8, showlegend=False, text=mols), row=1+j//2,col = 1+j%2)
            # res = stats.linregress(xVals, yVals)
            r_sq = self.calculate_r2(xVals, yVals)
            self.add_r_equation(fig, r_sq, j*0.55+0.33, 0.02)

        self._update_axes(fig, xtitle='')
        fig.update_xaxes(showgrid=False, gridcolor=self.GREY)
        fig.update_yaxes(range=[0, 13], showgrid=False) 
        fig.update_xaxes(range=[0, 13]) 
        fig.update_annotations(font_size=14)
        self._update_fig(fig)
        fig.update_layout(width=800, height=400, margin=dict(l=20, r=20, t=30, b=0))
        if interact:
            fig.show()
        else:
            self._save_fig(fig, f'correlate_methods'+suffix)

    # def correlate_extrapolation_summary(self, interact=False, suffix=''):
    #     names = [["MP2(D T Q 5)+DifD", "D-T-Q-CCSD"], ["MP2(T Q 5)+DifD", "T-Q-CCSD"], ["MP2(D T Q 5)+DifD(T)", "D-T-Q-CCSD(T)"], ["MP2(T Q 5)+DifD(T)", "T-Q-CCSD(T)"]]
    #     f = lambda x: ' vs. '.join(x)
    #     fig = make_subplots(rows=2, cols=2, subplot_titles=(f(names[0]), f(names[1]), f(names[2]), f(names[3])), vertical_spacing=0.10)
    #     bases = ['D', 'T', 'Q', '5']
    #     atoms = ['O']
    #     for j, extrapNames in enumerate(names):
    #         xVals, yVals = [], []
    #         for atom in atoms:
    #             for molecule in self.atomToMols[atom]:
    #                 data = self.extrapObj.atomToMolToCBS[atom][molecule]
    #                 xVals.append(data[extrapNames[0]])
    #                 yVals.append(data[extrapNames[1]])
    #         ref = [min(xVals), max(xVals)]
    #         fig.add_trace(self.create_scatter_trace(ref, ref, 'lines', label='y=x', color=self.correlate[0]), row=1+j//2,col = 1+j%2)
    #         fig.add_trace(self.create_scatter_trace(xVals, yVals, 'markers', label=f(extrapNames), color=self.correlate[1]), row=1+j//2,col = 1+j%2)

    #     self._update_axes(fig, xtitle='')
    #     fig.update_xaxes(showgrid=True, gridcolor=self.GREY)
    #     self._update_fig(fig, 'Accuracy of methods')
    #     fig.update_layout(width=1080, height=800)
    #     if interact:
    #         fig.show()
    #     else:
    #         self._save_fig(fig, 'correlate_{atoms[0]}'+suffix)

    def correlate_dstudy(self, dEffect, axrange=[-1, 1], interact=False, suffix='', extrapscale = None):
        self.TICK_SIZE = 12
        self.AXIS_TITLE_SIZE = 16
        self.ANNOTATION_SIZE = 8
        fig = make_subplots(rows=2, cols=2, subplot_titles=[f"<b>{t}</b>" for t in ("C-Series", "N-Series", "O-Series", "F-Series")], vertical_spacing=0.08, shared_xaxes=True, shared_yaxes=True, horizontal_spacing=0.07)
        atoms = ['C', 'N', 'O', 'F']
        # dEffect = [["MP2(D T Q 5)+DifD", "D-T-Q-CCSD", "MP2(D T Q 5)+DifD(T)", "D-T-Q-CCSD(T)"], ["MP2(T Q 5)+DifD", "T-Q-CCSD", "MP2(T Q 5)+DifD(T)", "T-Q-CCSD(T)"]]
        labels = ["with D", "without D"]
        for j, atom in enumerate(atoms):
            self.extrapObj.fill_dictionaries(self.extrapObj.all, [atom,])
            nameToErr = self.extrapObj.nameToErrors
            f = lambda x: x
            a = [f(d) for d in nameToErr[dEffect[0]]['raw']]
            b = [f(d) for d in nameToErr[dEffect[1]]['raw']]
            a, b = [], []
            for molecule in self.atomToMols[atom]:
                # if molecule in self.extrapObj.smallBasisException: continue
                data = self.extrapObj.atomToMolToCBS[atom][molecule]
                a.append(data[dEffect[0]]-self.extrapObj.parsed.molToExper[molecule])
                b.append(data[dEffect[1]]-self.extrapObj.parsed.molToExper[molecule])
            # o = stats.ttest_ind(a, b)
            # o = stats.ttest_rel(a, b)
            # o = stats.wilcoxon(a, b)
            # o = stats.mannwhitneyu(a, b)
            # ao = stats.shapiro(a)
            # bo = stats.shapiro(b)
            # print('normality test for', atom, ao.pvalue, bo.pvalue)
            # o = stats.f_oneway(a, b)
            ref = [min([min(a), min(b)]), max([max(a), max(b)])]
            fig.add_trace(self.create_scatter_trace(ref, ref, 'lines', label='y=x', color=extrapscale[3], showlegend=False, thick=2), row=1+j//2,col = 1+j%2)
            fig.add_trace(self.create_scatter_trace(a, b, 'markers', label='bestdata', color=extrapscale[1], markerBorW=1, markerSize=6, showlegend=False), row=1+j//2,col = 1+j%2)
            # fig.add_trace(self.create_scatter_trace(b, a, 'markers', label='bestdata', color=self.correlate[2], markerBorW=1, markerSize=8, showlegend=False), row=1+j//2,col = 1+j%2)
            # self.add_p_value(fig, self.number_in_scientific(o.pvalue, threshold=0), 0.01+0.55*(j%2), 0.40+0.55*(1-j//2))
         
        fig.update_layout(
            margin=dict(l=10, r=20, t=30, b=0),
            width=500, height=500
        )
        self._update_fig(fig)
        self._update_axes(fig, ytitle='Error from extrapolations<br>excluding double-zeta basis', ydtick=0.3, xtitle='Error from extrapolations<br>including double-zeta basis', xdtick=0.3)
        fig.update_yaxes(range=axrange, showgrid=False)
        fig.update_xaxes(range=axrange, showgrid=False, gridcolor=self.GREY, zeroline=True, zerolinecolor=self.BLACK)
        fig['layout']['xaxis1'].update(title='', showticklabels=False, tickvals=[])
        fig['layout']['xaxis2'].update(title='', showticklabels=False, tickvals=[])
        fig['layout']['yaxis2'].update(title='', tickvals=[])
        fig['layout']['yaxis4'].update(title='', tickvals=[])
        if interact: fig.show()
        else: self._save_fig(fig, 'correlate_dstudy'+suffix)

    def correlate_tstudy(self, tEffect, axrange=[-1, 1], interact=False, suffix='', extrapscale = None):
        self.TICK_SIZE = 12
        self.AXIS_TITLE_SIZE = 16
        self.ANNOTATION_SIZE = 8
        fig = make_subplots(rows=2, cols=2, subplot_titles=[f"<b>{t}</b>" for t in ("C-Series", "N-Series", "O-Series", "F-Series")], vertical_spacing=0.08, shared_xaxes=True, shared_yaxes=True, horizontal_spacing=0.07)
        atoms = ['C', 'N', 'O', 'F']
        # dEffect = [["MP2(D T Q 5)+DifD", "D-T-Q-CCSD", "MP2(D T Q 5)+DifD(T)", "D-T-Q-CCSD(T)"], ["MP2(T Q 5)+DifD", "T-Q-CCSD", "MP2(T Q 5)+DifD(T)", "T-Q-CCSD(T)"]]
        labels = ["without (T)", "with (T)"]
        for j, atom in enumerate(atoms):
            self.extrapObj.fill_dictionaries(self.extrapObj.all, [atom,])
            nameToErr = self.extrapObj.nameToErrors
            f = lambda x: x
            a = [f(d) for d in nameToErr[tEffect[0]]['raw']]
            b = [f(d) for d in nameToErr[tEffect[1]]['raw']]
            # print(atom, 'raw', len(a), len(b))
            a, b = [], []
            for molecule in self.atomToMols[atom]:
                # if molecule in self.extrapObj.smallBasisException: continue
                data = self.extrapObj.atomToMolToCBS[atom][molecule]
                a.append(data[tEffect[0]]-self.extrapObj.parsed.molToExper[molecule])
                b.append(data[tEffect[1]]-self.extrapObj.parsed.molToExper[molecule])
            # print(atom, 'molwise', len(a), len(b))
            # o = stats.ttest_ind(a, b)
            # o = stats.ttest_rel(a, b)
            # o = stats.wilcoxon(a, b)
            # o = stats.mannwhitneyu(a, b)
            # ao = stats.shapiro(a)
            # bo = stats.shapiro(b)
            # print('normality test for', atom, ao.pvalue, bo.pvalue)
            # o = stats.ttest_rel(a, b)
            # o = stats.levene(a, b)
            # o = stats.bartlett(a, b)
            # o = stats.f_oneway(a, b)
            ref = [min([min(a), min(b)]), max([max(a), max(b)])-0.1]
            fig.add_trace(self.create_scatter_trace(ref, ref, 'lines', label='y=x', color=extrapscale[3], showlegend=False, thick=2), row=1+j//2,col = 1+j%2)
            fig.add_trace(self.create_scatter_trace(a, b, 'markers', label='bestdata', color=extrapscale[1], markerBorW=1, markerSize=6, showlegend=False), row=1+j//2,col = 1+j%2)
            # fig.add_trace(self.create_scatter_trace(b, a, 'markers', label='bestdata', color=self.correlate[2], markerBorW=1, markerSize=8, showlegend=False), row=1+j//2,col = 1+j%2)
            # self.add_p_value(fig, self.number_in_scientific(o.pvalue, threshold=0, sigdigs=4), 0.01+0.55*(j%2), 0.40+0.55*(1-j//2))
         
        fig.update_layout(
            margin=dict(l=10, r=20, t=30, b=0),
            width=500, height=500
        )
        self._update_fig(fig)
        self._update_axes(fig, ytitle='Error from extrapolations<br>involving CCSD(T)', ydtick=0.3, xtitle='Error from extrapolations<br>involving CCSD', xdtick=0.3)
        fig.update_yaxes(range=axrange, showgrid=False)
        fig.update_xaxes(range=axrange, showgrid=False, gridcolor=self.GREY, zeroline=True, zerolinecolor=self.BLACK)
        fig['layout']['xaxis1'].update(title='', showticklabels=False, tickvals=[])
        fig['layout']['xaxis2'].update(title='', showticklabels=False, tickvals=[])
        fig['layout']['yaxis2'].update(title='', tickvals=[])
        fig['layout']['yaxis4'].update(title='', tickvals=[])
        if interact: fig.show()
        else: self._save_fig(fig, 'correlate_tstudy'+suffix)


    def summary(self):
        figure="test.pdf"
        fig=px.scatter(x=[0, 1, 2, 3, 4], y=[0, 1, 4, 9, 16])
        fig.write_image(figure, format="pdf")
        time.sleep(2)
        # self.correlate_extrapolation_sandbox(interact=True)
        atomSets = [['C'], ['N'], ['O'], ['F']]
        # for atoms in atomSets:
        #     self.correlate_extrapolation_within_series(atoms, interact=False)

        # self.correlate_extrapolation_within_series(suffix='-all')
        # self.correlate_extrapolation_within_series(['C', 'N', 'O', 'F'], suffix='-all-321G')
        names = {'do5':{
                    'C': [["MP2(D T Q 5)+DifD", "D-T-Q-CCSD"], ["MP2(D T Q 5)+Dif3-21G", "D-T-Q-CCSD"], ["MP2(D T Q 5)+DifD(T)", "D-T-Q-CCSD(T)"], ["MP2(D T Q 5)+Dif3-21G(T)", "D-T-Q-CCSD(T)"]],
                    # 'C': [["MP2(T Q 5)+DifD", "T-Q-CCSD"], ["MP2(T Q 5)+Dif3-21G", "T-Q-CCSD"], ["MP2(T Q 5)+DifD(T)", "T-Q-CCSD(T)"], ["MP2(T Q 5)+Dif3-21G(T)", "T-Q-CCSD(T)"]],
                    # 'N': [["MP2(D T Q 5)+DifD", "D-T-Q-CCSD"], ["MP2(D T Q 5)+Dif3-21G", "D-T-Q-CCSD"], ["MP2(D T Q 5)+DifD(T)", "D-T-Q-CCSD(T)"], ["MP2(D T Q 5)+Dif3-21G(T)", "D-T-Q-CCSD(T)"]],
                    'N': [["MP2(T Q 5)+DifD", "T-Q-CCSD"], ["MP2(T Q 5)+Dif3-21G", "T-Q-CCSD"], ["MP2(T Q 5)+DifD(T)", "T-Q-CCSD(T)"], ["MP2(T Q 5)+Dif3-21G(T)", "T-Q-CCSD(T)"]],
                    'O': [["MP2(T Q 5)+DifD", "T-Q-CCSD"], ["MP2(T Q 5)+Dif3-21G", "T-Q-CCSD"], ["MP2(T Q 5)+DifD(T)", "T-Q-CCSD(T)"], ["MP2(T Q 5)+Dif3-21G(T)", "T-Q-CCSD(T)"]],
                    'F': [["MP2(T Q 5)+DifD", "T-Q-CCSD"], ["MP2(T Q 5)+Dif3-21G", "T-Q-CCSD"], ["MP2(T Q 5)+DifD(T)", "T-Q-CCSD(T)"], ["MP2(T Q 5)+Dif3-21G(T)", "T-Q-CCSD(T)"]],
                    },
                'no5': {
                    'C': [["MP2(D T Q)+DifD", "D-T-Q-CCSD"], ["MP2(D T Q)+Dif3-21G", "D-T-Q-CCSD"], ["MP2(D T Q)+DifD(T)", "D-T-Q-CCSD(T)"], ["MP2(D T Q)+Dif3-21G(T)", "D-T-Q-CCSD(T)"]],
                    'N': [["MP2(D T Q)+DifD", "D-T-Q-CCSD"], ["MP2(D T Q)+Dif3-21G", "D-T-Q-CCSD"], ["MP2(D T Q)+DifD(T)", "D-T-Q-CCSD(T)"], ["MP2(D T Q)+Dif3-21G(T)", "D-T-Q-CCSD(T)"]],
                    'O': [["MP2(T Q)+DifD", "T-Q-CCSD"], ["MP2(T Q)+Dif3-21G", "T-Q-CCSD"], ["MP2(T Q)+DifD(T)", "T-Q-CCSD(T)"], ["MP2(T Q)+Dif3-21G(T)", "T-Q-CCSD(T)"]],
                    'F': [["MP2(T Q)+DifD", "T-Q-CCSD"], ["MP2(T Q)+Dif3-21G", "T-Q-CCSD"], ["MP2(T Q)+DifD(T)", "T-Q-CCSD(T)"], ["MP2(T Q)+Dif3-21G(T)", "T-Q-CCSD(T)"]],
                }}
        labels = ['MP2[∞]+DifD vs. ∞-CCSD', 'MP2[∞]+Dif3-21G vs. ∞-CCSD', 'MP2[∞]+DifD(T) vs. ∞-CCSD(T)', 'MP2[∞]+Dif3-21G(T) vs. ∞-CCSD(T)']
        # names = {'do5':{
        #             'C': [["MP2(D T Q 5)+DifD(T)", "D-T-Q-CCSD(T)"], ["MP2(D T Q 5)+Dif3-21G(T)", "D-T-Q-CCSD(T)"]],
        #             'N': [["MP2(D T Q 5)+DifD(T)", "D-T-Q-CCSD(T)"], ["MP2(D T Q 5)+Dif3-21G(T)", "D-T-Q-CCSD(T)"]],
        #             'O': [["MP2(T Q 5)+DifD", "T-Q-CCSD"], ["MP2(T Q 5)+Dif3-21G", "T-Q-CCSD"]],
        #             'F': [["MP2(T Q 5)+DifD(T)", "T-Q-CCSD(T)"], ["MP2(T Q 5)+Dif3-21G(T)", "T-Q-CCSD(T)"]],
        #             },
        #         'no5': {
        #             'C': [["MP2(D T Q)+DifD(T)", "D-T-Q-CCSD(T)"], ["MP2(D T Q)+Dif3-21G(T)", "D-T-Q-CCSD(T)"]],
        #             'N': [["MP2(D T Q)+DifD(T)", "D-T-Q-CCSD(T)"], ["MP2(D T Q)+Dif3-21G(T)", "D-T-Q-CCSD(T)"]],
        #             'O': [["MP2(T Q)+DifD", "T-Q-CCSD"], ["MP2(T Q)+Dif3-21G", "T-Q-CCSD"]],
        #             'F': [["MP2(T Q)+DifD(T)", "T-Q-CCSD(T)"], ["MP2(T Q)+Dif3-21G(T)", "T-Q-CCSD(T)"]],
        #         }}
        # labels = ['MP2(∞)+DifD vs. ∞-CCSD', 'MP2(∞)+Dif3-21G vs. ∞-CCSD']
        # self.correlate_extrapolation_summary(names, labels, nrows=1, ncols=2, width=600, height=300, do5=True, suffix='-with5', fRsqx=None, fRsqy=lambda x: 0.04)
        # self.correlate_extrapolation_summary(names, labels, nrows=1, ncols=2, width=600, height=300, do5=False, suffix='-no5', fRsqx=None, fRsqy=lambda x: 0.04)
        # self.correlate_extrapolation_summary(names, labels, nrows=2, ncols=2, do5=True, suffix='-with5-detail')
        # self.correlate_extrapolation_summary(names, labels, nrows=2, ncols=2, do5=False, suffix='-no5-detail')
        # self.correlate_methods(['C', 'N', 'O', 'F'], suffix='')

        for dEffect, suffix, axrange in (
            (["D-T-Q-CCSD", "T-Q-CCSD"], '-cc', [-0.79, 0.79]), 
            (["D-T-Q-CCSD(T)", "T-Q-CCSD(T)"], '-cc-triples', [-1.09, 0.59]),):
            self.correlate_dstudy(dEffect, axrange=axrange, suffix=suffix, extrapscale=self.ccscale)

        for tEffect, suffix, axrange in (
            (["T-Q-CCSD", "T-Q-CCSD(T)"], '-cc', [-0.79, 0.79]), 
            (["D-T-Q-CCSD", "D-T-Q-CCSD(T)"], '-cc-double', [-1.09, 0.59])):
            self.correlate_tstudy(tEffect, axrange=axrange, suffix=suffix, extrapscale=self.ccscale)

        for dEffect, suffix, axrange in (
            (["MP2(D T Q)+DifD", "MP2(T Q)+DifD"], '-mp', [-0.79, 0.79]), 
            (["MP2(D T Q)+DifD(T)", "MP2(T Q)+DifD(T)"], '-mp-triples', [-1.09, 0.59])):
            self.correlate_dstudy(dEffect, axrange=axrange, suffix=suffix, extrapscale=self.mpscale)

        for tEffect, suffix, axrange in (
            (["MP2(T Q)+DifD", "MP2(T Q)+DifD(T)"], '-mp', [-0.79, 0.79]), 
            (["MP2(D T Q)+DifD", "MP2(D T Q)+DifD(T)"], '-mp-double', [-1.09, 0.59])):
            self.correlate_tstudy(tEffect, axrange=axrange, suffix=suffix, extrapscale=self.mpscale)
        pass



if __name__ == "__main__":
    figures = Figures()
    figures.summary()
    # figures.interactive()