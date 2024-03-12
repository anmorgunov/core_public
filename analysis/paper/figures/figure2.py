import sys, constants
sys.path.append('../../core_excitations') # this allows us to import files from parent directories
import PRIVATE
import tools
import constants
import analysis.modules.table 
from analysis.modules.graph import Graph
import analysis.regression.regress
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
from scipy import stats

class Figures(Graph):

    def __init__(self):
        Graph.__init__(self)
        self.PAPER = PRIVATE.BASE_PATH + ['analysis', 'paper']
        self.regressObj = analysis.regression.regress.export()
        self.colorscale = ['#03045E', '#00B4D8', '#5A189A', '#E0AAFF', '#C77DFF']
        self.rounder = tools.rounder
        self.CBSnames = ['MP2(T Q 5)+DifD', 'T-Q-CCSD', 'MP2(T Q 5)+DifD(T)', 'T-Q-CCSD(T)', 'MP2(D T Q 5)+DifD', 'D-T-Q-CCSD', 'MP2(D T Q 5)+DifD(T)', 'D-T-Q-CCSD(T)']
    
    def create_scatter_trace(self, xvals, yvals, mode, label, color, thick=4, showlegend=True): 
        return go.Scatter(x=xvals, y=yvals, mode=mode, name=label, line=dict(color=color, width=thick), 
            marker=dict(color=color, size=10, line=dict(color='rgba(51, 51, 51, 1)', width=1)), connectgaps=False, showlegend=showlegend)

    def predict_values(self, xvals, slope, intercept):
        return [x*slope+intercept for x in xvals]
    
    def add_r_equation(self, figure, res, x, y, dy=0.05, num_digs=4):
        slope, inter, r = res.slope, res.intercept, res.rvalue
        annotations = []
        figure.add_annotation(dict(xref='paper', yref='paper', x=x, y=y+dy,
                              xanchor='left', yanchor='top',
                              text=f"R^2 = {self.rounder(num_digs)(r**2)}", #:.{num_digs}f
                              font=dict(family=self.FONT,
                                        size=self.ANNOTATION_SIZE,
                                        color=self.BLACK),
                              showarrow=False))
        figure.add_annotation(dict(xref='paper', yref='paper', x=x, y=y,
                              xanchor='left', yanchor='top',
                              text=f"y = {self.rounder(num_digs)(slope)}x + {self.rounder(num_digs-1)(inter)}", #:.{num_digs}f
                              font=dict(family=self.FONT,
                                        size=self.ANNOTATION_SIZE,
                                        color=self.BLACK),
                              showarrow=False))
        # for annotation in annotations:
        #     figure.add_annotation(annotation)
        # figure.update_layout(annotations=annotations)

    def single_method_experimental(self, atoms, method, basis):
        fig = go.Figure()
        xVals, yVals, labels, res = self.regressObj.experimental(atoms, method, basis)
        predicted = self.predict_values(xVals, res.slope, res.intercept)
        fig.add_trace(self.create_scatter_trace(xVals, yVals, 'markers', method, self.colorscale[1]))
        fig.add_trace(self.create_scatter_trace(xVals, predicted, 'lines', 'Best Fit', self.colorscale[0]))
        self.add_r_equation(fig, res, x=0.65, y=0.1)
        # fig.update_xaxes(range=[-1, 1])
        fig.update_layout(width=600, height=500)
        self._update_axes(fig, xtitle=f"Δ{method.upper()} CEBE (eV)", ytitle="Experimental CEBE (eV)")
        self._update_fig(fig, )
        fig.show()
        # self._save_fig(fig, 'extrap_histogram_all'+suffix)

    def two_methods(self, atoms, method1, method2, basis):
        fig = go.Figure()
        xVals, yVals, labels, res = self.regressObj.two_methods(atoms, method1, method2, basis)
        predicted = self.predict_values(xVals, res.slope, res.intercept)
        fig.add_trace(self.create_scatter_trace(xVals, yVals, 'markers', method1, self.colorscale[1]))
        fig.add_trace(self.create_scatter_trace(xVals, predicted, 'lines', 'Best Fit', self.colorscale[0]))
        self.add_r_equation(fig, res, x=0.65, y=0.1)
        # fig.update_xaxes(range=[-1, 1])
        fig.update_layout(width=600, height=500)
        self._update_axes(fig, xtitle=f"Δ{method1.upper()} CEBE (eV)", ytitle=f"Δ{method2.upper()} CEBE (eV)")
        self._update_fig(fig, )
        fig.show()
        # self._save_fig(fig, 'extrap_histogram_all'+suffix)

    def method_to_cbs(self, atoms, method, name, basis):
        fig = go.Figure()
        xVals, yVals, labels, res = self.regressObj.method_to_cbs(atoms, method, name, basis)
        predicted = self.predict_values(xVals, res.slope, res.intercept)
        fig.add_trace(self.create_scatter_trace(xVals, yVals, 'markers', method, self.colorscale[1]))
        fig.add_trace(self.create_scatter_trace(xVals, predicted, 'lines', 'Best Fit', self.colorscale[0]))
        self.add_r_equation(fig, res, x=0.65, y=0.1)
        # fig.update_xaxes(range=[-1, 1])
        fig.update_layout(width=600, height=500)
        self._update_axes(fig, xtitle=f"Δ{method.upper()} CEBE (eV)", ytitle=f"Δ{name} CEBE (eV)")
        self._update_fig(fig, )
        fig.show()

    
    def single_method_experimental_series(self, method, basis, suffix):
        fig = go.Figure()
        atomCombos = (('C',), ('N'), ('O'), ('F'))
        colorscale = ['#03045E', '#00B4D8', '#5A189A', '#C77DFF']
        fig = make_subplots(rows=2, cols=2, subplot_titles=("<b>C-Series</b>", "<b>N-Series</b>", "<b>O-Series</b>", "<b>F-Series</b>"), vertical_spacing=0.15)
        self.ANNOTATION_SIZE = 12
        top, down, left, right = 0.68, 0.055, 0.22, 0.77
        locs = [(left, top), (right, top), (left, down), (right, down)]
        for i, atoms in enumerate(atomCombos):
            showLegend = True if i == 0 else False
            xVals, yVals, labels, res = self.regressObj.experimental(atoms, method, basis)
            predicted = self.predict_values(xVals, res.slope, res.intercept)
            fig.add_trace(self.create_scatter_trace(xVals, yVals, 'markers', method, self.colorscale[1], showlegend=showLegend), row=1+i//2,col = 1+i%2)
            fig.add_trace(self.create_scatter_trace(xVals, predicted, 'lines', 'Best Fit', self.colorscale[0], showlegend=showLegend), row=1+i//2,col = 1+i%2)
            self.add_r_equation(fig, res, x=locs[i][0], y=locs[i][1], dy=0.03)
            # fig.update_xaxes/

        fig.update_layout(width=700, height=600, showlegend=True)
        self._update_fig(fig, title=f"<b>Linear Regression of CEBE obtained with {method} (x-axis) to experimental values (y-axis)</b>") #xtitle='Deviation from exp. (eV)'
        self._update_axes(fig)
        # fig.show()
        self._save_fig(fig, 'regress_exp'+suffix)

    def two_methods_series(self, method1, method2, basis, basis2=None, suffix=''):
        fig = go.Figure()
        atomCombos = (('C',), ('N'), ('O'), ('F'))
        colorscale = ['#03045E', '#00B4D8', '#5A189A', '#C77DFF']
        fig = make_subplots(rows=2, cols=2, subplot_titles=("<b>C-Series</b>", "<b>N-Series</b>", "<b>O-Series</b>", "<b>F-Series</b>"), vertical_spacing=0.15) #horizontal_spacing=0.2, vertical_spacing=0.2
        self.ANNOTATION_SIZE = 12
        top, down, left, right = 0.68, 0.055, 0.22, 0.77
        locs = [(left, top), (right, top), (left, down), (right, down)]
        for i, atoms in enumerate(atomCombos):
            showLegend = True if i == 0 else False
            xVals, yVals, labels, res = self.regressObj.two_methods(atoms, method1, method2, basis, basis2)
            predicted = self.predict_values(xVals, res.slope, res.intercept)
            fig.add_trace(self.create_scatter_trace(xVals, yVals, 'markers', method1, self.colorscale[1], showlegend=showLegend), row=1+i//2,col = 1+i%2)
            fig.add_trace(self.create_scatter_trace(xVals, predicted, 'lines', 'Best Fit', self.colorscale[0], showlegend=showLegend), row=1+i//2,col = 1+i%2)
            self.add_r_equation(fig, res, x=locs[i][0], y=locs[i][1], dy=0.03)
            # suf = '' if i == 0 else f"{i+1}"
            # fig['layout'][f'xaxis']['title'] = f'Δ{method1.upper()} CEBE (eV)'
            # fig['layout'][f'yaxis']['title'] = f'Δ{method2.upper()} CEBE (eV)'

        fig.update_layout(width=700, height=600, showlegend=True)
        self._update_fig(fig, title=f"<b>Linear Regression of CEBE obtained with {method1} (x-axis) to {method2} (y-axis)</b>") #xtitle='Deviation from exp. (eV)'
        self.AXIS_TITLE_SIZE = 12
        self._update_axes(fig)# xtitle=f'Δ{method1.upper()} CEBE (eV)', ytitle=f'Δ{method2.upper()} CEBE (eV)')
        # fig.show()
        self._save_fig(fig, 'regress_two'+suffix)

    def method_to_cbs_series(self, method, name, basis, suffix):
        fig = go.Figure()
        atomCombos = (('C',), ('N'), ('O'), ('F'))
        colorscale = ['#03045E', '#00B4D8', '#5A189A', '#C77DFF']
        fig = make_subplots(rows=2, cols=2, subplot_titles=("<b>C-Series</b>", "<b>N-Series</b>", "<b>O-Series</b>", "<b>F-Series</b>"), vertical_spacing=0.15) #horizontal_spacing=0.2, vertical_spacing=0.2
        self.ANNOTATION_SIZE = 12
        top, down, left, right = 0.68, 0.055, 0.22, 0.77
        locs = [(left, top), (right, top), (left, down), (right, down)]
        for i, atoms in enumerate(atomCombos):
            showLegend = True if i == 0 else False
            xVals, yVals, labels, res = self.regressObj.method_to_cbs(atoms, method, name, basis)
            predicted = self.predict_values(xVals, res.slope, res.intercept)
            fig.add_trace(self.create_scatter_trace(xVals, yVals, 'markers', method, self.colorscale[1], showlegend=showLegend), row=1+i//2,col = 1+i%2)
            fig.add_trace(self.create_scatter_trace(xVals, predicted, 'lines', 'Best Fit', self.colorscale[0], showlegend=showLegend), row=1+i//2,col = 1+i%2)
            self.add_r_equation(fig, res, x=locs[i][0], y=locs[i][1], dy=0.03)
            # suf = '' if i == 0 else f"{i+1}"
            # fig['layout'][f'xaxis']['title'] = f'Δ{method1.upper()} CEBE (eV)'
            # fig['layout'][f'yaxis']['title'] = f'Δ{method2.upper()} CEBE (eV)'

        fig.update_layout(width=700, height=600, showlegend=True)
        self._update_fig(fig, title=f"<b>Linear Regression of CEBE obtained with {method} (x-axis) to {name} (y-axis)</b>") #xtitle='Deviation from exp. (eV)'
        self.AXIS_TITLE_SIZE = 12
        self._update_axes(fig)# xtitle=f'Δ{method1.upper()} CEBE (eV)', ytitle=f'Δ{method2.upper()} CEBE (eV)')
        # fig.show()
        self._save_fig(fig, 'regress_cbs'+suffix)

    def summary(self):
        atoms = ('O', )
        method = 'MP3'
        method2 = 'CCSD(T)'
        name = 'MP2(T Q 5)+DifD'
        basis = 'T'
        for method in 'MP2 MP3'.split():
            for basis in 'T 5'.split():
                if basis == '5' and method == 'MP3': continue
                self.single_method_experimental_series(method, basis, f"{method}({basis})")

        method = 'MP3'
        for method2 in 'CCSD CCSD(T)'.split():
            self.two_methods_series(method, method2, 'T', 'Q', method2[4:])

        method = 'MP3'
        basis = 'T'
        for name in self.CBSnames:
            self.method_to_cbs_series(method, name, basis, name)


        # self.single_method_experimental(atoms, method, basis)
        # self.single_method_experimental(atoms, method, basis)
        # self.two_methods(atoms, method, method2, basis)
        # self.method_to_cbs(atoms, method, name, basis)
        



if __name__ == "__main__":
    figures = Figures()
    figures.summary()
