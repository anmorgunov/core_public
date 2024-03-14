import os
from .Modules.Figure import Styler
from .Modules.Parser import BasisStatsType, AtomBasisStatsType
from .Modules.Extrapolation import SchemeStatsType, AtomSchemeStatsType
from typing import List
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np

BASIS_TO_METHOD = {
    "D": "UHF MP2 CCSD CCSD(T)".split(),
    "T": "UHF MP2 CCSD CCSD(T)".split(),
    "Q": "UHF MP2 CCSD CCSD(T)".split(),
    "5": "UHF MP2".split(),
}
METHOD_TO_BASIS = {}
for basis, methods in BASIS_TO_METHOD.items():
    for method in methods:
        METHOD_TO_BASIS.setdefault(method, []).append(basis)

METHOD_TO_LABEL = {"UHF": "ΔHF", "MP2": "ΔMP2", "CCSD": "ΔCCSD", "CCSD(T)": "ΔCCSD(T)"}
METHOD_COLORS = ["#b2df8a", "#E0AAFF", "#48cae4", "#0077b6"]  # darker green 33a02c


def create_bar_trace(xVals, yVals, label, errors, color, showlegend=True):
    return go.Bar(
        name=label,
        x=xVals,
        y=yVals,
        marker_color=color,
        error_y=dict(type="data", array=errors, color="rgb(51, 51, 51)"),
        showlegend=showlegend,
    )


def method_error_bars_general(basisStats: BasisStatsType, save_path: str) -> go.Figure:
    fig = go.Figure()
    styler = Styler()
    styler.AXIS_TITLE_SIZE = 24
    styler.TICK_SIZE = 24
    styler.LEGEND_SIZE = 22
    bases = "D T Q 5".split()

    for i, method in enumerate(METHOD_TO_BASIS):
        yVals, yErrs = [], []
        for basis in bases:
            if basis not in METHOD_TO_BASIS[method]:
                continue
            meanAE = basisStats[basis][method]["MAE"]
            std = basisStats[basis][method]["STD(AE)"]
            yVals.append(meanAE)
            yErrs.append(std)
        fig.add_trace(
            create_bar_trace(
                bases, yVals, METHOD_TO_LABEL[method], yErrs, METHOD_COLORS[i]
            )
        )
    fig.update_yaxes(range=[0, 2.19])
    styler._update_axes(
        fig,
        ytitle="MAE compared with<br>experimental CEBEs (eV)",
        xtitle="Basis Set",
        ydtick=0.2,
    )
    styler._update_fig(fig)
    fig.update_layout(barmode="group", width=720, height=490)
    fig.update_layout(
        margin=dict(l=30, r=20, t=10, b=0),
        legend=dict(
            yanchor="bottom",
            y=0.91,
            xanchor="left",
            x=0.20,
            orientation="h",
            bgcolor="rgba(0,0,0,0)",
            title_font_family=styler.FONT,
            font=dict(size=styler.LEGEND_SIZE),
        ),
    )

    styler._save_fig(fig, save_path, "methods_bars_all")
    return fig


def method_error_bars_series(
    atomToBasisStats: AtomBasisStatsType, save_path: str
) -> go.Figure:
    styler = Styler()
    styler.AXIS_TITLE_SIZE = 14
    styler.TICK_SIZE = 12
    styler.LEGEND_SIZE = 16
    fig = make_subplots(
        rows=2,
        cols=2,
        subplot_titles=[
            f"<b>{t}</b>" for t in ("C-Series", "N-Series", "O-Series", "F-Series")
        ],
        vertical_spacing=0.08,
        horizontal_spacing=0.07,
    )
    bases = "D T Q 5".split()
    atoms = "c n o f".split()
    for j, atom in enumerate(atoms):
        showLegend = True if j == 0 else False
        basisToStats = atomToBasisStats[atom]
        for i, method in enumerate(METHOD_TO_BASIS):
            yVals, yErrs = [], []
            for basis in bases:
                if basis not in METHOD_TO_BASIS[method]:
                    continue
                meanAE = basisToStats[basis][method]["MAE"]
                std = basisToStats[basis][method]["STD(AE)"]
                yVals.append(meanAE)
                yErrs.append(std)
            fig.add_trace(
                create_bar_trace(
                    bases,
                    yVals,
                    METHOD_TO_LABEL[method],
                    yErrs,
                    METHOD_COLORS[i],
                    showLegend,
                ),
                row=1 + j // 2,
                col=1 + j % 2,
            )
    fig.update_yaxes(range=[0, 2.39])
    styler._update_axes(
        fig, ytitle="MAE compared with<br>experimental CEBEs (eV)", ydtick=0.4
    )
    styler._update_fig(fig)
    fig.update_layout(barmode="group")
    fig["layout"]["yaxis2"].update(title="")
    fig["layout"]["yaxis4"].update(title="")
    fig["layout"]["xaxis1"].update(title="", tickvals=[])
    fig["layout"]["xaxis2"].update(title="", tickvals=[])
    fig.update_yaxes(ticksuffix="")
    fig.update_layout(
        margin=dict(l=0, r=20, t=30, b=20),
        legend=dict(
            yanchor="bottom",
            y=-0.15,
            xanchor="left",
            x=0.15,
            orientation="h",
            bgcolor="rgba(0,0,0,0)",
            title_font_family=styler.FONT,
            font=dict(size=styler.LEGEND_SIZE),
        ),
    )

    styler._save_fig(fig, save_path, "methods_bars_series")
    return fig


def extrap_err_bars_dtstudy_general(
    schemeStats: SchemeStatsType,
    schemeNames: List[str],
    colorscale: List[str],
    save_path: str,
    suffix: str = "",
) -> go.Figure:
    styler = Styler()
    styler.TICK_SIZE = 12
    styler.AXIS_TITLE_SIZE = 14
    fig = go.Figure()

    labels = ["D & no (T)", "D & (T)", "no D & no (T)", "no D & (T)"]
    yVals, yErrs = [], []
    for name in schemeNames:
        yVals.append(schemeStats[name]["MAE"])
        yErrs.append(schemeStats[name]["STD(AE)"])
    fig.add_trace(
        create_bar_trace(labels, yVals, "general", yErrs, colorscale, showlegend=False)
    )

    fig.update_yaxes(range=[0, 0.69])
    fig.update_layout(margin=dict(l=10, r=20, t=30, b=0), width=720, height=490)
    styler._update_axes(fig, ytitle="Error compared to exp. (eV)", ydtick=0.1)
    styler._update_fig(fig)
    styler._save_fig(fig, save_path, "bars_dtstudy_general" + suffix)
    return fig


def extrap_err_bars_dtstudy(
    atomSchemeStats: AtomSchemeStatsType,
    schemeNames: List[str],
    colorscale: List[str],
    save_path: str,
    suffix: str = "",
) -> go.Figure:
    styler = Styler()
    styler.TICK_SIZE = 12
    styler.AXIS_TITLE_SIZE = 14
    fig = make_subplots(
        rows=2,
        cols=2,
        subplot_titles=[
            f"<b>{t}</b>" for t in ("C-Series", "N-Series", "O-Series", "F-Series")
        ],
        vertical_spacing=0.08,
        horizontal_spacing=0.07,
    )
    atoms = "c n o f".split()
    labels = ["D & no (T)", "D & (T)", "no D & no (T)", "no D & (T)"]
    for j, atom in enumerate(atoms):
        nameToErr = atomSchemeStats[atom]
        yVals, yErrs = [], []
        for name in schemeNames:
            yVals.append(nameToErr[name]["MAE"])
            yErrs.append(nameToErr[name]["STD(AE)"])
        if j // 2 == 0:
            fig.add_trace(
                create_bar_trace(
                    labels, yVals, "general", yErrs, colorscale, showlegend=False
                ),
                row=1 + j // 2,
                col=1 + j % 2,
            )
        else:
            fig.add_trace(
                create_bar_trace(
                    labels, yVals, "general", yErrs, colorscale, showlegend=False
                ),
                row=1 + j // 2,
                col=1 + j % 2,
            )

    fig.update_yaxes(range=[0, 0.69])
    fig.update_layout(margin=dict(l=10, r=20, t=30, b=0), width=720, height=490)
    styler._update_axes(
        fig, ytitle="MAE compared with<br>experimental CEBEs (eV)", ydtick=0.1
    )
    styler._update_fig(fig)
    fig["layout"]["xaxis1"].update(showticklabels=False, tickvals=[])
    fig["layout"]["xaxis2"].update(showticklabels=False, tickvals=[])
    fig["layout"]["yaxis2"].update(title="")
    fig["layout"]["yaxis4"].update(title="")

    styler._save_fig(fig, save_path, "bars_dtstudy" + suffix)
    return fig


def extrap_err_bars_summary(
    atomSchemeStats: AtomSchemeStatsType, save_path: str, include5: bool = True
) -> go.Figure:
    fig = go.Figure()
    styler = Styler()

    styler.TICK_SIZE = 20
    styler.AXIS_TITLE_SIZE = 22
    if include5:
        atomToNames = {
            "c": [
                "MP2[D T Q 5]",
                "MP2[D T Q 5]+DifSTO-3G(T)",
                "MP2[D T Q 5]+Dif3-21G(T)",
                "MP2[D T Q 5]+DifD(T)",
                "D-T-Q-CCSD(T)",
            ],
            "n": [
                "MP2[T Q 5]",
                "MP2[T Q 5]+DifSTO-3G",
                "MP2[T Q 5]+Dif3-21G",
                "MP2[T Q 5]+DifD",
                "T-Q-CCSD",
            ],
            "o": [
                "MP2[T Q 5]",
                "MP2[T Q 5]+DifSTO-3G",
                "MP2[T Q 5]+Dif3-21G",
                "MP2[T Q 5]+DifD",
                "T-Q-CCSD",
            ],
            "f": [
                "MP2[T Q 5]",
                "MP2[T Q 5]+DifSTO-3G(T)",
                "MP2[T Q 5]+Dif3-21G(T)",
                "MP2[T Q 5]+DifD(T)",
                "T-Q-CCSD(T)",
            ],
        }
    else:
        atomToNames = {
            "c": [
                "MP2[D T Q]",
                "MP2[D T Q]+DifSTO-3G(T)",
                "MP2[D T Q]+Dif3-21G(T)",
                "MP2[D T Q]+DifD(T)",
                "D-T-Q-CCSD(T)",
            ],
            "n": [
                "MP2[T Q]",
                "MP2[T Q]+DifSTO-3G",
                "MP2[T Q]+Dif3-21G",
                "MP2[T Q]+DifD",
                "T-Q-CCSD",
            ],
            "o": [
                "MP2[T Q]",
                "MP2[T Q]+DifSTO-3G",
                "MP2[T Q]+Dif3-21G",
                "MP2[T Q]+DifD",
                "T-Q-CCSD",
            ],
            "f": [
                "MP2[T Q]",
                "MP2[T Q]+DifSTO-3G(T)",
                "MP2[T Q]+Dif3-21G(T)",
                "MP2[T Q]+DifD(T)",
                "T-Q-CCSD(T)",
            ],
        }

    labels = [
        "MP2[∞]",
        "MP2[∞]<br>+δ(STO-3G)",
        "MP2[∞]<br>+δ(3-21G)",
        "MP2[∞]<br>+δ(D)",
        "CCSD[∞]",
    ]
    extrapscale = ["#c77dff"] + ["#9d4edd"] * 3 + ["#48cae4"] * 1

    yVals, yErrs = [], []

    for i in range(5):  # 5 is the number of schemeNames
        abs_errs = []
        for atom, names in atomToNames.items():
            abs_errs.extend(atomSchemeStats[atom][names[i]]["abs_errors"])
        yVals.append(np.mean(abs_errs))
        yErrs.append(np.std(abs_errs))

    fig.add_trace(
        create_bar_trace(labels, yVals, "general", yErrs, extrapscale, showlegend=False)
    )

    fig.update_yaxes(range=[0, 0.59])
    styler._update_axes(
        fig, ytitle="MAE compared with<br>experimental CEBEs (eV)", ydtick=0.1
    )  # ydtick=0.1
    styler._update_fig(fig)
    fig.update_layout(
        margin=dict(l=20, r=20, t=30, b=0),
    )

    suffix = "-w5" if include5 else "-no5"
    styler._save_fig(fig, save_path, "bars_summary" + suffix)
    return fig


def small_basis_study_subplots(
    atomSchemeStats: AtomSchemeStatsType, save_path: str
) -> go.Figure:
    styler = Styler()
    styler.TICK_SIZE = 12
    styler.AXIS_TITLE_SIZE = 14
    fig = make_subplots(
        rows=2,
        cols=2,
        subplot_titles=[
            f"<b>{t}</b>" for t in ("C-Series", "N-Series", "O-Series", "F-Series")
        ],
        vertical_spacing=0.08,
        horizontal_spacing=0.07,
    )
    atomToNames = {
        "c": [
            "MP2[D T Q]+DifSTO-3G(T)",
            "MP2[D T Q]+DifSTO-6G(T)",
            "MP2[D T Q]+Dif3-21G(T)",
            "MP2[D T Q]+Dif4-31G(T)",
            "MP2[D T Q]+Dif6-31G(T)",
            "MP2[D T Q]+DifD(T)",
        ],
        "n": [
            "MP2[T Q]+DifSTO-3G",
            "MP2[T Q]+DifSTO-6G",
            "MP2[T Q]+Dif3-21G",
            "MP2[T Q]+Dif4-31G",
            "MP2[T Q]+Dif6-31G",
            "MP2[T Q]+DifD",
        ],
        "o": [
            "MP2[T Q]+DifSTO-3G",
            "MP2[T Q]+DifSTO-6G",
            "MP2[T Q]+Dif3-21G",
            "MP2[T Q]+Dif4-31G",
            "MP2[T Q]+Dif6-31G",
            "MP2[T Q]+DifD",
        ],
        "f": [
            "MP2[T Q]+DifSTO-3G(T)",
            "MP2[T Q]+DifSTO-6G(T)",
            "MP2[T Q]+Dif3-21G(T)",
            "MP2[T Q]+Dif4-31G(T)",
            "MP2[T Q]+Dif6-31G(T)",
            "MP2[T Q]+DifD(T)",
        ],
    }
    labels = ["MP2[∞]", "+STO-3G", "+STO-6G", "+3-21G", "+4-31G", "+6-31G", "+D"]
    colorscale = ["#9d4edd"] * (len(atomToNames["c"]) + 1)
    for j, atom in enumerate(atomToNames):
        xVals, yVals, yErrs = [], [], []
        for i, name in enumerate(atomToNames[atom]):
            stats = atomSchemeStats[atom][name]
            xVals.append("+" + name.split("+")[1][3:])
            yVals.append(stats["MAE"])
            yErrs.append(stats["STD(AE)"])
        fig.add_trace(
            create_bar_trace(
                labels, yVals, "oxygen", yErrs, colorscale, showlegend=False
            ),
            row=1 + j // 2,
            col=1 + j % 2,
        )

    fig.update_layout(margin=dict(l=0, r=20, t=30, b=20))
    styler._update_axes(
        fig, ytitle="MAE compared with<br>experimental CEBEs (eV)", ydtick=0.1
    )
    styler._update_fig(fig)
    fig["layout"]["xaxis1"].update(showticklabels=False, tickvals=[])
    fig["layout"]["xaxis2"].update(showticklabels=False, tickvals=[])
    fig["layout"]["yaxis1"].update(
        range=[0, 0.65],
    )
    fig["layout"]["yaxis2"].update(range=[0, 0.65], title="")
    fig["layout"]["yaxis3"].update(
        range=[0, 0.65],
    )
    fig["layout"]["yaxis4"].update(range=[0, 0.65], title="")

    styler._save_fig(fig, save_path, "mp2_small")
    return fig


def big_basis_study_subplots(
    atomSchemeStats: AtomSchemeStatsType, save_path: str
) -> go.Figure:
    styler = Styler()
    styler.TICK_SIZE = 12
    styler.AXIS_TITLE_SIZE = 14
    fig = make_subplots(
        rows=2,
        cols=2,
        subplot_titles=[
            f"<b>{t}</b>" for t in ("C-Series", "N-Series", "O-Series", "F-Series")
        ],
        vertical_spacing=0.08,
        horizontal_spacing=0.07,
    )

    atomToNames = {
        "c": [
            "MP2[D T Q]+DifD(T)",
            "MP2[D T Q 5]+DifD(T)",
            "MP2[T Q]+DifD(T)",
            "MP2[T Q 5]+DifD(T)",
            "MP2[Q 5]+DifD(T)",
            "MP2[5]+DifD(T)",
        ],
        "n": [
            "MP2[D T Q]+DifD",
            "MP2[D T Q 5]+DifD",
            "MP2[T Q]+DifD",
            "MP2[T Q 5]+DifD",
            "MP2[Q 5]+DifD",
            "MP2[5]+DifD",
        ],
        "o": [
            "MP2[D T Q]+DifD",
            "MP2[D T Q 5]+DifD",
            "MP2[T Q]+DifD",
            "MP2[T Q 5]+DifD",
            "MP2[Q 5]+DifD",
            "MP2[5]+DifD",
        ],
        "f": [
            "MP2[D T Q]+DifD(T)",
            "MP2[D T Q 5]+DifD(T)",
            "MP2[T Q]+DifD(T)",
            "MP2[T Q 5]+DifD(T)",
            "MP2[Q 5]+DifD(T)",
            "MP2[5]+DifD(T)",
        ],
    }
    colorscale = ["#9d4edd"] * (len(atomToNames["c"]) + 1)
    labels = [
        "MP2[D T Q]",
        "MP2[D T Q 5]",
        "MP2[T Q]",
        "MP2[T Q 5]",
        "MP2[Q 5]",
        "MP2[5]",
    ]
    for j, atom in enumerate(atomToNames):
        yVals, yErrs = [], []
        for i, name in enumerate(atomToNames[atom]):
            stats = atomSchemeStats[atom][name]
            
            yVals.append(stats['MAE'])
            yErrs.append(stats['STD(AE)'])
        fig.add_trace(
            create_bar_trace(
                labels, yVals, "oxygen", yErrs, colorscale, showlegend=False
            ),
            row=1 + j // 2,
            col=1 + j % 2,
        )

    fig.update_layout(margin=dict(l=0, r=20, t=30, b=20))
    styler._update_axes(
        fig, ytitle="MAE compared with<br>experimental CEBEs (eV)", ydtick=0.1
    )
    styler._update_fig(fig)
    fig["layout"]["xaxis1"].update(showticklabels=False, tickvals=[])
    fig["layout"]["xaxis2"].update(showticklabels=False, tickvals=[])
    fig["layout"]["yaxis1"].update(
        range=[0, 0.65],
    )
    fig["layout"]["yaxis2"].update(range=[0, 0.65], title="")
    fig["layout"]["yaxis3"].update(
        range=[0, 0.65],
    )
    fig["layout"]["yaxis4"].update(range=[0, 0.65], title="")

    styler._save_fig(fig, save_path, "mp2_big")
    return fig

def extrap_err_bars_for_toc(atomSchemeStats: AtomSchemeStatsType, save_path: str) -> go.Figure:
    fig = go.Figure()
    styler = Styler()

    styler.TICK_SIZE = 20
    styler.AXIS_TITLE_SIZE = 22

    atomToNames = {
        'c': ["MP2[D T Q]", "MP2[D T Q]+DifD(T)", "D-T-Q-CCSD(T)"],
        'n': ["MP2[T Q]", "MP2[T Q]+DifD", "T-Q-CCSD"],
        'o': ["MP2[T Q]", "MP2[T Q]+DifD", "T-Q-CCSD"],
        'f': ["MP2[T Q]", "MP2[T Q]+DifD(T)", "T-Q-CCSD(T)"]
    }
    labels = ['MP2[∞]', 'MP2[∞]<br>+δ(D)', 'CCSD[∞]']
    extrapscale = ['#c77dff'] + ['#9d4edd']*1 + ['#48cae4']*1 
    
    yVals, yErrs = [], []
    for i in range(3):
        abs_errs = []
        for atom, names in atomToNames.items():
            abs_errs.extend(atomSchemeStats[atom][names[i]]['abs_errors'])
        yVals.append(np.mean(abs_errs))
        yErrs.append(np.std(abs_errs))

    fig.add_trace(create_bar_trace(labels, yVals, 'general', yErrs, extrapscale, showlegend=False))

    fig.update_yaxes(range=[0, 0.59]) 
    styler._update_axes(fig, ytitle='MAE compared with<br>experimental CEBEs (eV)', ydtick=0.1)
    styler._update_fig(fig)
    fig.update_layout(
        width=500, height=500,
        margin=dict(l=20, r=20, t=30, b=0),
    )

    styler._save_fig(fig, save_path, 'bars_toc')
    return fig