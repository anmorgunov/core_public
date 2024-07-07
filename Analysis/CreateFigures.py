from typing import List, Tuple, Union, cast

import numpy as np
import numpy.typing as npt
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from Analysis import figureSpecs
from Analysis.Modules.Extrapolation import AtomCBSType, AtomStatsType, SchemeStatsType
from Analysis.Modules.Figure import Styler
from Analysis.Modules.Parser import AtomBasisStatsType, BasisStatsType

Number = Union[int, float]

METHOD_TO_LABEL = {"UHF": "ΔHF", "MP2": "ΔMP2", "CCSD": "ΔCCSD", "CCSD(T)": "ΔCCSD(T)"}
METHOD_COLORS = ["#b2df8a", "#E0AAFF", "#48cae4", "#0077b6"]  # darker green 33a02c

SUBPLOT_TITLES = [
    f"<b>{t}</b>" for t in ("C-Series", "N-Series", "O-Series", "F-Series")
]


def _manual_delay() -> None:
    import time

    import plotly.express as px

    figure = "test.pdf"
    fig = px.scatter(x=[0, 1, 2, 3, 4], y=[0, 1, 4, 9, 16])
    fig.write_image(figure, format="pdf")
    time.sleep(2)


def create_bar_trace(
    xVals: List[str],
    yVals: List[float],
    label: str,
    errors: List[float],
    color: Union[str, List[str]],
    showlegend: bool = True,
) -> go.Bar:
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
    bases = "D T Q 5".split()
    labels = "cc-pCVDZ<br>/cc-pVDZ cc-pCVTZ<br>/cc-pVTZ cc-pCVQZ<br>/cc-pVQZ cc-pCV5Z<br>/cc-pV5Z".split()

    for i, method in enumerate(figureSpecs.METHOD_TO_BASIS):
        yVals, yErrs = [], []
        for basis in bases:
            if basis not in figureSpecs.METHOD_TO_BASIS[method]:
                continue
            meanAE = basisStats[basis][method]["MAE"]
            std = basisStats[basis][method]["STD(AE)"]
            yVals.append(meanAE)
            yErrs.append(std)
        fig.add_trace(
            create_bar_trace(
                labels, yVals, METHOD_TO_LABEL[method], yErrs, METHOD_COLORS[i]
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
    atomToBasisStats: AtomBasisStatsType,
    bases: List[str],
    x_labels: List[str],
    save_path: str,
    fname_suffix: str = "",
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
    atoms = "c n o f".split()
    for j, atom in enumerate(atoms):
        showLegend = True if j == 0 else False
        basisToStats = atomToBasisStats[atom]
        for i, method in enumerate(figureSpecs.METHOD_TO_BASIS):
            yVals, yErrs = [], []
            for basis in bases:
                if basis not in figureSpecs.METHOD_TO_BASIS[method]:
                    continue
                meanAE = basisToStats[basis][method]["MAE"]
                std = basisToStats[basis][method]["STD(AE)"]
                yVals.append(meanAE)
                yErrs.append(std)
            fig.add_trace(
                create_bar_trace(
                    x_labels,
                    yVals,
                    METHOD_TO_LABEL[method],
                    yErrs,
                    METHOD_COLORS[i],
                    showLegend,
                ),
                row=1 + j // 2,
                col=1 + j % 2,
            )
    fig.update_yaxes(range=[0, 2.25])
    styler._update_axes(
        fig, ytitle="MAE compared with<br>experimental CEBEs (eV)", ydtick=0.2
    )
    styler._update_fig(fig)
    fig.update_layout(barmode="group")
    fig["layout"]["yaxis2"].update(title="")
    fig["layout"]["yaxis4"].update(title="")
    fig["layout"]["xaxis1"].update(title="", tickvals=[])
    fig["layout"]["xaxis2"].update(title="", tickvals=[])
    fig.update_yaxes(ticksuffix="")
    fig.update_layout(
        margin=dict(l=0, r=20, t=30, b=80),
        legend=dict(
            yanchor="bottom",
            y=-0.2,
            xanchor="left",
            x=0.15,
            orientation="h",
            bgcolor="rgba(0,0,0,0)",
            title_font_family=styler.FONT,
            font=dict(size=styler.LEGEND_SIZE),
        ),
    )

    styler._save_fig(fig, save_path, "methods_bars_series" + fname_suffix)
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
    atomSchemeStats: AtomStatsType,
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
    atomSchemeStats: AtomStatsType, save_path: str, include5: bool = True
) -> go.Figure:
    fig = go.Figure()
    styler = Styler()
    styler.TICK_SIZE = 20
    styler.AXIS_TITLE_SIZE = 22
    if include5:
        atomToNames = figureSpecs.errBarrSummary_w5
    else:
        atomToNames = figureSpecs.errBarSummary_no5

    # fmt:off
    labels = ["MP2[∞]", "MP2[∞]<br>+δ(STO-3G)", "MP2[∞]<br>+δ(3-21G)", "MP2[∞]<br>+δ(D)", "CCSD[∞]"]  # fmt:on
    extrapscale = ["#c77dff"] + ["#9d4edd"] * 3 + ["#48cae4"] * 1

    yVals, yErrs = [], []

    for i in range(5):  # 5 is the number of schemeNames
        abs_errs: List[float] = []
        for atom, names in atomToNames.items():
            abs_errs.extend(atomSchemeStats[atom][names[i]]["abs_errors"])
        yVals.append(cast(float, np.mean(abs_errs)))
        yErrs.append(cast(float, np.std(abs_errs)))

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


def extrap_err_bars_atom_summary(
    atomSchemeStats: AtomStatsType, save_path: str, include5: bool = True
) -> go.Figure:
    fig = make_subplots(
        rows=2,
        cols=2,
        subplot_titles=SUBPLOT_TITLES,
        vertical_spacing=0.08,
        horizontal_spacing=0.07,
    )
    styler = Styler()
    styler.TICK_SIZE = 12
    styler.AXIS_TITLE_SIZE = 14
    if include5:
        atomToNames = figureSpecs.errBarrSummary_w5
    else:
        atomToNames = figureSpecs.errBarSummary_no5

    # fmt:off
    labels = ["MP2[∞]", "MP2[∞]<br>+δ(STO-3G)", "MP2[∞]<br>+δ(3-21G)", "MP2[∞]<br>+δ(D)", "CCSD[∞]"]  # fmt:on
    extrapscale = ["#c77dff"] + ["#9d4edd"] * 3 + ["#48cae4"] * 1

    for atom_idx, (atom, names) in enumerate(atomToNames.items()):
        yVals, yErrs = [], []
        for scheme_idx in range(5):  # 5 is the number of schemeNames
            errs = atomSchemeStats[atom][names[scheme_idx]]["abs_errors"]
            yVals.append(cast(float, np.mean(errs)))
            yErrs.append(cast(float, np.mean(errs)))

        fig.add_trace(
            create_bar_trace(
                labels, yVals, "general", yErrs, extrapscale, showlegend=False
            ),
            row=1 + atom_idx // 2,
            col=1 + atom_idx % 2,
        )

    fig.update_yaxes(range=[0, 0.59])
    styler._update_axes(
        fig, ytitle="MAE compared with<br>experimental CEBEs (eV)", ydtick=0.1
    )  # ydtick=0.1
    styler._update_fig(fig)
    fig.update_layout(margin=dict(l=20, r=20, t=30, b=0), width=850)
    fig["layout"]["xaxis1"].update(showticklabels=False, tickvals=[])
    fig["layout"]["xaxis2"].update(showticklabels=False, tickvals=[])
    fig["layout"]["yaxis2"].update(title="")
    fig["layout"]["yaxis4"].update(title="")
    suffix = "-w5" if include5 else "-no5"
    styler._save_fig(fig, save_path, "bars_atom_summary" + suffix)
    return fig


def benefit_of_extrapolation_over_special_basis_series(
    atomBasisStats: AtomBasisStatsType,
    atomSchemeStats: AtomStatsType,
    save_path: str,
    include5: bool = True,
) -> go.Figure:
    fig = make_subplots(
        rows=2,
        cols=2,
        subplot_titles=SUBPLOT_TITLES,
        vertical_spacing=0.08,
        horizontal_spacing=0.07,
    )
    styler = Styler()
    styler.TICK_SIZE = 12
    styler.AXIS_TITLE_SIZE = 14
    if include5:
        atomToNames = figureSpecs.extrapOrNot_w5
    else:
        atomToNames = figureSpecs.extrapOrNot_no5

    # fmt:off
    labels = ["MP2[∞]<br>+δ(D)", "CCSD[∞]", "CCSD[pcX-1]", "CCSD[ccX-DZ]"]  # fmt:on
    extrapscale = ["#9d4edd"] * 1 + ["#48cae4"] * 3

    for atom_idx, (atom, names) in enumerate(atomToNames.items()):
        yVals, yErrs = [], []
        for scheme_idx in range(2):  # 5 is the number of schemeNames
            errs = atomSchemeStats[atom][names[scheme_idx]]["abs_errors"]
            yVals.append(cast(float, np.mean(errs)))
            yErrs.append(cast(float, np.mean(errs)))

        for b_i, basis in enumerate(["pcX-1", "ccX-DZ"]):
            errDict = atomBasisStats[atom][basis][names[2 + b_i]]
            yVals.append(errDict["MAE"])
            yErrs.append(errDict["STD(AE)"])

        fig.add_trace(
            create_bar_trace(
                labels, yVals, "general", yErrs, extrapscale, showlegend=False
            ),
            row=1 + atom_idx // 2,
            col=1 + atom_idx % 2,
        )

    fig.update_yaxes(range=[0, 0.59])
    styler._update_axes(
        fig, ytitle="MAE compared with<br>experimental CEBEs (eV)", ydtick=0.1
    )  # ydtick=0.1
    styler._update_fig(fig)
    fig.update_layout(margin=dict(l=20, r=20, t=30, b=0), width=800)
    fig["layout"]["xaxis1"].update(showticklabels=False, tickvals=[])
    fig["layout"]["xaxis2"].update(showticklabels=False, tickvals=[])
    fig["layout"]["yaxis2"].update(title="")
    fig["layout"]["yaxis4"].update(title="")
    suffix = "-w5" if include5 else "-no5"
    styler._save_fig(fig, save_path, "bars_need_for_extrap_series" + suffix)
    return fig


def benefit_of_extrapolation_over_special_basis_overall(
    atomBasisStats: AtomBasisStatsType,
    atomSchemeStats: AtomStatsType,
    save_path: str,
    include5: bool = True,
) -> go.Figure:
    fig = go.Figure()
    styler = Styler()
    styler.TICK_SIZE = 12
    styler.AXIS_TITLE_SIZE = 14
    if include5:
        atomToNames = figureSpecs.extrapOrNot_w5
    else:
        atomToNames = figureSpecs.extrapOrNot_no5

    # fmt:off
    labels = ["MP2[∞]<br>+δ(D)", "CCSD[∞]", "CCSD[pcX-1]", "CCSD[ccX-DZ]"]  # fmt:on
    extrapscale = ["#9d4edd"] * 1 + ["#48cae4"] * 3

    yVals, yErrs = [], []
    for scheme_idx in range(2):  # 5 is the number of schemeNames
        errs: List[float] = []
        for atom, names in atomToNames.items():
            errs.extend(atomSchemeStats[atom][names[scheme_idx]]["abs_errors"])
        yVals.append(cast(float, np.mean(errs)))
        yErrs.append(cast(float, np.mean(errs)))

    for b_i, basis in enumerate(["pcX-1", "ccX-DZ"]):
        errs = []
        for atom, names in atomToNames.items():
            errs.extend(atomBasisStats[atom][basis][names[2 + b_i]]["abs_errors"])

        yVals.append(cast(float, np.mean(errs)))
        yErrs.append(cast(float, np.mean(errs)))

    fig.add_trace(
        create_bar_trace(labels, yVals, "general", yErrs, extrapscale, showlegend=False)
    )

    fig.update_yaxes(range=[0, 0.59])
    styler._update_axes(
        fig, ytitle="MAE compared with<br>experimental CEBEs (eV)", ydtick=0.1
    )  # ydtick=0.1
    styler._update_fig(fig)
    fig.update_layout(margin=dict(l=20, r=20, t=30, b=0))
    suffix = "-w5" if include5 else "-no5"
    styler._save_fig(fig, save_path, "bars_need_for_extrap_overall" + suffix)
    return fig


def small_basis_study_subplots(
    atomSchemeStats: AtomStatsType, save_path: str
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
    atomToNames = figureSpecs.smallBasisStudy
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
    atomSchemeStats: AtomStatsType, save_path: str
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

    atomToNames = figureSpecs.bigBasisStudy
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

    styler._save_fig(fig, save_path, "mp2_big")
    return fig


def extrap_err_bars_for_toc(
    atomSchemeStats: AtomStatsType, save_path: str
) -> go.Figure:
    fig = go.Figure()
    styler = Styler()

    styler.TICK_SIZE = 20
    styler.AXIS_TITLE_SIZE = 22

    atomToNames = {
        "c": ["MP2[D T Q]", "MP2[D T Q]+DifD(T)", "D-T-Q-CCSD(T)"],
        "n": ["MP2[T Q]", "MP2[T Q]+DifD", "T-Q-CCSD"],
        "o": ["MP2[T Q]", "MP2[T Q]+DifD", "T-Q-CCSD"],
        "f": ["MP2[T Q]", "MP2[T Q]+DifD(T)", "T-Q-CCSD(T)"],
    }
    labels = ["MP2[∞]", "MP2[∞]<br>+δ(D)", "CCSD[∞]"]
    extrapscale = ["#c77dff"] + ["#9d4edd"] * 1 + ["#48cae4"] * 1

    yVals: List[float] = []
    yErrs: List[float] = []
    for i in range(3):
        abs_errs: List[float] = []
        for atom, names in atomToNames.items():
            abs_errs.extend(atomSchemeStats[atom][names[i]]["abs_errors"])
        yVals.append(cast(float, np.mean(abs_errs)))
        yErrs.append(cast(float, np.std(abs_errs)))

    fig.add_trace(
        create_bar_trace(labels, yVals, "general", yErrs, extrapscale, showlegend=False)
    )

    fig.update_yaxes(range=[0, 0.59])
    styler._update_axes(
        fig, ytitle="MAE compared with<br>experimental CEBEs (eV)", ydtick=0.1
    )
    styler._update_fig(fig)
    fig.update_layout(
        width=500,
        height=500,
        margin=dict(l=20, r=20, t=30, b=0),
    )

    styler._save_fig(fig, save_path, "bars_toc")
    return fig


def create_scatter_trace(
    xvals: Union[List[float], npt.NDArray[np.float64]],
    yvals: Union[List[float], npt.NDArray[np.float64]],
    mode: str,
    label: str,
    color: str,
    thick: int = 4,
    showlegend: bool = True,
    markerSize: int = 10,
    markerBorW: int = 1,
    text: Union[str, List[str]] = "",
) -> go.Scatter:
    return go.Scatter(
        x=xvals,
        y=yvals,
        mode=mode,
        name=label,
        line=dict(color=color, width=thick, dash="dot"),
        marker=dict(
            color=color,
            size=markerSize,
            line=dict(color="rgba(51, 51, 51, 1)", width=markerBorW),
        ),
        text=text,
        connectgaps=False,
        showlegend=showlegend,
    )


def calculate_r2(xs: List[Number], ys: List[Number]) -> float:
    y_bar = np.mean(ys)
    residual_squares = cast(float, np.sum([(y - x) ** 2 for x, y in zip(xs, ys)]))
    total_squares = cast(float, np.sum([(y - y_bar) ** 2 for y in ys]))
    return 1 - residual_squares / total_squares


def calculate_statistics(xs: List[Number], ys: List[Number]) -> Tuple[float, float]:
    rmsd = np.sqrt(np.mean([(y - x) ** 2 for x, y in zip(xs, ys)]))
    mae = np.mean([abs(y - x) for x, y in zip(xs, ys)])
    return cast(float, rmsd), cast(float, mae)


CORRELATE_COLORS = ["#03045e", "#00b4d8", "#f72585"]


def correlate_extrapolation_summary(
    atomToMolCBS: AtomCBSType,
    save_path: str,
    include5: bool = True,
) -> go.Figure:
    n_rows, n_cols = 2, 2
    # fmt:off
    labels = ["MP2[∞]+DifD vs. ∞-CCSD", "MP2[∞]+Dif3-21G vs. ∞-CCSD", "MP2[∞]+DifD(T) vs. ∞-CCSD(T)", "MP2[∞]+Dif3-21G(T) vs. ∞-CCSD(T)"]  # fmt:on
    styler = Styler()
    styler.TITLE_SIZE = 12
    atoms = ["c", "n", "o", "f"]
    if include5:
        atomToNames = figureSpecs.corrSummary["do5"]
        suffix = "-with5-detail"
    else:
        atomToNames = figureSpecs.corrSummary["no5"]
        suffix = "-no5-detail"

    fig = make_subplots(
        rows=n_rows, cols=n_cols, subplot_titles=labels, vertical_spacing=0.12
    )

    for j in range(n_rows * n_cols):
        xVals, yVals = [], []
        for atom in atoms:
            atomX, atomY, mols = [], [], []
            xmin, ymin = float("inf"), float("inf")
            for molecule, nameToCBS in atomToMolCBS[atom].items():
                nameX = atomToNames[atom][j][0]
                nameY = atomToNames[atom][j][1]
                if molecule in figureSpecs.corrSmBasisException:
                    continue
                if (x := nameToCBS[nameX]["cbs+corr"]) is None:
                    raise TypeError(
                        f"Expected a number for extrapolation {nameX}, received None"
                    )
                if (y := nameToCBS[nameY]["cbs"]) is None:
                    raise TypeError(
                        f"Expected a number for extrapolation {nameY}, received None"
                    )
                atomX.append(x)
                atomY.append(y)
                mols.append(molecule)
                if x < xmin:
                    xmin = x
                    ymin = y
            shiftedX = [x - xmin + 0.5 for x in atomX]
            shiftedY = [y - ymin + 0.5 for y in atomY]
            xVals.extend(shiftedX)
            yVals.extend(shiftedY)

        ref = [min(xVals), max(xVals)]
        fig.add_trace(
            create_scatter_trace(
                ref,
                ref,
                "lines",
                label="y=x",
                color=CORRELATE_COLORS[0],
                showlegend=False,
            ),
            row=1 + j // 2,
            col=1 + j % 2,
        )

        fig.add_trace(
            create_scatter_trace(
                xVals,
                yVals,
                "markers",
                label="bestdata",
                color=CORRELATE_COLORS[1],
                markerBorW=1,
                markerSize=8,
                showlegend=False,
                text=mols,
            ),
            row=1 + j // 2,
            col=1 + j % 2,
        )
        r_sq = calculate_r2(xVals, yVals)
        rmsd, mae = calculate_statistics(xVals, yVals)
        print(f"{include5=}. {j=}. R^2: {r_sq:.4f}, RMSD: {rmsd:.3f}, MAE: {mae:.3f}")

        # the below is really only about the placement of the R^2 text
        def fRsqx(j: int) -> float:
            return (j % 2) * 0.55 + 0.29

        def fRsqy(j: int) -> float:
            return -0.005 + 0.557 * (1 - j // 2)

        styler.add_r_equation(fig, r_sq, fRsqx(j), fRsqy(j))

    styler._update_axes(fig, xtitle="")
    fig.update_xaxes(showgrid=False, gridcolor=styler.GREY)
    fig.update_yaxes(range=[0, 13], showgrid=False)
    fig.update_xaxes(range=[0, 13])
    fig.update_annotations(font_size=14)
    styler._update_fig(fig)
    fig.update_layout(width=600, height=570, margin=dict(l=20, r=20, t=30, b=0))

    styler._save_fig(fig, save_path, "correlate_summary" + suffix)
    return fig


def correlate_study(
    atomSchemeStats: AtomStatsType,
    dEffect: List[str],
    xtitle: str,
    ytitle: str,
    suffix: str,
    colors: List[str],
    save_path: str,
    axrange: List[Number] = [-1, 1],
) -> go.Figure:
    styler = Styler()
    styler.TICK_SIZE = 12
    styler.AXIS_TITLE_SIZE = 16
    styler.ANNOTATION_SIZE = 8
    fig = make_subplots(
        rows=2,
        cols=2,
        subplot_titles=[
            f"<b>{t}</b>" for t in ("C-Series", "N-Series", "O-Series", "F-Series")
        ],
        vertical_spacing=0.08,
        shared_xaxes=True,
        shared_yaxes=True,
        horizontal_spacing=0.07,
    )
    atoms = "c n o f".split()

    for j, atom in enumerate(atoms):
        nameToErr = atomSchemeStats[atom]
        a = nameToErr[dEffect[0]]["errors"]
        b = nameToErr[dEffect[1]]["errors"]

        ref = [min([min(a), min(b)]), max([max(a), max(b)])]
        fig.add_trace(
            create_scatter_trace(
                ref,
                ref,
                "lines",
                label="y=x",
                color=colors[3],
                showlegend=False,
                thick=2,
            ),
            row=1 + j // 2,
            col=1 + j % 2,
        )
        fig.add_trace(
            create_scatter_trace(
                xvals=a,
                yvals=b,
                mode="markers",
                label="bestdata",
                color=colors[1],
                markerBorW=1,
                markerSize=6,
                showlegend=False,
            ),
            row=1 + j // 2,
            col=1 + j % 2,
        )

    fig.update_layout(margin=dict(l=10, r=20, t=30, b=0), width=500, height=500)
    styler._update_fig(fig)
    styler._update_axes(fig, ytitle=ytitle, ydtick=0.3, xtitle=xtitle, xdtick=0.3)
    fig.update_yaxes(range=axrange, showgrid=False)
    fig.update_xaxes(
        range=axrange,
        showgrid=False,
        gridcolor=styler.GREY,
        zeroline=True,
        zerolinecolor=styler.BLACK,
    )
    fig["layout"]["xaxis1"].update(title="", showticklabels=False, tickvals=[])
    fig["layout"]["xaxis2"].update(title="", showticklabels=False, tickvals=[])
    fig["layout"]["yaxis2"].update(title="", tickvals=[])
    fig["layout"]["yaxis4"].update(title="", tickvals=[])

    styler._save_fig(fig, save_path, "correlate" + suffix)
    return fig
