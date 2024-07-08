import itertools
import pprint
from pathlib import Path
from typing import Dict

from analysis import constants, create_figures, create_tables, timing
from analysis.modules import extrapolation, parsers

pp = pprint.PrettyPrinter(indent=2)

BASE_PATH = Path(__file__).resolve().parent
DATA_PATH = BASE_PATH / "Data"
TABLE_PATH = DATA_PATH / "paper-tables"
FIGURE_PATH = DATA_PATH / "paper-figures"

# --------------- PERFORM PARSING ---------------
pobj = parsers.ParsedData(
    experimental_wb=str(DATA_PATH / "experimental.xlsx"),
    calculations_wb=str(DATA_PATH / "CEBE_Data.xlsx"),
    algorithms=["mom"],
    calculations_folder=DATA_PATH / "calculations",
)
# pobj.debug = True
pobj.process(save=True)

# EXAMPLE: pp.pprint(pobj.algoToDelta['mom']['o']['h2o']['D'])

# the calculations_wb file can become huge quite quickly, so if you want to check on results
# of certain molecules of interest, you can use the method below
# mols = set("h2o co co2".split())  # specify molecules you're interested in
# fname = "o-series"  # specify the file name to which this will be saved
# extracted_path = DATA_PATH / "extracts" / f"{fname}.xlsx"
# pobj.extract_molecules(mols, str(extracted_path))

# --------------- PERFORM EXTRAPOLATION ---------------

# # First, let's filter out the molecules we want to use
atomToMols = constants.ATOM_TO_MOLS
excludeAtomToMols = {
    "c": set("c-o ch3-c-o2h".split()),
    "n": set(
        "mnh2pyridi-n-e onh2pyridi-n-e o-n-h2pyridine m-n-h2pyridine pfpyridine".split()
    ),
    "o": set("cf3co-o-h cf3c-o-oh ch3co-o-ch3 ch3c-o-och3 ".split()),
    "f": set("cf3ocf3".split()),
}
atomToMols = {
    atom.lower(): atomToMols[atom] - excludeAtomToMols[atom] for atom in atomToMols
}
filtered_by_mols = pobj.filter_data_by_molecules(pobj.algoToDelta, atomToMols)
filteredData = pobj.filter_by_presence_of_experimental(filtered_by_mols)

# filtered_atomToMols = {atom: list(mols) for atom, mols in filteredData['mom'].items()}
pobj.calculate_errors(filteredData)
pobj.calculate_series_statistics()
pobj.calculate_overall_statistics()


eobj = extrapolation.WholeDataset(include_pcX=True, include_ccX=True)
eobj._create_extrapolation_schemes()
eobj.debug = False
# eobj.extrapolate_all_delta_data(filteredData, schemes=eobj.schemeDict["HF"])
eobj.extrapolate_all_delta_data(filteredData, schemes=eobj.schemeDict["CCSD"])
eobj.extrapolate_all_delta_data(filteredData, schemes=eobj.schemeDict["MP2"])
eobj.extrapolate_all_delta_data(filteredData, schemes=eobj.schemeDict["MP2_pure"])

eobj.calculate_errors(pobj.molToExper)
eobj.calculate_series_statistics()
eobj.calculate_overall_statistics()

# --------------- CREATE TABLES ---------------

table1 = create_tables.SingleZetaResults(
    atomToData=filteredData["mom"],
    molToExper=pobj.molToExper,
    save_folder=TABLE_PATH / "single-zeta",
)
table1.all_results()

table2 = create_tables.MethodSummary(
    pobj.algoToStats["mom"],
    pobj.algoToAtomStats["mom"],
    save_folder=TABLE_PATH / "method-summaries",
    show_sample_size=True,
    isPublication=True,
)
table2.all_results()

table3 = create_tables.ExtrapSchemeSummary(
    eobj.algoToStats["mom"],
    eobj.algoToAtomStats["mom"],
    save_folder=TABLE_PATH / "extrap-summaries",
    show_sample_size=False,
    isPublication=True,
)

eobj._schemeIterKeys = ["HF", "CCSD", "MP2"]
eobj._schemeIterKeys = ["CCSD"]
table3.prefix = "cc_"
table3.results_for_schemes(scheme_factory=eobj.scheme_generator)
eobj._schemeIterKeys = ["MP2_pure"]
table3.prefix = "mp_pure_"
table3.results_for_schemes(scheme_factory=eobj.scheme_generator)
eobj._schemeIterKeys = ["MP2"]
table3.prefix = "mp_comp_"
table3.results_for_schemes(scheme_factory=eobj.scheme_generator)

table4 = create_tables.UsedGeometries(
    geom_wb=str(DATA_PATH / "geometriesDB.xlsx"),
    save_folder=TABLE_PATH / "geometries",
)
relevantMols = {
    atom: set(molData.keys()) for atom, molData in filteredData["mom"].items()
}

table4.main(relevantMols)

# # --------------- CREATE FIGURES ---------------
create_figures._manual_delay()

figures_path = str(FIGURE_PATH)
create_figures.method_error_bars_general(pobj.algoToStats["mom"], figures_path)
create_figures.method_error_bars_series(
    atomToBasisStats=pobj.algoToAtomStats["mom"],
    bases="D T Q 5".split(),
    x_labels="cc-pCVDZ<br>/cc-pVDZ cc-pCVTZ<br>/cc-pVTZ cc-pCVQZ<br>/cc-pVQZ cc-pCV5Z<br>/cc-pV5Z".split(),
    save_path=figures_path,
    fname_suffix="_ccReg",
)
create_figures.method_error_bars_series(
    atomToBasisStats=pobj.algoToAtomStats["mom"],
    bases="pcX-1 pcX-2 pcX-3 pcX-4".split(),
    x_labels="pcX-1<br>/pc-1 pcX-2<br>/pc-2 pcX-3<br>/pc-3 pcX-4<br>/pc-4".split(),
    save_path=figures_path,
    fname_suffix="_pcX",
)
create_figures.method_error_bars_series(
    atomToBasisStats=pobj.algoToAtomStats["mom"],
    bases="ccX-DZ ccX-TZ ccX-QZ ccX-5Z".split(),
    x_labels="ccX-DZ/<br>cc-pVDZ ccX-TZ/<br>cc-pVTZ ccX-QZ/<br>cc-pVQZ ccX-5Z/<br>cc-pV5Z".split(),
    save_path=figures_path,
    fname_suffix="_ccX",
)

ccsdNames = "D-T-Q-CCSD D-T-Q-CCSD(T) T-Q-CCSD T-Q-CCSD(T)".split()
ccsdColors = ["#90caf9", "#64b5f6", "#2196f3", "#1976d2"]

create_figures.extrap_err_bars_dtstudy(
    eobj.algoToAtomStats["mom"], ccsdNames, ccsdColors, figures_path
)
create_figures.extrap_err_bars_dtstudy_general(
    eobj.algoToStats["mom"], ccsdNames, ccsdColors, figures_path
)


mp2Colors = ["#e0aaff", "#c77dff", "#9d4edd", "#7b2cbf"]
mp2_no5 = (
    "MP2[D T Q]+DifD | MP2[D T Q]+DifD(T) | MP2[T Q]+DifD | MP2[T Q]+DifD(T)".split(
        " | "
    )
)
mp2_w5 = "MP2[D T Q 5]+DifD | MP2[D T Q 5]+DifD(T) | MP2[T Q 5]+DifD | MP2[T Q 5]+DifD(T)".split(
    " | "
)
create_figures.extrap_err_bars_dtstudy(
    eobj.algoToAtomStats["mom"], mp2_no5, mp2Colors, figures_path, suffix="-mp2-no5"
)
create_figures.extrap_err_bars_dtstudy(
    eobj.algoToAtomStats["mom"], mp2_w5, mp2Colors, figures_path, suffix="-mp2-w5"
)

eobj.extrapolate_all_delta_data(filteredData, schemes=eobj.schemeDict["MP2_EXT"])
eobj.calculate_errors(pobj.molToExper)
eobj.calculate_series_statistics()
eobj.calculate_overall_statistics()

for include_pentuple in (True, False):
    create_figures.extrap_err_bars_summary(
        eobj.algoToAtomStats["mom"], figures_path, include_pentuple
    )
    create_figures.extrap_err_bars_atom_summary(
        eobj.algoToAtomStats["mom"], figures_path, include_pentuple
    )
    create_figures.benefit_of_extrapolation_over_special_basis_series(
        pobj.algoToAtomStats["mom"],
        eobj.algoToAtomStats["mom"],
        figures_path,
        include_pentuple,
    )
    create_figures.benefit_of_extrapolation_over_special_basis_overall(
        pobj.algoToAtomStats["mom"],
        eobj.algoToAtomStats["mom"],
        figures_path,
        include_pentuple,
    )

create_figures.small_basis_study_subplots(eobj.algoToAtomStats["mom"], figures_path)
create_figures.big_basis_study_subplots(eobj.algoToAtomStats["mom"], figures_path)
create_figures.extrap_err_bars_for_toc(eobj.algoToAtomStats["mom"], figures_path)


# ----- Correlation Figures

for include_pentuple in (True, False):
    create_figures.correlate_extrapolation_summary(
        eobj.algoToCBS["mom"], figures_path, include_pentuple
    )
# include5=True. j=0. R^2: 0.9996, RMSD: 0.037, MAE: 0.027
# include5=True. j=1. R^2: 0.9971, RMSD: 0.106, MAE: 0.074
# include5=True. j=2. R^2: 0.9997, RMSD: 0.033, MAE: 0.023
# include5=True. j=3. R^2: 0.9974, RMSD: 0.097, MAE: 0.063
# include5=False. j=0. R^2: 0.9997, RMSD: 0.036, MAE: 0.025
# include5=False. j=1. R^2: 0.9973, RMSD: 0.101, MAE: 0.072
# include5=False. j=2. R^2: 0.9997, RMSD: 0.031, MAE: 0.020
# include5=False. j=3. R^2: 0.9976, RMSD: 0.093, MAE: 0.061

ccColors = ["#90caf9", "#64b5f6", "#2196f3", "#1976d2"]
mpColors = ["#e0aaff", "#c77dff", "#9d4edd", "#7b2cbf"]
dStudy_xTitle = "Error from extrapolations<br>including double-zeta basis"
dStudy_yTitle = "Error from extrapolations<br>excluding double-zeta basis"
for dEffect, suffix, axrange in (
    (["D-T-Q-CCSD", "T-Q-CCSD"], "_dstudy-cc", [-0.79, 0.79]),
    (["D-T-Q-CCSD(T)", "T-Q-CCSD(T)"], "_dstudy-cc-triples", [-1.09, 0.59]),
    (["MP2[D T Q]+DifD", "MP2[T Q]+DifD"], "_dstudy-mp", [-0.79, 0.79]),
    (["MP2[D T Q]+DifD(T)", "MP2[T Q]+DifD(T)"], "_dstudy-mp-triples", [-1.09, 0.59]),
    (["T-Q-CCSD", "T-Q-CCSD(T)"], "_tstudy-cc", [-0.79, 0.79]),
    (["D-T-Q-CCSD", "D-T-Q-CCSD(T)"], "_tstudy-cc-double", [-1.09, 0.59]),
    (["MP2[T Q]+DifD", "MP2[T Q]+DifD(T)"], "_tstudy-mp", [-0.79, 0.79]),
    (["MP2[D T Q]+DifD", "MP2[D T Q]+DifD(T)"], "_tstudy-mp-double", [-1.09, 0.59]),
):
    if "mp" in suffix:
        colors = mpColors
    else:
        colors = ccColors
    create_figures.correlate_study(
        eobj.algoToAtomStats["mom"],
        dEffect,
        xtitle=dStudy_xTitle,
        ytitle=dStudy_yTitle,
        axrange=axrange,
        suffix=suffix,
        colors=colors,
        save_path=figures_path,
    )


# # --------------- ANALYZE TIMING ---------------

TimingObj = timing.Timer(
    calc_path=DATA_PATH / "calculations" / "mom",
    atom_mols=constants.filtered_atom_to_mols,
)
TimingObj.main(save_path=DATA_PATH)

# # ---------------- Reference Values Table ----------------
atomToExper: Dict[str, Dict[str, float]] = {}
for atom, molData in filteredData["mom"].items():
    for molecule in molData.keys():
        atomToExper.setdefault(atom, {}).setdefault(
            constants.FNAME_TO_MOLS[molecule]["latex"], pobj.molToExper[molecule]
        )


def create_reference_values_table(
    data: Dict[str, Dict[str, float]], columns: int = 3
) -> str:
    subcolumn_headers = ["Molecule", "CEBE (Exp)"] * columns
    header_format = " & ".join(subcolumn_headers) + " \\\\"

    latex_table = [
        "\\begin{table}",
        "\\centering",
        "\\begin{tabular}{" + "lc" * columns + "}",
        "\\toprule",
        header_format,
        "\\midrule",
    ]

    row_format = " & ".join(["{} & {:.2f}"] * columns) + " \\\\"

    items = [
        ("\ch{%s}" % molecule, cebe)
        for atom, molToCebe in data.items()
        for molecule, cebe in molToCebe.items()
    ]
    items_sorted = sorted(items, key=lambda x: x[-1])

    for row in itertools.zip_longest(
        *[iter(items_sorted)] * columns, fillvalue=("", 0)
    ):
        latex_table.append(row_format.format(*itertools.chain.from_iterable(row)))
    latex_table.extend(
        [
            "\\bottomrule",
            "\\end{tabular}",
            "\\caption{Energies of molecules by atom}",
            "\\end{table}",
        ]
    )

    return "\n".join(latex_table)


ref_table = create_reference_values_table(atomToExper, columns=3)
with open(DATA_PATH / "reference_values_table.tex", "w") as f:
    f.write(ref_table)
# # ---------------------------------------
