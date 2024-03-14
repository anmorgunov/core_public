from Analysis import constants
import Analysis

from Analysis.Modules import Parser
from Analysis.Modules import Extrapolation
from Analysis import CreateTables
from Analysis import CreateFigures
from Analysis import Timing

import os
import pprint

pp = pprint.PrettyPrinter(indent=2)


def print(args):
    pp.pprint(args)


BASE_PATH = os.path.dirname(__file__)
DATA_PATH = os.path.join(BASE_PATH, "Data")

# --------------- PERFORM PARSING ---------------
experimental_wb = os.path.join(DATA_PATH, "experimental.xlsx")
calculations_folder = os.path.join(DATA_PATH, "calculations")
calculations_wb = os.path.join(DATA_PATH, "CEBE_Data.xlsx")
algorithms = ["mom"]
pobj = Parser.ParsedData(
    experimental_wb, calculations_folder, algorithms, calculations_wb
)
# pobj.debug = True
pobj.main(save=False)
# pp.pprint(pobj.algoToData['mom']['o']['h2o']['D'])


# the calculations_wb file can become huge quite quickly, so if you want to check on results
# of certain molecules of interest, you can use the method below
mols = set("h2o co co2".split())  # specify molecules you're interested in
fname = "o-series"  # specify the file name to which this will be saved
extracted_path = os.path.join(BASE_PATH, "Data", "extracts", f"{fname}.xlsx")
# pobj.extract_molecules(mols, extracted_path)

# --------------- PERFORM EXTRAPOLATION ---------------

eobj = Extrapolation.WholeDataset()
# eobj.parse_scheme('T-Q-CCSD')
# eobj.parse_scheme('T-Q-CCSD(T)')
# eobj.parse_scheme('MP2[D T Q 5]')
# eobj.parse_scheme('MP2[D T Q]+DifD')
# eobj.parse_scheme('MP2[D T Q 5]+DifD(T)')

# r1 = eobj.extrapolate_molecule_given_scheme(pobj.algoToData['mom']['o']['h2o'], 'T-Q-CCSD')

# r2 = eobj.extrapolate_molecule_given_scheme(pobj.algoToData['mom']['o']['h2o'], 'T-Q-CCSD(T)')

# r3 = eobj.extrapolate_molecule_given_scheme(pobj.algoToData['mom']['o']['h2o'], 'MP2[T Q 5]')

# r4 = eobj.extrapolate_molecule_given_scheme(pobj.algoToData['mom']['o']['h2o'], 'MP2[T Q 5]+DifD')

# print(f"{r1=}\n{r2=}\n{r3=}\n{r4=}")

# First, let's filter out the molecules we want to use
atomToMols = {
    atom.lower(): set(molStr.split()) for atom, molStr in constants.ATOM_TO_MOLS.items()
}

filteredData = pobj.filter_data_by_molecules(pobj.algoToData, atomToMols)
filteredData = pobj.filter_by_presence_of_experimental(filteredData)
# pp.pprint(atomToMols)

# eobj.extrapolate_all_data(
#     filteredData, schemes=["T-Q-CCSD", "T-Q-5-HF", "MP2[T Q]+DifD(T)"]
# )
eobj._create_extrapolation_schemes()
eobj.debug = False
# eobj.parse_scheme(scheme='MP2[T Q 5]+DifSTO-3G(T)')
eobj.extrapolate_all_data(filteredData, schemes=eobj.schemeDict["HF"])
eobj.extrapolate_all_data(filteredData, schemes=eobj.schemeDict["CCSD"])
eobj.extrapolate_all_data(filteredData, schemes=eobj.schemeDict["MP2"])
# pp.pprint(pobj.molToExper)
eobj.calculate_errors(pobj.molToExper)
# print(eobj.algorithmToError['mom']['o']['h2o'])
eobj.calculate_series_statistics()
# print(eobj.algoToAtomStats['mom']['o']['T-Q-CCSD'])
eobj.calculate_overall_statistics()
# print(eobj.algoToStats['mom']['T-Q-CCSD'])
# print(eobj.algoToStats['mom']['MP2[D T Q 5]'])
# print(eobj.algoToStats['mom']['MP2[D T Q 5]+DifD(T)'])
# print(eobj.smallBasisException)

pobj.calculate_errors(filteredData)
pobj.calculate_series_statistics()
pobj.calculate_overall_statistics()

# --------------- CREATE TABLES ---------------

table1 = CreateTables.SingleZetaResults(filteredData["mom"], pobj.molToExper)
table1.all_results(save_folder=os.path.join(DATA_PATH, "paper-tables", "single-zeta"))

table2 = CreateTables.MethodSummary(
    pobj.algoToStats["mom"],
    pobj.algoToAtomStats["mom"],
    save_folder=os.path.join(DATA_PATH, "paper-tables", "method-summaries"),
    show_sample_size=False,
    isPublication=True,
)
table2.all_results()

table3 = CreateTables.ExtrapSchemeSummary(
    eobj.algoToStats["mom"],
    eobj.algoToAtomStats["mom"],
    save_folder=os.path.join(DATA_PATH, "paper-tables", "extrap-summaries"),
    show_sample_size=False,
    isPublication=True,
)
# print(eobj.algoToAtomStats['mom']['o'])
eobj._schemeIterKeys = ["HF", "CCSD", "MP2"]
table3.results_for_schemes(scheme_factory=eobj.scheme_generator)

table4 = CreateTables.UsedGeometries(
    geom_wb=os.path.join(DATA_PATH, "geometriesDB.xlsx"),
    save_folder=os.path.join(DATA_PATH, "paper-tables", "geometries"),
)
relevantMols = {
    atom: set(molData.keys()) for atom, molData in filteredData["mom"].items()
}

table4.main(relevantMols)

# --------------- CREATE FIGURES ---------------
# CreateFigures._manual_delay()

figures_path = os.path.join(DATA_PATH, "paper-figures")
# CreateFigures.method_error_bars_general(pobj.algoToStats["mom"], figures_path)
# CreateFigures.method_error_bars_series(pobj.algoToAtomStats["mom"], figures_path)

ccsdNames = "D-T-Q-CCSD D-T-Q-CCSD(T) T-Q-CCSD T-Q-CCSD(T)".split()
ccsdColors = ["#90caf9", "#64b5f6", "#2196f3", "#1976d2"]

# CreateFigures.extrap_err_bars_dtstudy(
#     eobj.algoToAtomStats["mom"], ccsdNames, ccsdColors, figures_path
# )
# CreateFigures.extrap_err_bars_dtstudy_general(
#     eobj.algoToStats["mom"], ccsdNames, ccsdColors, figures_path
# )


mp2Colors = ["#e0aaff", "#c77dff", "#9d4edd", "#7b2cbf"]
mp2_no5 = (
    "MP2[D T Q]+DifD | MP2[D T Q]+DifD(T) | MP2[T Q]+DifD | MP2[T Q]+DifD(T)".split(
        " | "
    )
)
mp2_w5 = "MP2[D T Q 5]+DifD | MP2[D T Q 5]+DifD(T) | MP2[T Q 5]+DifD | MP2[T Q 5]+DifD(T)".split(
    " | "
)
# CreateFigures.extrap_err_bars_dtstudy(
#     eobj.algoToAtomStats["mom"], mp2_no5, mp2Colors, figures_path, suffix="-mp2-no5"
# )
# CreateFigures.extrap_err_bars_dtstudy(
#     eobj.algoToAtomStats["mom"], mp2_w5, mp2Colors, figures_path, suffix="-mp2-w5"
# )

eobj.extrapolate_all_data(filteredData, schemes=eobj.schemeDict["MP2_EXT"])
eobj.calculate_errors(pobj.molToExper)
eobj.calculate_series_statistics()
eobj.calculate_overall_statistics()

# for include_pentuple in (True, False):
#     CreateFigures.extrap_err_bars_summary(eobj.algoToAtomStats["mom"], figures_path, include_pentuple)

# CreateFigures.small_basis_study_subplots(eobj.algoToAtomStats["mom"], figures_path)
# CreateFigures.big_basis_study_subplots(eobj.algoToAtomStats["mom"], figures_path)
# CreateFigures.extrap_err_bars_for_toc(eobj.algoToAtomStats["mom"], figures_path)


# ----- Correlation Figures

# for include_pentuple in (True, False):
#     CreateFigures.correlate_extrapolation_summary(
#         eobj.algoToCBS["mom"], figures_path, include_pentuple
#     )

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
    # CreateFigures.correlate_study(
    #     eobj.algoToAtomStats["mom"],
    #     dEffect,
    #     xtitle=dStudy_xTitle,
    #     ytitle=dStudy_yTitle,
    #     axrange=axrange,
    #     suffix=suffix,
    #     colors=colors,
    #     save_path=figures_path,
    # )


# --------------- ANALYZE TIMING ---------------

TimingObj = Timing.Timer(
    calc_path=os.path.join(DATA_PATH, "calculations", "mom"), atoms={"c", "n", "o", "f"}
)
TimingObj.main(save_path=DATA_PATH)
