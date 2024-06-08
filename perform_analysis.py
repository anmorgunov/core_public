import pprint
from pathlib import Path

from Analysis.Modules import Parser

pp = pprint.PrettyPrinter(indent=2)


def print(args):
    pp.pprint(args)


BASE_PATH = Path(__file__).resolve().parent
DATA_PATH = BASE_PATH / "Data"
CALCS_FOLDER = BASE_PATH / "calculations"

# --------------- PERFORM PARSING ---------------
pobj = Parser.ParsedData(
    experimental_wb=str(DATA_PATH / "experimental.xlsx"),
    calculations_wb=str(DATA_PATH / "CEBE_Data.xlsx"),
    algorithms=["mom"],
    calculations_folder=DATA_PATH / "calculations",
)
pobj.debug = True
pobj.process(save=False)
# pp.pprint(pobj.algoToData['mom']['o']['h2o']['D'])


# the calculations_wb file can become huge quite quickly, so if you want to check on results
# of certain molecules of interest, you can use the method below
mols = set("h2o co co2".split())  # specify molecules you're interested in
fname = "o-series"  # specify the file name to which this will be saved
# extracted_path = os.path.join(BASE_PATH, "Data", "extracts", f"{fname}.xlsx")
# pobj.extract_molecules(mols, extracted_path)

# --------------- PERFORM EXTRAPOLATION ---------------

# eobj = Extrapolation.WholeDataset()

# # First, let's filter out the molecules we want to use
# ATOM_TO_MOLS = constants.ATOM_TO_MOLS
# excludeAtomToMols = {
#     "C": set("c-o ch3-c-o2h".split()),
#     "N": set(
#         "mnh2pyridi-n-e onh2pyridi-n-e o-n-h2pyridine m-n-h2pyridine pfpyridine".split()
#     ),
#     "O": set("cf3co-o-h cf3c-o-oh ch3co-o-ch3 ch3c-o-och3 ".split()),
#     "F": set("cf3ocf3".split()),
# }
# atomToMols = {
#     atom.lower(): set(ATOM_TO_MOLS[atom]) - excludeAtomToMols[atom]
#     for atom in ATOM_TO_MOLS
# }
# # atomToMols = ATOM_TO_MOLS
# # atomToMols = constants.ATOM_PARTITIONED['big']
# # atomToMols['o'] = atomToMols['o']-constants.ATOM_PARTITIONED['small']['o']

# filteredData = pobj.filter_data_by_molecules(pobj.algoToData, atomToMols)
# filteredData = pobj.filter_by_presence_of_experimental(filteredData)

# # ---------------------------------------
# # atomToExper = {}
# # for atom, molData in filteredData["mom"].items():
# #     for molecule in molData.keys():
# #         atomToExper.setdefault(atom, {}).setdefault(constants.FNAME_TO_MOLS[molecule]["latex"], pobj.molToExper[molecule])

# # import itertools
# # def create_reference_values_table(data, columns=3):
# #     subcolumn_headers = ["Molecule", "CEBE (Exp)"] * columns
# #     header_format = " & ".join(subcolumn_headers) + " \\\\"

# #     latex_table = [
# #         "\\begin{table}",
# #         "\\centering",
# #         "\\begin{tabular}{" + "lc" * columns + "}",
# #         "\\toprule",
# #         header_format,
# #         "\\midrule"
# #     ]

# #     row_format = " & ".join(["{} & {:.2f}"] * columns) + " \\\\"

# #     items = [("\ch{%s}"%molecule, cebe) for atom, molToCebe in data.items() for molecule, cebe in molToCebe.items()]
# #     items_sorted = sorted(items, key=lambda x: x[-1])

# #     for row in itertools.zip_longest(*[iter(items_sorted)]*columns, fillvalue=("", 0)):
# #         latex_table.append(row_format.format(*itertools.chain.from_iterable(row)))
# #     latex_table.extend([
# #         "\\bottomrule",
# #         "\\end{tabular}",
# #         "\\caption{Energies of molecules by atom}",
# #         "\\end{table}"
# #     ])

# #     return "\n".join(latex_table)

# # ref_table = create_reference_values_table(atomToExper, columns=3)
# # with open(os.path.join(DATA_PATH, "reference_values_table.tex"), "w") as f:
# #     f.write(ref_table)
# # ---------------------------------------

# pobj.calculate_errors(filteredData)
# pobj.calculate_series_statistics()
# pobj.calculate_overall_statistics()


# # eobj._create_extrapolation_schemes()
# # eobj.debug = False
# # eobj.extrapolate_all_data(filteredData, schemes=eobj.schemeDict["HF"])
# # eobj.extrapolate_all_data(filteredData, schemes=eobj.schemeDict["CCSD"])
# # eobj.extrapolate_all_data(filteredData, schemes=eobj.schemeDict["MP2"])

# # eobj.calculate_errors(pobj.molToExper)
# # eobj.calculate_series_statistics()
# # eobj.calculate_overall_statistics()


# # --------------- CREATE TABLES ---------------

# table1 = CreateTables.SingleZetaResults(filteredData["mom"], pobj.molToExper)
# table1.all_results(save_folder=os.path.join(DATA_PATH, "paper-tables", "single-zeta"))

# table2 = CreateTables.MethodSummary(
#     pobj.algoToStats["mom"],
#     pobj.algoToAtomStats["mom"],
#     save_folder=os.path.join(DATA_PATH, "paper-tables", "method-summaries"),
#     show_sample_size=True,
#     isPublication=False,
# )
# table2.all_results()

# # # # breakpoint()

# # table3 = CreateTables.ExtrapSchemeSummary(
# #     eobj.algoToStats["mom"],
# #     eobj.algoToAtomStats["mom"],
# #     save_folder=os.path.join(DATA_PATH, "paper-tables", "extrap-summaries"),
# #     show_sample_size=True,
# #     isPublication=False,
# # )
# # # # print(eobj.algoToAtomStats['mom']['o'])
# # eobj._schemeIterKeys = ["HF", "CCSD", "MP2"]
# # table3.results_for_schemes(scheme_factory=eobj.scheme_generator)

# # table4 = CreateTables.UsedGeometries(
# #     geom_wb=os.path.join(DATA_PATH, "geometriesDB.xlsx"),
# #     save_folder=os.path.join(DATA_PATH, "paper-tables", "geometries"),
# # )
# # relevantMols = {
# #     atom: set(molData.keys()) for atom, molData in filteredData["mom"].items()
# # }

# # table4.main(relevantMols)

# # --------------- CREATE FIGURES ---------------
# # CreateFigures._manual_delay()

# # figures_path = os.path.join(DATA_PATH, "paper-figures")
# # # CreateFigures.method_error_bars_general(pobj.algoToStats["mom"], figures_path)
# # CreateFigures.method_error_bars_series(
# #     atomToBasisStats=pobj.algoToAtomStats["mom"],
# #     bases="D T Q 5".split(),
# #     save_path=figures_path,
# #     fname_suffix="_ccReg"
# # )
# # CreateFigures.method_error_bars_series(
# #     atomToBasisStats=pobj.algoToAtomStats["mom"],
# #     bases="pcX-1 pcX-2 pcX-3 pcX-4".split(),
# #     save_path=figures_path,
# #     fname_suffix="_pcX"
# # )
# # CreateFigures.method_error_bars_series(
# #     atomToBasisStats=pobj.algoToAtomStats["mom"],
# #     bases="ccX-DZ ccX-TZ ccX-QZ ccX-5Z".split(),
# #     save_path=figures_path,
# #     fname_suffix="_ccX"
# # )
# # ccsdNames = "D-T-Q-CCSD D-T-Q-CCSD(T) T-Q-CCSD T-Q-CCSD(T)".split()
# # ccsdColors = ["#90caf9", "#64b5f6", "#2196f3", "#1976d2"]

# # CreateFigures.extrap_err_bars_dtstudy(
# #     eobj.algoToAtomStats["mom"], ccsdNames, ccsdColors, figures_path
# # )
# # CreateFigures.extrap_err_bars_dtstudy_general(
# #     eobj.algoToStats["mom"], ccsdNames, ccsdColors, figures_path
# # )


# mp2Colors = ["#e0aaff", "#c77dff", "#9d4edd", "#7b2cbf"]
# mp2_no5 = (
#     "MP2[D T Q]+DifD | MP2[D T Q]+DifD(T) | MP2[T Q]+DifD | MP2[T Q]+DifD(T)".split(
#         " | "
#     )
# )
# mp2_w5 = "MP2[D T Q 5]+DifD | MP2[D T Q 5]+DifD(T) | MP2[T Q 5]+DifD | MP2[T Q 5]+DifD(T)".split(
#     " | "
# )
# CreateFigures.extrap_err_bars_dtstudy(
#     eobj.algoToAtomStats["mom"], mp2_no5, mp2Colors, figures_path, suffix="-mp2-no5"
# )
# CreateFigures.extrap_err_bars_dtstudy(
#     eobj.algoToAtomStats["mom"], mp2_w5, mp2Colors, figures_path, suffix="-mp2-w5"
# )

# eobj.extrapolate_all_data(filteredData, schemes=eobj.schemeDict["MP2_EXT"])
# eobj.calculate_errors(pobj.molToExper)
# eobj.calculate_series_statistics()
# eobj.calculate_overall_statistics()

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

# ccColors = ["#90caf9", "#64b5f6", "#2196f3", "#1976d2"]
# mpColors = ["#e0aaff", "#c77dff", "#9d4edd", "#7b2cbf"]
# dStudy_xTitle = "Error from extrapolations<br>including double-zeta basis"
# dStudy_yTitle = "Error from extrapolations<br>excluding double-zeta basis"
# for dEffect, suffix, axrange in (
#     (["D-T-Q-CCSD", "T-Q-CCSD"], "_dstudy-cc", [-0.79, 0.79]),
#     (["D-T-Q-CCSD(T)", "T-Q-CCSD(T)"], "_dstudy-cc-triples", [-1.09, 0.59]),
#     (["MP2[D T Q]+DifD", "MP2[T Q]+DifD"], "_dstudy-mp", [-0.79, 0.79]),
#     (["MP2[D T Q]+DifD(T)", "MP2[T Q]+DifD(T)"], "_dstudy-mp-triples", [-1.09, 0.59]),
#     (["T-Q-CCSD", "T-Q-CCSD(T)"], "_tstudy-cc", [-0.79, 0.79]),
#     (["D-T-Q-CCSD", "D-T-Q-CCSD(T)"], "_tstudy-cc-double", [-1.09, 0.59]),
#     (["MP2[T Q]+DifD", "MP2[T Q]+DifD(T)"], "_tstudy-mp", [-0.79, 0.79]),
#     (["MP2[D T Q]+DifD", "MP2[D T Q]+DifD(T)"], "_tstudy-mp-double", [-1.09, 0.59]),
# ):
#     if "mp" in suffix:
#         colors = mpColors
#     else:
#         colors = ccColors
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

# TimingObj = Timing.Timer(
#     calc_path=DATA_PATH / "calculations" / "mom", atoms={"o"}#{"c", "n", "o", "f"}
# )
# TimingObj.main(save_path=DATA_PATH)
