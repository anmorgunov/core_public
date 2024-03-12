
import Analysis

from Analysis.Modules import Parser

import os
import pprint
pp = pprint.PrettyPrinter(indent=2)

BASE_PATH = os.path.dirname(__file__)

# --------------- PERFORM PARSING ---------------
experimental_wb = os.path.join(BASE_PATH, "Data", "experimental.xlsx")
calculations_folder = os.path.join(BASE_PATH, "Data", "calculations")
calculations_wb = os.path.join(BASE_PATH, "Data", "CEBE_Data.xlsx")
methods = ['mom']
pobj = Parser.ParsedData(experimental_wb, calculations_folder, methods, calculations_wb)
# pobj.debug = True
pobj.main(save=True)
# pp.pprint(pobj.methodToData['mom']['o']['h2o']['D'])


# the calculations_wb file can become huge quite quickly, so if you want to check on results
# of certain molecules of interest, you can use the method below
mols = set("h2o co co2".split()) # specify molecules you're interested in
fname = 'o-series' # specify the file name to which this will be saved
extracted_path = os.path.join(BASE_PATH, 'Data', 'extracts', f'{fname}.xlsx')
pobj.extract_molecules(mols, extracted_path)


