import sys, Analysis.constants as constants
sys.path.append('../../core_excitations') # this allows us to import files from parent directories
import PRIVATE
import tools
import analysis.modules.table 
import analysis.regression.evaluate
import analysis.extrapolation.extrapolations

class Tables:

    def __init__(self):
        self.latT = analysis.modules.table.LaTeXTable()
        self.evalObj = analysis.regression.evaluate.export()
        # self.extrapObj = analysis.extrapolation.extrapolations.export('O')


    def all_results(self):
        pass


if __name__ == "__main__":
    tables = Tables()



    