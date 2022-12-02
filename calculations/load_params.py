# import tools 
# import PRIVATE 

class ParamParser:
    """This is an object that's responsible for parsing launch_parameters.txt. Currently it's only used as a secondary tool in "generate_bash.py"
    """

    def __init__(self):
        """Initialize an object and read lines from "launch_parameters.txt"
        """
        # with open(tools.join_path(PRIVATE.BASE_PATH + ["calculations", 'launch_parameters.txt']), 'r') as f:
        with open('launch_parameters.txt', 'r') as f:
            self.lines = f.readlines()
        self.atomToData = {}
        self.launchParams = []

        
    def _parse(self):
        """The main workhorse for parsing launch_parameters.py
        """        
        endIndex = None
        for i, rline in enumerate(self.lines): #for each line in launch_parameters.txt
            line = rline.strip() # remove newline character
            if line == 'END': # if we're at the end of a block, just continue iterating
                continue

            if not line or line[0] == '#': #if an empty line or line starts with a # we ignore it
                continue

            if line == 'START': #if we found a line that says "START", we start parsing
                for j in range(i, len(self.lines)): #first, let's find how long is the block
                    if self.lines[j].rstrip() == 'END':
                        endIndex = j
                        break
                assert endIndex is not None, "Launch_parameters.txt contains a block without END" # the block must end somewhere

                launchParamSet = {}
                for paramIndex in range(i+1, endIndex): # parse each line between START and END 
                    keyword, value = self.lines[paramIndex].rstrip().split(': ')
                    launchParamSet[keyword] = value
                
                assert 'atom' in launchParamSet, 'You must specify an atom on which a core orbital will be ionized'
                assert 'mols' in launchParamSet or 'mols(loc)' in launchParamSet, 'You must provide at least one molecule for calculations'
                self.launchParams.append(launchParamSet)

    def main(self):
        """As of now, this is essentially a proxy for getting self.launchParams

        Returns:
            list: a list of dictionaries that specify parameters for a job 
        """
        self._parse() 
        # return self.atomToData
        return self.launchParams

if __name__ == "__main__":
    sg = ParamParser()
    # sg.main()
    # print(sg.atomToData)
    # print(sg.launchParams)