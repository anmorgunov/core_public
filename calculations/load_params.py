from typing import Dict, List


class ParamParser:
    """This is an object that's responsible for parsing launch_parameters.txt. Currently it's only used as a secondary tool in "generate_bash.py" """

    def __init__(self) -> None:
        """Initialize an object and read lines from "launch_parameters.txt" """
        with open("launch_parameters.txt", "r") as f:
            self.lines = f.readlines()
        self.launch_params: List[Dict[str, str]] = []

    def _parse(self) -> None:
        """The main workhorse for parsing launch_parameters.py"""
        endIndex = None
        for i, rline in enumerate(self.lines):
            line = rline.strip()
            if line == "END":
                continue

            if not line or line[0] == "#":
                continue

            if line == "START":
                for j in range(i, len(self.lines)):
                    # first, let's find how long is the block
                    if self.lines[j].rstrip() == "END":
                        endIndex = j
                        break
                if endIndex is None:
                    raise SyntaxError(
                        "launch_parameters.txt contains a block without END"
                    )
                param_dict = {}
                for paramIndex in range(i + 1, endIndex):
                    keyword, value = self.lines[paramIndex].rstrip().split(": ")
                    param_dict[keyword] = value

                if "atom" not in param_dict:
                    raise SyntaxError(
                        "You must specify an atom on which a core orbital will be ionized"
                    )
                if "mols" not in param_dict and "mols(loc)" not in param_dict:
                    raise SyntaxError(
                        "You must provide at least one molecule for calculations"
                    )
                self.launch_params.append(param_dict)

    def main(self) -> List[Dict[str, str]]:
        """As of now, this is essentially a proxy for getting self.launchParams

        Returns:
            list: a list of dictionaries that specify parameters for a job
        """
        self._parse()
        return self.launch_params


if __name__ == "__main__":
    parser = ParamParser()
    parser.main()
    print(parser.launch_params)
